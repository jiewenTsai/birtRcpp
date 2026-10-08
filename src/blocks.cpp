// [[Rcpp::depends(RcppArmadillo)]]
#include "blocks.h"
#include "rng.h"

Ald ald_constants(double q) { return Ald{(1.0 - 2.0 * q) / (q * (1.0 - q)), 2.0 / (q * (1.0 - q))}; }

Prior make_prior(const Rcpp::NumericVector& v) {
  if (v.size() != 18) Rcpp::stop("prior vector must have 18 elements");
  return Prior{v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15], v[16], v[17]};
}
Prior julia_prior() {
  return Prior{1.0, 1.0, 0.0, 1.0, -4.0, 4.0, NA_REAL, NA_REAL, 0.0, 0.0, 0.001, 0.001, 1.0, 1.0, 0.0, 0.001, 0.001, 1.0};
}

double draw_var(double alpha, double beta, double fam, double p1, double p2, double current) {
  if (fam == 0.0) return rinvgamma(p1 + alpha, p2 + beta);
  // half-t(df = p1, scale = p2) on sqrt(var) (Huang & Wand, 2013): var | u ~ IG(df / 2, df / u),
  // u ~ IG(1 / 2, 1 / scale^2). Gibbs steps u | var ~ IG((df + 1) / 2, df / var + 1 / scale^2), then
  // var | u, . ~ IG(df / 2 + alpha, df / u + beta); u is redrawn from the current var, so it is not stored.
  const double u = rinvgamma(0.5 * (p1 + 1.0), p1 / current + 1.0 / (p2 * p2));
  return rinvgamma(0.5 * p1 + alpha, p1 / u + beta);
}

// omega_jk | . ~ PG(1, a_k (theta_j - b_k))
arma::mat draw_omega(const arma::vec& a, const arma::vec& b, const arma::vec& theta) {
  arma::mat omega(theta.n_elem, a.n_elem);
  for (arma::uword k = 0; k < a.n_elem; ++k)
    for (arma::uword j = 0; j < theta.n_elem; ++j) omega(j, k) = rpg1(a(k) * (theta(j) - b(k)));
  return omega;
}

// a_k | . : the Polya-Gamma likelihood is normal in a, exp(nl a - pl a^2 / 2) with pl = sum_j omega (theta - b)^2,
// nl = sum_j kappa (theta - b). Prior N+(a_mu, a_sd^2): N+(num / prec, 1 / prec), prec = 1 / a_sd^2 + pl,
// num = a_mu / a_sd^2 + nl.
arma::vec draw_a(const arma::mat& kappa, const arma::mat& omega, const arma::vec& theta, const arma::vec& b, const Prior& P) {
  arma::vec a(b.n_elem);
  for (arma::uword k = 0; k < b.n_elem; ++k) {
    arma::vec d = theta - b(k);
    double prec = 1.0 / (P.a_sd * P.a_sd) + arma::dot(omega.col(k), d % d);
    double num = P.a_mu / (P.a_sd * P.a_sd) + arma::dot(kappa.col(k), d);
    a(k) = rtnorm(num / prec, std::sqrt(1.0 / prec), 0.0, R_PosInf);
  }
  return a;
}

// b_k | . ~ N(num / prec, 1 / prec) on [b_lo, b_hi]: prec = 1 / b_sd^2 + a^2 sum_j omega,
//                                                    num = b_mu / b_sd^2 + a sum_j (a omega theta - kappa)
arma::vec draw_b(const arma::mat& kappa, const arma::mat& omega, const arma::vec& theta, const arma::vec& a, const Prior& P) {
  arma::vec b(a.n_elem);
  for (arma::uword k = 0; k < a.n_elem; ++k) {
    double prec = 1.0 / (P.b_sd * P.b_sd) + a(k) * a(k) * arma::accu(omega.col(k));
    double num = P.b_mu / (P.b_sd * P.b_sd) + a(k) * arma::accu(a(k) * omega.col(k) % theta - kappa.col(k));
    b(k) = rtnorm(num / prec, std::sqrt(1.0 / prec), P.b_lo, P.b_hi);
  }
  return b;
}

// theta_j: prec += sum_k a^2 omega, num += sum_k a (kappa + a b omega)
void irt_part(const arma::mat& kappa, const arma::mat& omega, const arma::vec& a, const arma::vec& b,
              arma::vec& prec, arma::vec& num) {
  prec = omega * (a % a);
  num = (kappa + omega.each_row() % (a % b).t()) * a;
}

// lambda_k | . ~ N on (lo_l, Inf): prec = 1/sd_l^2 + sum_j W, num = mu_l/sd_l^2 + sum_j W (L - m - k1 nu)
arma::vec draw_lambda(const arma::mat& L, const arma::mat& m, const arma::mat& k1nu, const arma::mat& W,
                      double mu_l, double sd_l, double lo_l) {
  arma::vec lambda(L.n_cols);
  for (arma::uword k = 0; k < L.n_cols; ++k) {
    double prec = 1.0 / (sd_l * sd_l) + arma::accu(W.col(k));
    double num = mu_l / (sd_l * sd_l) + arma::dot(W.col(k), L.col(k) - m.col(k) - k1nu.col(k));
    lambda(k) = rtnorm(num / prec, std::sqrt(1.0 / prec), lo_l, R_PosInf);
  }
  return lambda;
}

// sigma2_k | . : likelihood sigma2^(-N/2) exp(-sum_j r^2 / (2 sigma2)), r = L - lambda - m (draw_var)
arma::vec draw_sigma2_normal(const arma::mat& L, const arma::mat& m, const arma::vec& lambda, const Prior& P, const arma::vec& current) {
  arma::vec s(L.n_cols);
  for (arma::uword k = 0; k < L.n_cols; ++k) {
    arma::vec r = L.col(k) - lambda(k) - m.col(k);
    s(k) = draw_var(0.5 * L.n_rows, 0.5 * arma::dot(r, r), P.s2_fam, P.s2_p1, P.s2_p2, current(k));
  }
  return s;
}

// ALD scale (Kozumi & Kobayashi, 2011): likelihood s^(-3N/2) exp(-(sum (r - k1 nu)^2 / (2 k2 nu) + sum nu) / s)
arma::vec draw_sigma2_ald(const arma::mat& L, const arma::mat& m, const arma::vec& lambda,
                          const arma::mat& nu, const Ald& c, const Prior& P, const arma::vec& current) {
  arma::vec s(L.n_cols);
  for (arma::uword k = 0; k < L.n_cols; ++k) {
    arma::vec r = L.col(k) - lambda(k) - m.col(k) - c.k1 * nu.col(k);
    const double A = arma::accu(r % r / (2.0 * c.k2 * nu.col(k))), Bn = arma::accu(nu.col(k));
    s(k) = P.s2_fam == 0.0 ? rinvgamma(P.s2_p1 + 1.5 * L.n_rows, P.s2_p2 + A + Bn)
                           : draw_var(1.5 * L.n_rows, A + Bn, P.s2_fam, P.s2_p1, P.s2_p2, current(k));
  }
  return s;
}

double draw_location_shift(arma::vec& lambda, arma::vec& zeta, const arma::vec& mz, const arma::vec& wz,
                           double mu_l, double sd_l, double lo_l) {
  double prec = lambda.n_elem / (sd_l * sd_l) + arma::accu(wz);
  double num = arma::accu(mu_l - lambda) / (sd_l * sd_l) + arma::dot(wz, mz - zeta);
  double c = num / prec + R::norm_rand() / std::sqrt(prec);
  if (arma::min(lambda) + c <= lo_l) return 0.0;     // lambda truncated at lo_l: keep the state
  lambda += c; zeta += c;
  return c;
}

// 1 / nu ~ inverse Gaussian(sqrt(psi / chi), psi)
double draw_nu(double r, double s, const Ald& c) {
  double psi = (c.k1 * c.k1 + 2.0 * c.k2) / (c.k2 * s), chi = r * r / (c.k2 * s);
  double mu = std::min(std::max(std::sqrt(psi / chi), 1e-10), 1e10);
  return std::min(std::max(1.0 / rinvgauss(mu, psi), 1e-10), 1e10);
}

arma::vec draw_theta_collapsed(const arma::vec& prec, const arma::vec& num, const arma::mat& W, const arma::mat& y,
                               const arma::vec& rho, const arma::vec& m0, const arma::vec& g, const arma::vec& tau) {
  // joint precision Q and linear term h of (theta_j, zeta_j); theta's marginal is the Schur complement
  arma::vec Wy = arma::sum(W % y, 1);
  arma::vec Qzz = arma::sum(W, 1) + 1.0 / tau, Qtz = W * rho - g / tau;
  arma::vec Qtt = prec + W * (rho % rho) + g % g / tau;
  arma::vec ht = num + (W % y) * rho - g % m0 / tau, hz = Wy + m0 / tau;
  arma::vec p = Qtt - Qtz % Qtz / Qzz, h = ht - Qtz % hz / Qzz;
  arma::vec theta(prec.n_elem);
  for (arma::uword j = 0; j < theta.n_elem; ++j) theta(j) = h(j) / p(j) + R::norm_rand() / std::sqrt(p(j));
  return theta;
}

void draw_rt_items_collapsed(const arma::mat& L, const arma::mat& k1nu, const arma::mat& W, const arma::vec& theta,
                             const arma::vec& mz, const arma::vec& tau, bool with_rho, double mu_l, double sd_l,
                             double rho_sd, arma::vec& lambda, arma::vec& rho, double lo_l) {
  const arma::uword K = L.n_cols;
  // Cov(u_j)^-1 = D_j - g_j d_j d_j' (Sherman-Morrison), d_j = W_j, g_j = 1 / (1 / tau_j + sum_k d_jk)
  arma::vec gj = 1.0 / (1.0 / tau + arma::sum(W, 1));
  arma::mat Z = L - k1nu; Z.each_col() += mz;
  arma::mat WZ = W % Z;
  arma::vec wz = arma::sum(WZ, 1);
  auto quad = [&](const arma::vec& f) {               // sum_j f_j A_j
    arma::mat Wf = W; Wf.each_col() %= gj % f;
    return arma::mat(arma::diagmat(W.t() * f) - W.t() * Wf);
  };
  auto lin = [&](const arma::vec& f) { return arma::vec(WZ.t() * f - W.t() * (gj % f % wz)); };
  const arma::vec one(L.n_rows, arma::fill::ones);
  const arma::uword D = with_rho ? 2 * K : K;
  arma::mat Q(D, D, arma::fill::zeros); arma::vec h(D);
  Q.submat(0, 0, K - 1, K - 1) = quad(one) + arma::eye(K, K) / (sd_l * sd_l);
  h.head(K) = lin(one) + mu_l / (sd_l * sd_l);
  if (with_rho) {                                     // design [I, -theta_j I]; rho_k ~ N(0, rho_sd^2)
    arma::mat S1 = quad(theta);
    Q.submat(0, K, K - 1, D - 1) = -S1; Q.submat(K, 0, D - 1, K - 1) = -S1;
    Q.submat(K, K, D - 1, D - 1) = quad(theta % theta) + arma::eye(K, K) / (rho_sd * rho_sd);
    h.tail(K) = -lin(theta);
  }
  // exact draw (rejection on lambda > lo_l when truncated); if that keeps failing, one Gibbs sweep through the block
  for (int tries = 0; tries < 50; ++tries) {
    arma::vec x = rmvn_prec(Q, h);
    if (x.head(K).min() > lo_l) { lambda = x.head(K); if (with_rho) rho = x.tail(K); return; }
  }
  arma::vec x = with_rho ? arma::join_cols(lambda, rho) : lambda;
  for (arma::uword i = 0; i < D; ++i) {
    double m = (h(i) - arma::dot(Q.col(i), x) + Q(i, i) * x(i)) / Q(i, i);
    x(i) = rtnorm(m, 1.0 / std::sqrt(Q(i, i)), i < K ? lo_l : R_NegInf, R_PosInf);
  }
  lambda = x.head(K); if (with_rho) rho = x.tail(K);
}

// One Metropolis-Hastings step for the two parameters x = (x0, x1) of an item given theta, with the latent
// responses (Polya-Gamma or Albert-Chib variables) integrated out. eta_j = x0 (theta_j - x1) (form 0),
// x0 theta_j - x1 (form 1) or x0 theta_j + x1 (form 2); logit (link 0) or probit (link 1); missing responses
// (NaN) contribute nothing. The proposal is a Fisher-scoring step, N(x + G(x)^-1 grad f(x), G(x)^-1), with f
// the log posterior and G the Fisher information plus the prior's metric: a position-dependent proposal, so
// the reverse proposal density enters the acceptance ratio. Proposals outside the prior's support are rejected.
static void item_terms(const arma::vec& y, const arma::uvec& obs, const arma::vec& theta, const arma::vec& x, int link, int form,
                       const ItemPrior& prior, double& f, arma::vec& g, arma::mat& G) {
  const double a = x(0), c = x(1), dc = form == 0 ? -a : (form == 1 ? -1.0 : 1.0);
  f = 0.0; g.zeros(2); G.zeros(2, 2);
  for (arma::uword j : obs) {
    const double eta = form == 0 ? a * (theta(j) - c) : a * theta(j) + dc * c, da = form == 0 ? theta(j) - c : theta(j);
    double ll, r, w;
    if (link == 0) {
      const double p = 1.0 / (1.0 + std::exp(-eta));
      ll = y(j) * eta - (eta > 0 ? eta + std::log1p(std::exp(-eta)) : std::log1p(std::exp(eta)));
      r = y(j) - p; w = p * (1.0 - p);
    } else {
      const double l1 = R::pnorm(eta, 0.0, 1.0, 1, 1), l0 = R::pnorm(eta, 0.0, 1.0, 0, 1), ld = R::dnorm(eta, 0.0, 1.0, 1);
      ll = y(j) > 0.5 ? l1 : l0;
      r = y(j) > 0.5 ? std::exp(ld - l1) : -std::exp(ld - l0);
      w = std::exp(2.0 * ld - l1 - l0);
    }
    f += ll; g(0) += r * da; g(1) += r * dc;
    G(0, 0) += w * da * da; G(0, 1) += w * da * dc; G(1, 1) += w * dc * dc;
  }
  G(1, 0) = G(0, 1);
  prior.terms(x, f, g, G);
}

// log density of N(x; mu, G^-1) up to the constant
static double ldnorm_prec(const arma::vec& x, const arma::vec& mu, const arma::mat& G) {
  arma::vec d = x - mu;
  return 0.5 * std::log(arma::det(G)) - 0.5 * arma::dot(d, G * d);
}

bool item_mh_step(const arma::vec& y, const arma::vec& theta, arma::vec& x, int link, int form, const ItemPrior& prior) {
  const arma::uvec obs = arma::find_finite(y);
  double f0, f1; arma::vec g0, g1; arma::mat G0, G1;
  item_terms(y, obs, theta, x, link, form, prior, f0, g0, G0);
  arma::vec mu0 = x + arma::solve(G0, g0);
  arma::mat U = arma::chol(G0);                        // G0 = U' U, so U^-1 z ~ N(0, G0^-1)
  arma::vec z = {R::norm_rand(), R::norm_rand()};
  arma::vec xp = mu0 + arma::solve(arma::trimatu(U), z);
  if (!prior.support(xp)) return false;
  item_terms(y, obs, theta, xp, link, form, prior, f1, g1, G1);
  arma::vec mu1 = xp + arma::solve(G1, g1);
  if (std::log(R::unif_rand()) < f1 - f0 + ldnorm_prec(x, mu1, G1) - ldnorm_prec(xp, mu0, G0)) { x = xp; return true; }
  return false;
}

// The samplers' joint step works in the intercept parametrisation eta = a theta + d (d = -a b, form 2): the
// logistic likelihood is log-concave in (a, d), so the posterior is close to an ellipse there even for very
// easy or very hard items, where it is a curved ridge in (a, b). The priors stay on (a, b) (a ~ N+, b
// truncated normal, see Prior); with b = -d / a the density of (a, d) carries the Jacobian 1 / a.
int draw_ab_joint(const arma::mat& Y, const arma::vec& theta, arma::vec& a, arma::vec& b, const Prior& P) {
  const double va = P.a_sd * P.a_sd, vb = P.b_sd * P.b_sd;
  const ItemPrior prior{
    [&P](const arma::vec& x) { const double b = -x(1) / x(0); return x(0) > 0.0 && b >= P.b_lo && b <= P.b_hi; },
    [&P, va, vb](const arma::vec& x, double& f, arma::vec& g, arma::mat& G) {
      const double a = x(0), d = x(1), b = -d / a, eb = (b - P.b_mu) / vb;
      f += -0.5 * (a - P.a_mu) * (a - P.a_mu) / va - 0.5 * (b - P.b_mu) * eb - std::log(a);
      g(0) += -(a - P.a_mu) / va - eb * d / (a * a) - 1.0 / a; g(1) += eb / a;     // db/da = d / a^2, db/dd = -1 / a
      G(0, 0) += 1.0 / va; G(1, 1) += 1.0 / (vb * a * a);
    }};
  int accepted = 0;
  for (arma::uword k = 0; k < a.n_elem; ++k) {
    arma::vec x = {a(k), -a(k) * b(k)};                 // (a, d)
    if (item_mh_step(Y.col(k), theta, x, 0, 2, prior)) { a(k) = x(0); b(k) = -x(1) / x(0); ++accepted; }
  }
  return accepted;
}

// exported for testing: draws of (a, b) of one item given a fixed theta, by the Polya-Gamma Gibbs
// steps (joint = false) or the joint Metropolis-Hastings step (joint = true)
// [[Rcpp::export(.item_chain)]]
arma::mat item_chain(const arma::vec& y, const arma::vec& theta, int n_iter, bool joint) {
  arma::mat Y(y), kappa = Y - 0.5, out(n_iter, 2);
  arma::vec a = {1.0}, b = {0.0};
  for (int it = 0; it < n_iter; ++it) {
    const Prior P = julia_prior();
    if (joint) draw_ab_joint(Y, theta, a, b, P);
    else {
      arma::mat omega = draw_omega(a, b, theta);
      b = draw_b(kappa, omega, theta, a, P);
      a = draw_a(kappa, omega, theta, b, P);
    }
    out(it, 0) = a(0); out(it, 1) = b(0);
  }
  return out;
}
