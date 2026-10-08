// R interface to the building blocks of the samplers (see R/blocks.R). Each function draws one block
// from its full conditional with R's RNG, so a sampler written in R from these blocks reproduces the
// C++ samplers step by step.
// [[Rcpp::depends(RcppArmadillo)]]
#include "blocks.h"
#include "rng.h"

// [[Rcpp::export(.blk_gauss)]]
arma::vec blk_gauss(const arma::vec& prec, const arma::vec& num) {
  if (prec.n_elem != num.n_elem) Rcpp::stop("prec and num must have the same length");
  arma::vec x(prec.n_elem);
  for (arma::uword j = 0; j < x.n_elem; ++j) x(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
  return x;
}

// [[Rcpp::export(.blk_invgamma)]]
arma::vec blk_invgamma(const arma::vec& shape, const arma::vec& scale) {
  arma::vec x(shape.n_elem);
  for (arma::uword j = 0; j < x.n_elem; ++j) x(j) = rinvgamma(shape(j), scale(j));
  return x;
}

// [[Rcpp::export(.blk_tnorm)]]
arma::vec blk_tnorm(const arma::vec& mean, const arma::vec& sd, double lo, double hi) {
  arma::vec x(mean.n_elem);
  for (arma::uword j = 0; j < x.n_elem; ++j) x(j) = rtnorm(mean(j), sd(j), lo, hi);
  return x;
}

// [[Rcpp::export(.blk_omega)]]
arma::mat blk_omega(const arma::vec& a, const arma::vec& b, const arma::vec& theta) { return draw_omega(a, b, theta); }

// [[Rcpp::export(.blk_pg)]]
arma::mat blk_pg(const arma::mat& eta) {
  arma::mat w(eta.n_rows, eta.n_cols);
  for (arma::uword k = 0; k < eta.n_cols; ++k)
    for (arma::uword j = 0; j < eta.n_rows; ++j) w(j, k) = rpg1(eta(j, k));
  return w;
}

// Albert & Chib (1993): latent normal responses z ~ N(eta, 1) truncated to z > 0 when y = 1, z < 0 when
// y = 0 (unrestricted when y is missing)
// [[Rcpp::export(.blk_probit_z)]]
arma::mat blk_probit_z(const arma::mat& Y, const arma::mat& eta) {
  if (Y.n_rows != eta.n_rows || Y.n_cols != eta.n_cols) Rcpp::stop("Y and eta must have the same dimensions");
  arma::mat z(Y.n_rows, Y.n_cols);
  for (arma::uword k = 0; k < Y.n_cols; ++k)
    for (arma::uword j = 0; j < Y.n_rows; ++j) {
      const double y = Y(j, k), m = eta(j, k);
      z(j, k) = ISNAN(y) ? m + R::norm_rand() : (y > 0.5 ? rtnorm(m, 1.0, 0.0, R_PosInf) : rtnorm(m, 1.0, R_NegInf, 0.0));
    }
  return z;
}

// [[Rcpp::export(.blk_irt_part)]]
Rcpp::List blk_irt_part(const arma::mat& Y, const arma::mat& omega, const arma::vec& a, const arma::vec& b) {
  arma::vec prec, num;
  irt_part(Y - 0.5, omega, a, b, prec, num);
  return Rcpp::List::create(Rcpp::Named("prec") = prec, Rcpp::Named("num") = num);
}

// per-item prior means of b (b_mu empty: the mean of the prior for every item)
static Prior item_prior(const Prior& P, const arma::vec& b_mu, arma::uword k) {
  Prior Q = P; if (b_mu.n_elem) Q.b_mu = b_mu(k); return Q;
}

// [[Rcpp::export(.blk_items_pg)]]
Rcpp::List blk_items_pg(const arma::mat& Y, const arma::mat& omega, const arma::vec& theta, arma::vec a,
                        const Rcpp::NumericVector& prior, bool one_pl, const arma::vec& b_mu) {
  const Prior P = make_prior(prior); const arma::mat kappa = Y - 0.5;
  arma::vec b(a.n_elem);
  if (!b_mu.n_elem) b = draw_b(kappa, omega, theta, a, P);
  else for (arma::uword k = 0; k < a.n_elem; ++k)
    b(k) = draw_b(kappa.col(k), omega.col(k), theta, a.subvec(k, k), item_prior(P, b_mu, k))(0);
  if (!one_pl) a = draw_a(kappa, omega, theta, b, P);
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("b") = b);
}

// [[Rcpp::export(.blk_items_joint)]]
Rcpp::List blk_items_joint(const arma::mat& Y, const arma::vec& theta, arma::vec a, arma::vec b, const Rcpp::NumericVector& prior,
                           const arma::vec& b_mu) {
  const Prior P = make_prior(prior); int acc = 0;
  if (!b_mu.n_elem) acc = draw_ab_joint(Y, theta, a, b, P);
  else for (arma::uword k = 0; k < a.n_elem; ++k) {
    arma::vec ak = a.subvec(k, k), bk = b.subvec(k, k);
    acc += draw_ab_joint(Y.col(k), theta, ak, bk, item_prior(P, b_mu, k)); a(k) = ak(0); b(k) = bk(0);
  }
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("b") = b, Rcpp::Named("accepted") = acc);
}

// [[Rcpp::export(.blk_lambda)]]
arma::vec blk_lambda(const arma::mat& L, const arma::mat& m, const arma::mat& k1nu, const arma::mat& W, double mu_l, double sd_l,
                     double lo_l) {
  return draw_lambda(L, m, k1nu, W, mu_l, sd_l, lo_l);
}

// [[Rcpp::export(.blk_sigma2)]]
arma::vec blk_sigma2(const arma::mat& L, const arma::mat& m, const arma::vec& lambda, const Rcpp::NumericVector& prior,
                     Rcpp::Nullable<Rcpp::NumericMatrix> nu, double q, const arma::vec& current) {
  const Prior P = make_prior(prior);
  if (nu.isNull()) return draw_sigma2_normal(L, m, lambda, P, current);
  return draw_sigma2_ald(L, m, lambda, Rcpp::as<arma::mat>(nu.get()), ald_constants(q), P, current);
}

// [[Rcpp::export(.blk_nu)]]
arma::mat blk_nu(const arma::mat& r, const arma::mat& s, double q) {
  const Ald c = ald_constants(q); arma::mat nu(r.n_rows, r.n_cols);
  for (arma::uword k = 0; k < r.n_cols; ++k)
    for (arma::uword j = 0; j < r.n_rows; ++j) nu(j, k) = draw_nu(r(j, k), s(j, k), c);
  return nu;
}

// [[Rcpp::export(.blk_theta_collapsed)]]
arma::vec blk_theta_collapsed(const arma::vec& prec, const arma::vec& num, const arma::mat& W, const arma::mat& y,
                              const arma::vec& rho, const arma::vec& m0, const arma::vec& g, const arma::vec& tau) {
  return draw_theta_collapsed(prec, num, W, y, rho, m0, g, tau);
}

// [[Rcpp::export(.blk_rt_items_collapsed)]]
Rcpp::List blk_rt_items_collapsed(const arma::mat& L, const arma::mat& k1nu, const arma::mat& W, const arma::vec& theta,
                                  const arma::vec& mz, const arma::vec& tau, bool with_rho, double mu_l, double sd_l,
                                  double rho_sd, arma::vec lambda, arma::vec rho, double lo_l) {
  draw_rt_items_collapsed(L, k1nu, W, theta, mz, tau, with_rho, mu_l, sd_l, rho_sd, lambda, rho, lo_l);
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda, Rcpp::Named("rho") = rho);
}

// [[Rcpp::export(.blk_var)]]
double blk_var(double alpha, double beta, double fam, double p1, double p2, double current) {
  return draw_var(alpha, beta, fam, p1, p2, current);
}

// Multivariate normal draws in precision form, one per slice: x_n ~ N(Q_n^-1 h_n, Q_n^-1), optionally
// truncated to the box [lower, upper] by rejection (up to max_tries; otherwise the slice is NA)
// [[Rcpp::export(.blk_mvn)]]
arma::mat blk_mvn(const arma::cube& Q, const arma::mat& h, const arma::vec& lower, const arma::vec& upper, int max_tries) {
  const arma::uword p = Q.n_rows, n = Q.n_slices;
  arma::mat x(p, n);
  for (arma::uword s = 0; s < n; ++s) {
    arma::mat U = arma::chol(arma::symmatu(Q.slice(s)));             // Q = U' U
    arma::vec m = arma::solve(arma::trimatu(U), arma::solve(arma::trimatl(U.t()), h.col(s)));
    bool ok = false; arma::vec y(p);
    for (int t = 0; t < max_tries && !ok; ++t) {
      arma::vec z(p); for (arma::uword k = 0; k < p; ++k) z(k) = R::norm_rand();
      y = m + arma::solve(arma::trimatu(U), z);
      ok = arma::all(y >= lower) && arma::all(y <= upper);
    }
    x.col(s) = ok ? y : arma::vec(p, arma::fill::value(NA_REAL));
  }
  return x;
}

// (a_i, b_i) of every item with the latent responses integrated out: item_mh_step() (see blocks.cpp).
// link 0 = logit, 1 = probit; form 0: eta = a (theta - b), 1: eta = a theta - b. Prior N((a, b); mean_i, P^-1)
// with a > 0.
// [[Rcpp::export(.blk_items_mh)]]
Rcpp::List blk_items_mh(const arma::mat& Y, const arma::vec& theta, arma::vec a, arma::vec b, const arma::mat& mean,
                        const arma::mat& P, int link, int form) {
  const arma::uword K = Y.n_cols; arma::uvec acc(K, arma::fill::zeros);
  for (arma::uword i = 0; i < K; ++i) {
    const arma::vec m = mean.row(i).t();
    const ItemPrior prior{
      [](const arma::vec& x) { return x(0) > 0.0; },
      [&m, &P](const arma::vec& x, double& f, arma::vec& g, arma::mat& G) {
        const arma::vec dv = x - m; f -= 0.5 * arma::dot(dv, P * dv); g -= P * dv; G += P;
      }};
    arma::vec x = {a(i), b(i)};
    if (item_mh_step(Y.col(i), theta, x, link, form, prior)) { a(i) = x(0); b(i) = x(1); acc(i) = 1; }
  }
  return Rcpp::List::create(Rcpp::Named("a") = a, Rcpp::Named("b") = b, Rcpp::Named("accepted") = acc);
}
