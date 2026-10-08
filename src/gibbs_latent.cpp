// Latent regression of speed on ability and covariates (dissertation Chapter 4).
//   logit P(y_jk = 1) = a_k (theta_j - b_k),   log T_jk = lambda_k - zeta_j + e_jk, e ~ N(0, sigma2_k)
//   theta_j ~ N(0, 1)
//   zeta_j = x_j' beta + gamma theta_j + u_j,  u ~ N(0, s) (mean regression)
//                                         or  u ~ ALD(q, scale s): u = k1 nu + sqrt(k2 s nu) z, nu ~ Exp(mean s)
//   (beta, gamma) ~ N(0, beta_sd^2 I), s with the speed-variance prior (priors: Prior). X may contain an intercept column.
// The update order follows ExtendedRtIrtModeling.jl (GibbsRtIrtLatentQr).
// collapse = true: (beta, gamma), theta and lambda are drawn with zeta integrated out, each followed
// (eventually) by a draw of zeta from its full conditional (steps 2c, 5c-8c); the location move is skipped.
// [[Rcpp::depends(RcppArmadillo)]]
#include "blocks.h"
#include "rng.h"

// [[Rcpp::export(.gibbs_latent)]]
Rcpp::List gibbs_latent(const arma::mat& Y, const arma::mat& L, const arma::mat& X, bool ald, double q, bool one_pl,
                        bool collapse, bool ab_joint, int n_iter, int n_burn, Rcpp::List init,
                       Rcpp::NumericVector prior) {
  const Prior Pr = make_prior(prior);
  const arma::uword N = Y.n_rows, K = Y.n_cols, P = X.n_cols;
  const arma::mat kappa = Y - 0.5;
  const double mu_l = ISNAN(Pr.l_mu) ? arma::mean(arma::vectorise(L)) : Pr.l_mu,
               sd_l = ISNAN(Pr.l_sd) ? arma::stddev(arma::vectorise(L)) : Pr.l_sd;
  const Ald c = ald ? ald_constants(q) : Ald{0.0, 1.0};

  arma::vec theta = init["theta"], zeta = init["zeta"], a = init["a"], b = init["b"];
  arma::vec lambda = init["lambda"], sigma2 = init["sigma2"], beta = init["beta"];   // beta = (beta_x, gamma)
  double s = init["s"];
  arma::vec nu(N, arma::fill::ones), w(N), shift(N, arma::fill::zeros), prec, num;
  arma::mat omega(N, K);
  auto update_weights = [&]() { w = ald ? 1.0 / (c.k2 * s * nu) : arma::vec(N, arma::fill::value(1.0 / s)); shift = c.k1 * nu; };
  update_weights();

  arma::mat draws(n_iter, 4 * K + P + 2);
  arma::vec loglik(n_iter), theta_sum(N, arma::fill::zeros), theta_ss(N, arma::fill::zeros);
  arma::vec zeta_sum(N, arma::fill::zeros), zeta_ss(N, arma::fill::zeros), nu_sum(N, arma::fill::zeros);

  // RT person fit (Marianti et al., 2014): the measurement model of log T is normal here
  arma::vec lt_sum(N, arma::fill::zeros), lt_hit(N, arma::fill::zeros);
  const double lt_crit = R::qchisq(0.95, K, 1, 0);
  int ab_accepted = 0;
  for (int it = 0; it < n_iter; ++it) {
    Rcpp::checkUserInterrupt();
    arma::mat Z = arma::join_rows(X, theta);           // design of the structural regression
    // 1. ALD mixing variables: residual r = zeta - x' beta - gamma theta
    if (ald) {
      arma::vec r = zeta - Z * beta;
      for (arma::uword j = 0; j < N; ++j) nu(j) = draw_nu(r(j), s, c);
      update_weights();
    }
    if (collapse) {
      // 2c. (beta, gamma) with zeta integrated out: the RT part of zeta_j is the pseudo-observation
      //     zhat_j = sum_k (lambda_k - log T_jk) / sigma2_k / P_r ~ N(zeta_j, 1 / P_r), P_r = sum_k 1 / sigma2_k,
      //     so zhat_j - k1 nu_j ~ N(z_j' beta, 1 / w_j + 1 / P_r); then zeta from its full conditional
      const double Prt = arma::accu(1.0 / sigma2);
      arma::mat yz = -L;
      yz.each_row() += lambda.t();
      arma::vec zhat = yz * (1.0 / sigma2) / Prt, wc = 1.0 / (1.0 / w + 1.0 / Prt);
      beta = rmvn_prec(arma::eye(P + 1, P + 1) / (Pr.beta_sd * Pr.beta_sd) + Z.t() * (Z.each_col() % wc), Z.t() * (wc % (zhat - shift)));
      arma::vec mz = Z * beta + shift, pz = w + Prt;
      for (arma::uword j = 0; j < N; ++j) zeta(j) = (w(j) * mz(j) + Prt * zhat(j)) / pz(j) + R::norm_rand() / std::sqrt(pz(j));
    } else {
    // 2. (beta, gamma) jointly: weighted regression of zeta - k1 nu on Z, prior N(0, I)
    beta = rmvn_prec(arma::eye(P + 1, P + 1) / (Pr.beta_sd * Pr.beta_sd) + Z.t() * (Z.each_col() % w), Z.t() * (w % (zeta - shift)));
    }
    // 3. residual scale s (normal variance or ALD scale)
    arma::vec r = zeta - Z * beta;
    if (ald) {
      const double A = arma::accu(arma::square(r - shift) / (2.0 * c.k2 * nu)), Bn = arma::accu(nu);
      s = Pr.vs_fam == 0.0 ? rinvgamma(Pr.vs_p1 + 1.5 * N, Pr.vs_p2 + A + Bn) : draw_var(1.5 * N, A + Bn, Pr.vs_fam, Pr.vs_p1, Pr.vs_p2, s);
    } else s = draw_var(0.5 * N, 0.5 * arma::dot(r, r), Pr.vs_fam, Pr.vs_p1, Pr.vs_p2, s);
    update_weights();
    // 4. accuracy: Polya-Gamma variables, difficulty, discrimination
    // collapse: first (a, b) jointly with omega integrated out (moves along the a-b ridge), then the
    // Polya-Gamma steps as usual (they always move, so the chain cannot stall in the tails)
    if (ab_joint && !one_pl) ab_accepted += draw_ab_joint(Y, theta, a, b, Pr);
    omega = draw_omega(a, b, theta);
    b = draw_b(kappa, omega, theta, a, Pr);
    if (!one_pl) a = draw_a(kappa, omega, theta, b, Pr);
    if (collapse) {
      const arma::vec tau = 1.0 / w, m0 = X * beta.head(P) + shift, g(N, arma::fill::value(beta(P)));
      arma::mat Wt = arma::repmat((1.0 / sigma2).t(), N, 1), y = -L;
      y.each_row() += lambda.t();                      // y = lambda - log T = zeta + e
      // 5c. theta_j with zeta_j ~ N(x'beta + gamma theta + k1 nu, 1 / w) integrated out
      irt_part(kappa, omega, a, b, prec, num);
      theta = draw_theta_collapsed(prec + 1.0, num, Wt, y, arma::zeros(K), m0, g, tau);
      // 6c. lambda, zeta integrated out
      arma::vec rho0;
      draw_rt_items_collapsed(L, arma::zeros(N, K), Wt, theta, m0 + g % theta, tau, false, mu_l, sd_l, Pr.rho_sd, lambda, rho0, Pr.l_lo);
      // 7c. zeta_j from its full conditional (as step 8)
      arma::vec mz = m0 + g % theta;
      arma::mat yz = -L;
      yz.each_row() += lambda.t();
      prec = w + arma::accu(1.0 / sigma2);
      num = w % mz + yz * (1.0 / sigma2);
      for (arma::uword j = 0; j < N; ++j) zeta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
      // 8c. sigma2_k
      sigma2 = draw_sigma2_normal(L, arma::repmat(-zeta, 1, K), lambda, Pr, sigma2);
    } else {
    // 5. theta_j: N(0, 1) prior + accuracy part + structural part (theta is a regressor of zeta)
    double gamma = beta(P);
    arma::vec rest = X * beta.head(P) + shift;
    irt_part(kappa, omega, a, b, prec, num);
    prec += 1.0 + gamma * gamma * w;
    num += gamma * w % (zeta - rest);
    for (arma::uword j = 0; j < N; ++j) theta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
    // 6.-7. lambda_k and sigma2_k (normal log RT, m = -zeta)
    arma::mat m = arma::repmat(-zeta, 1, K), Wt = arma::repmat((1.0 / sigma2).t(), N, 1);
    lambda = draw_lambda(L, m, arma::zeros(N, K), Wt, mu_l, sd_l, Pr.l_lo);
    sigma2 = draw_sigma2_normal(L, m, lambda, Pr, sigma2);
    // 8. zeta_j: structural prior N(x'beta + gamma theta + k1 nu, 1 / w) + RT part: lambda - log T = zeta + e
    arma::vec mz = X * beta.head(P) + beta(P) * theta + shift;
    arma::mat yz = -L;
    yz.each_row() += lambda.t();
    prec = w + arma::accu(1.0 / sigma2);
    num = w % mz + yz * (1.0 / sigma2);
    for (arma::uword j = 0; j < N; ++j) zeta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
    // 9. location move: common shift of (lambda, zeta); zeta_j has prior N(mz_j, 1 / w_j)
    draw_location_shift(lambda, zeta, mz, w, mu_l, sd_l, Pr.l_lo);
    }

    // complete-data log-likelihood (for DIC) and storage
    const arma::vec mz = X * beta.head(P) + beta(P) * theta + shift;
    arma::mat eta = theta * a.t(); eta.each_row() -= (a % b).t();
    arma::mat res = L; res.each_row() -= lambda.t(); res.each_col() += zeta;
    loglik(it) = arma::accu(Y % eta - arma::log1p(arma::exp(eta))) +
                 -0.5 * arma::accu(arma::square(res) * (1.0 / sigma2)) - 0.5 * N * arma::accu(arma::log(2 * M_PI * sigma2)) +
                 arma::accu(-0.5 * theta % theta) - 0.5 * N * std::log(2 * M_PI) +
                 arma::accu(0.5 * arma::log(w / (2 * M_PI)) - 0.5 * w % arma::square(zeta - mz));
    draws.row(it) = arma::join_cols(arma::join_cols(a, b, lambda, sigma2), beta, arma::vec{s}).t();
    if (it >= n_burn) {
      theta_sum += theta; theta_ss += theta % theta; zeta_sum += zeta; zeta_ss += zeta % zeta;
      if (ald) nu_sum += nu;
      arma::vec lt = arma::square(res) * (1.0 / sigma2); lt_sum += lt; lt_hit += arma::conv_to<arma::vec>::from(lt > lt_crit);
    }
  }
  double n_keep = n_iter - n_burn;
  return Rcpp::List::create(
    Rcpp::Named("ab_accept") = ab_joint && !one_pl ? ab_accepted / double(n_iter * K) : NA_REAL,
    Rcpp::Named("draws") = draws, Rcpp::Named("loglik") = loglik,
    Rcpp::Named("theta_mean") = theta_sum / n_keep, Rcpp::Named("theta_sd") = arma::sqrt(theta_ss / n_keep - arma::square(theta_sum / n_keep)),
    Rcpp::Named("zeta_mean") = zeta_sum / n_keep, Rcpp::Named("zeta_sd") = arma::sqrt(zeta_ss / n_keep - arma::square(zeta_sum / n_keep)),
    Rcpp::Named("nu_mean") = nu_sum / n_keep,
    Rcpp::Named("lt_mean") = lt_sum / n_keep, Rcpp::Named("lt_exceed") = lt_hit / n_keep,
    Rcpp::Named("last") = Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("zeta") = zeta, Rcpp::Named("a") = a,
                                             Rcpp::Named("b") = b, Rcpp::Named("lambda") = lambda, Rcpp::Named("sigma2") = sigma2,
                                             Rcpp::Named("beta") = beta, Rcpp::Named("s") = s));
}
