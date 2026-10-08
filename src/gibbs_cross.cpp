// Cross-relation RT-IRT (dissertation Chapter 5), mean (normal) or quantile (ALD) version.
//   logit P(y_jk = 1) = a_k (theta_j - b_k)
//   log T_jk = lambda_k - zeta_j - rho_k theta_j + e_jk
//   e_jk ~ N(0, sigma2_k)  or  ALD(q): e = k1 nu + sqrt(k2 sigma2_k nu) z, nu ~ Exp(mean sigma2_k)
//   theta_j ~ N(0, 1), zeta_j ~ N(0, s) with s = 1 (cov2one) or s ~ IG; rho_k ~ N(0, rho_sd^2) (priors: Prior)
// One call runs one chain; the update order follows ExtendedRtIrtModeling.jl (GibbsRtIrtCrossQr).
// collapse = true: theta and (lambda, rho) are drawn with zeta integrated out, zeta last (steps 5c-8c);
// this covers the directions of the two group moves, which are then skipped.
// [[Rcpp::depends(RcppArmadillo)]]
#include "blocks.h"
#include "rng.h"

// [[Rcpp::export(.gibbs_cross)]]
Rcpp::List gibbs_cross(const arma::mat& Y, const arma::mat& L, bool ald, double q, bool cov2one, bool one_pl,
                       bool collapse, bool ab_joint, int n_iter, int n_burn, Rcpp::List init,
                       Rcpp::NumericVector prior) {
  const Prior Pr = make_prior(prior);
  const arma::uword N = Y.n_rows, K = Y.n_cols;
  const arma::mat kappa = Y - 0.5;
  const double mu_l = ISNAN(Pr.l_mu) ? arma::mean(arma::vectorise(L)) : Pr.l_mu,
               sd_l = ISNAN(Pr.l_sd) ? arma::stddev(arma::vectorise(L)) : Pr.l_sd;
  const Ald c = ald ? ald_constants(q) : Ald{0.0, 1.0};

  arma::vec theta = init["theta"], zeta = init["zeta"], a = init["a"], b = init["b"];
  arma::vec lambda = init["lambda"], sigma2 = init["sigma2"], rho = init["rho"];
  double s = init["s"];
  arma::mat nu(N, K, arma::fill::ones), omega(N, K), W(N, K), k1nu(N, K, arma::fill::zeros);
  auto update_weights = [&]() {                       // W = 1 / Var(e_jk), k1nu = k1 nu_jk
    for (arma::uword k = 0; k < K; ++k) {
      W.col(k) = ald ? 1.0 / (c.k2 * sigma2(k) * nu.col(k)) : arma::vec(N, arma::fill::value(1.0 / sigma2(k)));
      k1nu.col(k) = c.k1 * nu.col(k);
    }
  };
  update_weights();

  arma::mat draws(n_iter, 5 * K + 1);
  arma::vec loglik(n_iter), theta_sum(N, arma::fill::zeros), theta_ss(N, arma::fill::zeros);
  arma::vec zeta_sum(N, arma::fill::zeros), zeta_ss(N, arma::fill::zeros);
  arma::mat nu_sum(N, K, arma::fill::zeros);
  arma::vec prec, num;

  // RT person fit (Marianti et al., 2014); chi-square with K df for normal log T only
  arma::vec lt_sum(N, arma::fill::zeros), lt_hit(N, arma::fill::zeros);
  const double lt_crit = R::qchisq(0.95, K, 1, 0);
  int ab_accepted = 0;
  for (int it = 0; it < n_iter; ++it) {
    Rcpp::checkUserInterrupt();
    // 1. ALD mixing variables: residual r = log T - lambda + zeta + rho theta
    if (ald) {
      for (arma::uword k = 0; k < K; ++k)
        for (arma::uword j = 0; j < N; ++j)
          nu(j, k) = draw_nu(L(j, k) - lambda(k) + zeta(j) + rho(k) * theta(j), sigma2(k), c);
      update_weights();
    }
    // 2. rho_k: y = lambda - zeta + k1 nu - log T = rho_k theta + error
    if (!collapse) for (arma::uword k = 0; k < K; ++k) {
      arma::vec y = lambda(k) - zeta + k1nu.col(k) - L.col(k);
      double p = 1.0 / (Pr.rho_sd * Pr.rho_sd) + arma::dot(W.col(k), theta % theta);
      rho(k) = arma::dot(W.col(k), theta % y) / p + R::norm_rand() / std::sqrt(p);
    }
    // 3. speed variance
    if (!cov2one) s = draw_var(0.5 * N, 0.5 * arma::dot(zeta, zeta), Pr.vs_fam, Pr.vs_p1, Pr.vs_p2, s);
    // 4. accuracy: Polya-Gamma variables, item difficulty and discrimination
    // collapse: first (a, b) jointly with omega integrated out (moves along the a-b ridge), then the
    // Polya-Gamma steps as usual (they always move, so the chain cannot stall in the tails)
    if (ab_joint && !one_pl) ab_accepted += draw_ab_joint(Y, theta, a, b, Pr);
    omega = draw_omega(a, b, theta);
    b = draw_b(kappa, omega, theta, a, Pr);
    if (!one_pl) a = draw_a(kappa, omega, theta, b, Pr);
    if (collapse) {
      const arma::vec zero(N, arma::fill::zeros), tau(N, arma::fill::value(s));
      // 5c. theta_j with zeta_j ~ N(0, s) integrated out: y = lambda + k1 nu - log T = zeta + rho theta + e
      irt_part(kappa, omega, a, b, prec, num);
      arma::mat y = k1nu - L;
      y.each_row() += lambda.t();
      theta = draw_theta_collapsed(prec + 1.0, num, W, y, rho, zero, zero, tau);
      // 6c. (lambda, rho) jointly, zeta integrated out
      draw_rt_items_collapsed(L, k1nu, W, theta, zero, tau, true, mu_l, sd_l, Pr.rho_sd, lambda, rho, Pr.l_lo);
      // 7c. zeta_j from its full conditional (as step 8)
      arma::mat yz = k1nu - L - theta * rho.t();
      yz.each_row() += lambda.t();
      prec = 1.0 / s + arma::sum(W, 1);
      num = arma::sum(W % yz, 1);
      for (arma::uword j = 0; j < N; ++j) zeta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
      // 8c. sigma2_k
      arma::mat m = -theta * rho.t();
      m.each_col() -= zeta;
      sigma2 = ald ? draw_sigma2_ald(L, m, lambda, nu, c, Pr, sigma2) : draw_sigma2_normal(L, m, lambda, Pr, sigma2);
      update_weights();
    } else {
    // 5. theta_j: N(0, 1) prior + accuracy part + RT part (theta enters log T through rho_k)
    irt_part(kappa, omega, a, b, prec, num);
    arma::mat y = k1nu - L;                              // lambda - log T + k1 nu - zeta
    y.each_row() += lambda.t();
    y.each_col() -= zeta;
    prec += 1.0 + W * (rho % rho);
    num += (W % y) * rho;
    for (arma::uword j = 0; j < N; ++j) theta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
    // 6. lambda_k: log T = lambda + m + k1 nu + e with m = -zeta - rho theta
    arma::mat m = -theta * rho.t();
    m.each_col() -= zeta;
    lambda = draw_lambda(L, m, k1nu, W, mu_l, sd_l, Pr.l_lo);
    // 7. sigma2_k (normal residual variance or ALD scale)
    sigma2 = ald ? draw_sigma2_ald(L, m, lambda, nu, c, Pr, sigma2) : draw_sigma2_normal(L, m, lambda, Pr, sigma2);
    update_weights();
    // 8. zeta_j: N(0, s) prior + RT part: lambda - rho theta + k1 nu - log T = zeta + error
    arma::mat yz = k1nu - L - theta * rho.t();
    yz.each_row() += lambda.t();
    prec = 1.0 / s + arma::sum(W, 1);
    num = arma::sum(W % yz, 1);
    for (arma::uword j = 0; j < N; ++j) zeta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
    // 9. location move: common shift of (lambda, zeta); zeta_j has prior N(0, s)
    draw_location_shift(lambda, zeta, arma::zeros(N), arma::vec(N, arma::fill::value(1.0 / s)), mu_l, sd_l, Pr.l_lo);
    // 10. common cross-relation move: (rho_k - d, zeta_j + d theta_j) leaves log T unchanged, so d is
    //     identified only by the priors rho_k ~ N(0, rho_sd^2), zeta_j ~ N(0, s) (independent of theta).
    //     Exact Gibbs draw of d (shear transformation, Jacobian 1): prec = K / rho_sd^2 + sum theta^2 / s.
    {
      double vr = Pr.rho_sd * Pr.rho_sd, p = K / vr + arma::dot(theta, theta) / s;
      double d = (arma::accu(rho) / vr - arma::dot(zeta, theta) / s) / p + R::norm_rand() / std::sqrt(p);
      rho -= d; zeta += d * theta;
    }
    }

    // complete-data log-likelihood (for DIC) and storage
    arma::mat eta = theta * a.t(); eta.each_row() -= (a % b).t();
    arma::mat mu = -theta * rho.t() + k1nu; mu.each_row() += lambda.t(); mu.each_col() -= zeta;
    loglik(it) = arma::accu(Y % eta - arma::log1p(arma::exp(eta))) +
                 arma::accu(0.5 * arma::log(W) - 0.5 * W % arma::square(L - mu)) - 0.5 * N * K * std::log(2 * M_PI) +
                 arma::accu(-0.5 * theta % theta) - 0.5 * N * std::log(2 * M_PI) +
                 arma::accu(-0.5 * zeta % zeta / s) - 0.5 * N * std::log(2 * M_PI * s);
    draws.row(it) = arma::join_cols(arma::join_cols(a, b, lambda, sigma2), rho, arma::vec{s}).t();
    if (it >= n_burn) {
      theta_sum += theta; theta_ss += theta % theta; zeta_sum += zeta; zeta_ss += zeta % zeta;
      if (ald) nu_sum += nu;
      arma::vec lt = arma::sum(W % arma::square(L - mu), 1); lt_sum += lt; lt_hit += arma::conv_to<arma::vec>::from(lt > lt_crit);
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
                                             Rcpp::Named("rho") = rho, Rcpp::Named("s") = s));
}
