// Joint RT-IRT with a latent regression of ability and speed (dissertation Chapter 3); with no
// covariates it is the null (hierarchical) model of van der Linden (2007).
//   logit P(y_jk = 1) = a_k (theta_j - b_k),   log T_jk = lambda_k - zeta_j + e_jk, e ~ N(0, sigma2_k)
//   theta_j = x_j' beta1 + e1_j,                    e1 ~ N(0, 1)
//   zeta_j  = x_j' beta2 + c e1_j + e2_j,           e2 ~ N(0, v)
// so Var(theta | x) = 1, Var(zeta | x) = c^2 + v and Corr = c / sqrt(c^2 + v).
//   cov2one = false: c ~ N(0, c_sd^2), v with the speed-variance prior (speed variance free)
//   cov2one = true:  Var(zeta | x) = 1, c = rho ~ U(-1, 1), v = 1 - rho^2 (random-walk Metropolis on atanh rho)
//   beta1, beta2 ~ N(0, beta_sd^2 I) (priors: Prior). X may contain an intercept column or have no columns.
// collapse = true: theta and lambda are drawn with zeta integrated out, zeta last (steps 5c-8c);
// the location move is then skipped.
// [[Rcpp::depends(RcppArmadillo)]]
#include "blocks.h"
#include "rng.h"

// log-likelihood of rho for e2 = r2 - rho e1 ~ N(0, 1 - rho^2), plus the Jacobian of rho = tanh(z)
static double ll_rho(double rho, const arma::vec& e1, const arma::vec& r2) {
  double v = 1.0 - rho * rho;
  return -0.5 * e1.n_elem * std::log(v) - arma::accu(arma::square(r2 - rho * e1)) / (2.0 * v) + std::log(v);
}

// [[Rcpp::export(.gibbs_rtirt)]]
Rcpp::List gibbs_rtirt(const arma::mat& Y, const arma::mat& L, const arma::mat& X, bool cov2one, bool one_pl,
                       bool collapse, bool ab_joint, int n_iter, int n_burn, Rcpp::List init,
                       Rcpp::NumericVector prior) {
  const Prior Pr = make_prior(prior);
  const arma::uword N = Y.n_rows, K = Y.n_cols, P = X.n_cols;
  const arma::mat kappa = Y - 0.5, XX = X.t() * X, I_P = arma::eye(P, P) / (Pr.beta_sd * Pr.beta_sd);
  const double mu_l = ISNAN(Pr.l_mu) ? arma::mean(arma::vectorise(L)) : Pr.l_mu,
               sd_l = ISNAN(Pr.l_sd) ? arma::stddev(arma::vectorise(L)) : Pr.l_sd;

  arma::vec theta = init["theta"], zeta = init["zeta"], a = init["a"], b = init["b"];
  arma::vec lambda = init["lambda"], sigma2 = init["sigma2"], beta1 = init["beta1"], beta2 = init["beta2"];
  double cc = init["c"], v = init["v"];
  arma::mat omega(N, K);
  arma::vec prec, num;
  int accepted = 0;

  arma::mat draws(n_iter, 4 * K + 2 * P + 2);
  arma::vec loglik(n_iter), theta_sum(N, arma::fill::zeros), theta_ss(N, arma::fill::zeros);
  arma::vec zeta_sum(N, arma::fill::zeros), zeta_ss(N, arma::fill::zeros);

  // RT person fit (Marianti et al., 2014): l^t_j = sum_k (log T - lambda + zeta)^2 / sigma2_k ~ chi2(K)
  arma::vec lt_sum(N, arma::fill::zeros), lt_hit(N, arma::fill::zeros);
  const double lt_crit = R::qchisq(0.95, K, 1, 0);
  int ab_accepted = 0;
  for (int it = 0; it < n_iter; ++it) {
    Rcpp::checkUserInterrupt();
    // 1. beta1: theta - X beta1 ~ N(0, 1) and (zeta - X beta2 - c theta) + c X beta1 ~ N(0, v)
    if (P > 0) {
      arma::vec r = zeta - X * beta2 - cc * theta;
      beta1 = rmvn_prec(I_P + XX * (1.0 + cc * cc / v), X.t() * theta - (cc / v) * (X.t() * r));
    }
    arma::vec e1 = theta - X * beta1;
    // 2. beta2: zeta - c e1 - X beta2 ~ N(0, v)
    if (P > 0) beta2 = rmvn_prec(I_P + XX / v, X.t() * (zeta - cc * e1) / v);
    arma::vec r2 = zeta - X * beta2;
    // 3. person covariance: regression of r2 on e1
    if (!cov2one) {
      double p = 1.0 / (Pr.c_sd * Pr.c_sd) + arma::dot(e1, e1) / v;
      cc = arma::dot(e1, r2) / v / p + R::norm_rand() / std::sqrt(p);
      v = draw_var(0.5 * N, 0.5 * arma::accu(arma::square(r2 - cc * e1)), Pr.vs_fam, Pr.vs_p1, Pr.vs_p2, v);
    } else {
      for (int rep = 0; rep < 5; ++rep) {
        double prop = std::tanh(std::atanh(cc) + 0.05 * R::norm_rand());
        if (std::log(R::unif_rand()) < ll_rho(prop, e1, r2) - ll_rho(cc, e1, r2)) { cc = prop; ++accepted; }
      }
      v = 1.0 - cc * cc;
    }
    // 4. accuracy: Polya-Gamma variables, difficulty, discrimination
    // collapse: first (a, b) jointly with omega integrated out (moves along the a-b ridge), then the
    // Polya-Gamma steps as usual (they always move, so the chain cannot stall in the tails)
    if (ab_joint && !one_pl) ab_accepted += draw_ab_joint(Y, theta, a, b, Pr);
    omega = draw_omega(a, b, theta);
    b = draw_b(kappa, omega, theta, a, Pr);
    if (!one_pl) a = draw_a(kappa, omega, theta, b, Pr);
    arma::vec mu1 = X * beta1, mu2 = X * beta2;
    if (collapse) {
      const arma::vec tau(N, arma::fill::value(v)), m0 = mu2 - cc * mu1, g(N, arma::fill::value(cc));
      arma::mat Wt = arma::repmat((1.0 / sigma2).t(), N, 1), y = -L;
      y.each_row() += lambda.t();                      // y = lambda - log T = zeta + e
      // 5c. theta_j with zeta_j | theta_j ~ N(mu2 + c (theta - mu1), v) integrated out
      irt_part(kappa, omega, a, b, prec, num);
      theta = draw_theta_collapsed(prec + 1.0, num + mu1, Wt, y, arma::zeros(K), m0, g, tau);
      // 6c. lambda, zeta integrated out
      arma::vec rho0, mz = m0 + cc * theta;
      draw_rt_items_collapsed(L, arma::zeros(N, K), Wt, theta, mz, tau, false, mu_l, sd_l, Pr.rho_sd, lambda, rho0, Pr.l_lo);
      // 7c. zeta_j from its full conditional (as step 8)
      arma::mat yz = -L;
      yz.each_row() += lambda.t();
      prec = 1.0 / v + arma::accu(1.0 / sigma2) * arma::ones(N);
      num = mz / v + yz * (1.0 / sigma2);
      for (arma::uword j = 0; j < N; ++j) zeta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
      // 8c. sigma2_k
      sigma2 = draw_sigma2_normal(L, arma::repmat(-zeta, 1, K), lambda, Pr, sigma2);
    } else {
    // 5. theta_j: prior N(mu1, 1) and zeta | theta ~ N(mu2 + c (theta - mu1), v), plus accuracy part
    irt_part(kappa, omega, a, b, prec, num);
    prec += 1.0 + cc * cc / v;
    num += mu1 + (cc / v) * (zeta - mu2 + cc * mu1);
    for (arma::uword j = 0; j < N; ++j) theta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
    // 6.-7. lambda_k and sigma2_k (normal log RT, m = -zeta)
    arma::mat m = arma::repmat(-zeta, 1, K), Wt = arma::repmat((1.0 / sigma2).t(), N, 1);
    lambda = draw_lambda(L, m, arma::zeros(N, K), Wt, mu_l, sd_l, Pr.l_lo);
    sigma2 = draw_sigma2_normal(L, m, lambda, Pr, sigma2);
    // 8. zeta_j: prior N(mu2 + c (theta - mu1), v) + RT part: lambda - log T = zeta + e
    arma::mat yz = -L;
    yz.each_row() += lambda.t();
    prec = 1.0 / v + arma::accu(1.0 / sigma2) * arma::ones(N);
    num = (mu2 + cc * (theta - mu1)) / v + yz * (1.0 / sigma2);
    for (arma::uword j = 0; j < N; ++j) zeta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));
    // 9. location move: common shift of (lambda, zeta); zeta_j has prior N(mu2 + c (theta - mu1), v)
    draw_location_shift(lambda, zeta, mu2 + cc * (theta - mu1), arma::vec(N, arma::fill::value(1.0 / v)), mu_l, sd_l, Pr.l_lo);
    }

    // complete-data log-likelihood (for DIC) and storage
    arma::mat eta = theta * a.t(); eta.each_row() -= (a % b).t();
    arma::mat res = L; res.each_row() -= lambda.t(); res.each_col() += zeta;
    arma::vec d1 = theta - X * beta1, d2 = zeta - X * beta2 - cc * d1;
    loglik(it) = arma::accu(Y % eta - arma::log1p(arma::exp(eta))) +
                 -0.5 * arma::accu(arma::square(res) * (1.0 / sigma2)) - 0.5 * N * arma::accu(arma::log(2 * M_PI * sigma2)) +
                 arma::accu(-0.5 * d1 % d1) - 0.5 * N * std::log(2 * M_PI) +
                 arma::accu(-0.5 * d2 % d2 / v) - 0.5 * N * std::log(2 * M_PI * v);
    arma::vec person = {cc / std::sqrt(cc * cc + v), cc * cc + v};
    draws.row(it) = arma::join_cols(arma::join_cols(a, b, lambda, sigma2), beta1, beta2, person).t();
    if (it >= n_burn) {
      theta_sum += theta; theta_ss += theta % theta; zeta_sum += zeta; zeta_ss += zeta % zeta;
      arma::vec lt = arma::square(res) * (1.0 / sigma2); lt_sum += lt; lt_hit += arma::conv_to<arma::vec>::from(lt > lt_crit);
    }
  }
  double n_keep = n_iter - n_burn;
  return Rcpp::List::create(
    Rcpp::Named("ab_accept") = ab_joint && !one_pl ? ab_accepted / double(n_iter * K) : NA_REAL,
    Rcpp::Named("draws") = draws, Rcpp::Named("loglik") = loglik,
    Rcpp::Named("theta_mean") = theta_sum / n_keep, Rcpp::Named("theta_sd") = arma::sqrt(theta_ss / n_keep - arma::square(theta_sum / n_keep)),
    Rcpp::Named("zeta_mean") = zeta_sum / n_keep, Rcpp::Named("zeta_sd") = arma::sqrt(zeta_ss / n_keep - arma::square(zeta_sum / n_keep)),
    Rcpp::Named("accept") = cov2one ? accepted / (5.0 * n_iter) : NA_REAL,
    Rcpp::Named("lt_mean") = lt_sum / n_keep, Rcpp::Named("lt_exceed") = lt_hit / n_keep,
    Rcpp::Named("last") = Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("zeta") = zeta, Rcpp::Named("a") = a,
                                             Rcpp::Named("b") = b, Rcpp::Named("lambda") = lambda, Rcpp::Named("sigma2") = sigma2,
                                             Rcpp::Named("beta1") = beta1, Rcpp::Named("beta2") = beta2, Rcpp::Named("c") = cc,
                                             Rcpp::Named("v") = v));
}
