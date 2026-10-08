// 2PL IRT with a latent regression of ability (no response times).
//   logit P(y_jk = 1) = a_k (theta_j - b_k),   theta_j = x_j' beta + e_j, e ~ N(0, 1),   beta ~ N(0, I)
// collapse = true: beta is drawn with theta integrated out (given the Polya-Gamma variables, the
// accuracy part of theta_j is a normal pseudo-observation num_j / prec_j with variance 1 / prec_j).
// [[Rcpp::depends(RcppArmadillo)]]
#include "blocks.h"
#include "rng.h"

// [[Rcpp::export(.gibbs_mlirt)]]
Rcpp::List gibbs_mlirt(const arma::mat& Y, const arma::mat& X, bool one_pl, bool collapse, bool ab_joint, int n_iter, int n_burn, Rcpp::List init,
                       Rcpp::NumericVector prior) {
  const Prior Pr = make_prior(prior);
  const arma::uword N = Y.n_rows, K = Y.n_cols, P = X.n_cols;
  const arma::mat kappa = Y - 0.5;
  arma::vec theta = init["theta"], a = init["a"], b = init["b"], beta = init["beta"], prec, num;
  arma::mat omega(N, K), draws(n_iter, 2 * K + P);
  arma::vec loglik(n_iter), theta_sum(N, arma::fill::zeros), theta_ss(N, arma::fill::zeros);

  int ab_accepted = 0;
  for (int it = 0; it < n_iter; ++it) {
    Rcpp::checkUserInterrupt();
    // 1. beta: theta - X beta ~ N(0, 1), prior N(0, I)
    if (P > 0 && !collapse) beta = rmvn_prec(arma::eye(P, P) / (Pr.beta_sd * Pr.beta_sd) + X.t() * X, X.t() * theta);
    // 2. accuracy: Polya-Gamma variables, difficulty, discrimination
    // collapse: first (a, b) jointly with omega integrated out (moves along the a-b ridge), then the
    // Polya-Gamma steps as usual (they always move, so the chain cannot stall in the tails)
    if (ab_joint && !one_pl) ab_accepted += draw_ab_joint(Y, theta, a, b, Pr);
    omega = draw_omega(a, b, theta);
    b = draw_b(kappa, omega, theta, a, Pr);
    if (!one_pl) a = draw_a(kappa, omega, theta, b, Pr);
    irt_part(kappa, omega, a, b, prec, num);
    // 1c. beta, theta integrated out: num / prec ~ N(x' beta, 1 + 1 / prec)
    if (P > 0 && collapse) {
      arma::vec wt = prec / (1.0 + prec);
      beta = rmvn_prec(arma::eye(P, P) / (Pr.beta_sd * Pr.beta_sd) + X.t() * (X.each_col() % wt), X.t() * (wt % (num / prec)));
    }
    // 3. theta_j: prior N(x' beta, 1) + accuracy part
    arma::vec mu = X * beta;
    prec += 1.0; num += mu;
    for (arma::uword j = 0; j < N; ++j) theta(j) = num(j) / prec(j) + R::norm_rand() / std::sqrt(prec(j));

    arma::mat eta = theta * a.t(); eta.each_row() -= (a % b).t();
    loglik(it) = arma::accu(Y % eta - arma::log1p(arma::exp(eta))) + arma::accu(-0.5 * arma::square(theta - mu)) - 0.5 * N * std::log(2 * M_PI);
    draws.row(it) = arma::join_cols(a, b, beta).t();
    if (it >= n_burn) { theta_sum += theta; theta_ss += theta % theta; }
  }
  double n_keep = n_iter - n_burn;
  return Rcpp::List::create(
    Rcpp::Named("ab_accept") = ab_joint && !one_pl ? ab_accepted / double(n_iter * K) : NA_REAL,
    Rcpp::Named("draws") = draws, Rcpp::Named("loglik") = loglik,
    Rcpp::Named("theta_mean") = theta_sum / n_keep,
    Rcpp::Named("theta_sd") = arma::sqrt(theta_ss / n_keep - arma::square(theta_sum / n_keep)),
    Rcpp::Named("last") = Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("a") = a, Rcpp::Named("b") = b,
                                             Rcpp::Named("beta") = beta));
}
