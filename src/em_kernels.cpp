// Kernels of the EM algorithm for nodes that differ between persons (latent regressions: theta_jq =
// x_j' beta + z_q). With nodes common to all persons the same quantities are matrix products in R
// (calc_loglik_nodes(), calc_artificial_data()). Persons j are rows, items i columns.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "rng.h"

static inline double softplus(double x) { return x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)); }

// N x Q log-likelihood of the 2PL responses at the nodes TH (N x Q): sum_i y eta - log(1 + e^eta),
// eta = a_i (theta - b_i); missing responses (NaN) are skipped
// [[Rcpp::export(.em_loglik_nodes)]]
arma::mat em_loglik_nodes(const arma::mat& Y, const arma::vec& a, const arma::vec& b, const arma::mat& TH) {
  const arma::uword N = Y.n_rows, K = Y.n_cols, Q = TH.n_cols;
  arma::mat L(N, Q, arma::fill::zeros);
  for (arma::uword q = 0; q < Q; ++q)
    for (arma::uword i = 0; i < K; ++i)
      for (arma::uword j = 0; j < N; ++j) {
        const double y = Y(j, i); if (std::isnan(y)) continue;
        const double eta = a(i) * (TH(j, q) - b(i));
        L(j, q) += y * eta - softplus(eta);
      }
  return L;
}

// Weighted item sums over persons and nodes (weights W, N x Q): with score = false the Polya-Gamma
// minorizer sums S_r = sum W E[omega | eta] theta^r (r = 0, 1, 2), with score = true the score sums
// S_r = sum W (y - P) theta^r (r = 0, 1)
// [[Rcpp::export(.em_item_sums)]]
Rcpp::List em_item_sums(const arma::mat& Y, const arma::vec& a, const arma::vec& b, const arma::mat& TH, const arma::mat& W, bool score) {
  const arma::uword N = Y.n_rows, K = Y.n_cols, Q = TH.n_cols;
  arma::vec S0(K, arma::fill::zeros), S1(K, arma::fill::zeros), S2(K, arma::fill::zeros);
  for (arma::uword q = 0; q < Q; ++q)
    for (arma::uword i = 0; i < K; ++i)
      for (arma::uword j = 0; j < N; ++j) {
        const double y = Y(j, i); if (std::isnan(y)) continue;
        const double th = TH(j, q), eta = a(i) * (th - b(i));
        double g;
        if (score) g = y - 1.0 / (1.0 + std::exp(-eta));
        else g = pg_mean1(eta);
        g *= W(j, q);
        S0(i) += g; S1(i) += th * g; if (!score) S2(i) += th * th * g;
      }
  return Rcpp::List::create(Rcpp::Named("S0") = S0, Rcpp::Named("S1") = S1, Rcpp::Named("S2") = S2);
}

// E[omega], omega ~ PG(1, eta), elementwise (keeps the dimensions of eta)
// [[Rcpp::export(pg_mean)]]
Rcpp::NumericVector pg_mean_vec(const Rcpp::NumericVector& eta) {
  Rcpp::NumericVector out = Rcpp::clone(eta);
  for (R_xlen_t k = 0; k < out.size(); ++k) out[k] = pg_mean1(eta[k]);
  return out;
}
