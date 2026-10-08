// Closed-form integral over speed for rtirt_cross quantile models (see R/marginal.R):
// log int prod_k ALD(b_ik + zeta; q, s_k) N(zeta; 0, v) dzeta for each row i of B.
// The ALD log-likelihood is piecewise linear in zeta with kinks at -b_ik, so the integral is
// a sum over K + 1 segments of exp(alpha zeta + beta) N(zeta; 0, v), i.e. normal-CDF differences.
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// log(Phi(b) - Phi(a)) for a < b, stable in both tails
static double log_pdiff(double a, double b) {
  if (a > 0) {                                     // upper tail: use survival functions
    double la = R::pnorm(a, 0.0, 1.0, 0, 1), lb = R::pnorm(b, 0.0, 1.0, 0, 1);
    return la + std::log1p(-std::exp(lb - la));
  }
  double la = R::pnorm(a, 0.0, 1.0, 1, 1), lb = R::pnorm(b, 0.0, 1.0, 1, 1);
  return lb + std::log1p(-std::exp(la - lb));
}

// [[Rcpp::export(.log_ald_rt)]]
Rcpp::NumericVector log_ald_rt_cpp(const arma::mat& B, const arma::vec& s, double v, double q) {
  const arma::uword N = B.n_rows, K = B.n_cols;
  const double sd = std::sqrt(v);
  double inv_sum = 0.0, cst = 0.0;
  for (arma::uword k = 0; k < K; ++k) { inv_sum += 1.0 / s(k); cst += std::log(q * (1.0 - q) / s(k)); }
  Rcpp::NumericVector out(N);
  std::vector<arma::uword> ord(K);
  std::vector<double> L(K + 1);
  for (arma::uword i = 0; i < N; ++i) {
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](arma::uword x, arma::uword y) { return -B(i, x) < -B(i, y); });
    double bsum = 0.0;
    for (arma::uword k = 0; k < K; ++k) bsum += B(i, k) / s(k);
    // segment j: zeta in (kink_(j), kink_(j+1)); items with kinks above zeta (sorted positions >= j) have u < 0
    double above_inv = inv_sum, above_b = bsum, mx = -INFINITY;
    for (arma::uword j = 0; j <= K; ++j) {
      if (j > 0) { arma::uword k = ord[j - 1]; above_inv -= 1.0 / s(k); above_b -= B(i, k) / s(k); }
      double al = -q * inv_sum + above_inv, be = -q * bsum + above_b, m = al * v;
      double lo = (j == 0) ? -INFINITY : -B(i, ord[j - 1]), hi = (j == K) ? INFINITY : -B(i, ord[j]);
      double za = (lo - m) / sd, zb = (hi - m) / sd;
      L[j] = (zb > za) ? be + al * al * v / 2.0 + log_pdiff(za, zb) : -INFINITY;
      if (L[j] > mx) mx = L[j];
    }
    double acc = 0.0;
    for (arma::uword j = 0; j <= K; ++j) acc += std::exp(L[j] - mx);
    out[i] = cst + mx + std::log(acc);
  }
  return out;
}
