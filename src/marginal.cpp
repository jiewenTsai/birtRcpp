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

// mean and second moment of N(m, v) truncated to (lo, hi), given log(Phi(zb) - Phi(za))
static void trunc_moments(double m, double sd, double lo, double hi, double za, double zb, double logZ,
                          double& mean, double& second) {
  if (hi - lo < 1e-10 * sd) {                       // a segment of (almost) zero width
    mean = (std::isfinite(lo) && std::isfinite(hi)) ? 0.5 * (lo + hi) : (std::isfinite(lo) ? lo : hi);
    second = mean * mean;
    return;
  }
  double ra = std::isfinite(za) ? std::exp(R::dnorm(za, 0.0, 1.0, 1) - logZ) : 0.0;
  double rb = std::isfinite(zb) ? std::exp(R::dnorm(zb, 0.0, 1.0, 1) - logZ) : 0.0;
  double zra = std::isfinite(za) ? za * ra : 0.0, zrb = std::isfinite(zb) ? zb * rb : 0.0;
  mean = m + sd * (ra - rb);
  double var = sd * sd * (1.0 + zra - zrb - (ra - rb) * (ra - rb));
  if (!(mean >= lo && mean <= hi) || !(var >= 0.0)) {  // cancellation in a far tail: the segment midpoint
    mean = std::min(std::max(mean, lo), hi);
    var = 0.0;
  }
  second = var + mean * mean;
}

// The integral of .log_ald_rt() and its derivatives (Fisher identity: posterior expectations of
// speed, which is a mixture of truncated normals over the K + 1 segments):
//   d/dB_ik = -(q - P(zeta < -B_ik)) / s_k
//   d/ds_k  = -1 / s_k + E[rho_q(B_ik + zeta)] / s_k^2
//   d/dv    = -1 / (2 v) + E[zeta^2] / (2 v^2)
// plus the posterior mean and second moment of speed.
// [[Rcpp::export(.ald_rt_grad)]]
Rcpp::List ald_rt_grad_cpp(const arma::mat& B, const arma::vec& s, double v, double q) {
  const arma::uword N = B.n_rows, K = B.n_cols;
  const double sd = std::sqrt(v);
  double inv_sum = 0.0, cst = 0.0;
  for (arma::uword k = 0; k < K; ++k) { inv_sum += 1.0 / s(k); cst += std::log(q * (1.0 - q) / s(k)); }
  Rcpp::NumericVector ll(N), gv(N), ez(N), ez2(N);
  arma::mat gb(N, K), gs(N, K);
  std::vector<arma::uword> ord(K);
  std::vector<double> L(K + 1), mean(K + 1), second(K + 1);
  for (arma::uword i = 0; i < N; ++i) {
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](arma::uword x, arma::uword y) { return -B(i, x) < -B(i, y); });
    double bsum = 0.0;
    for (arma::uword k = 0; k < K; ++k) bsum += B(i, k) / s(k);
    double above_inv = inv_sum, above_b = bsum, mx = -INFINITY;
    for (arma::uword j = 0; j <= K; ++j) {
      if (j > 0) { arma::uword k = ord[j - 1]; above_inv -= 1.0 / s(k); above_b -= B(i, k) / s(k); }
      double al = -q * inv_sum + above_inv, be = -q * bsum + above_b, m = al * v;
      double lo = (j == 0) ? -INFINITY : -B(i, ord[j - 1]), hi = (j == K) ? INFINITY : -B(i, ord[j]);
      double za = (lo - m) / sd, zb = (hi - m) / sd;
      if (zb > za) {
        double logZ = log_pdiff(za, zb);
        L[j] = be + al * al * v / 2.0 + logZ;
        trunc_moments(m, sd, lo, hi, za, zb, logZ, mean[j], second[j]);
      } else {
        L[j] = -INFINITY;
        mean[j] = second[j] = 0.0;
      }
      if (L[j] > mx) mx = L[j];
    }
    double acc = 0.0;
    for (arma::uword j = 0; j <= K; ++j) acc += std::exp(L[j] - mx);
    ll[i] = cst + mx + std::log(acc);
    // segment probabilities, cumulated from below: item ord[j] has its kink at the top of segment j
    double cum_p = 0.0, cum_m = 0.0, e1 = 0.0, e2 = 0.0;
    std::vector<double> below_p(K), below_m(K);
    for (arma::uword j = 0; j <= K; ++j) {
      double pj = std::exp(L[j] - mx) / acc;
      cum_p += pj;
      cum_m += pj * mean[j];
      e1 += pj * mean[j];
      e2 += pj * second[j];
      if (j < K) { below_p[ord[j]] = cum_p; below_m[ord[j]] = cum_m; }
    }
    ez[i] = e1;
    ez2[i] = e2;
    gv[i] = -0.5 / v + e2 / (2.0 * v * v);
    for (arma::uword k = 0; k < K; ++k) {
      double b = B(i, k);
      double check = q * (b + e1) - b * below_p[k] - below_m[k];   // E[rho_q(b + zeta)]
      gb(i, k) = -(q - below_p[k]) / s(k);
      gs(i, k) = -1.0 / s(k) + check / (s(k) * s(k));
    }
  }
  return Rcpp::List::create(Rcpp::Named("ll") = ll, Rcpp::Named("gb") = gb, Rcpp::Named("gs") = gs,
                            Rcpp::Named("gv") = gv, Rcpp::Named("ez") = ez, Rcpp::Named("ez2") = ez2);
}
