// [[Rcpp::depends(RcppArmadillo)]]
#include "rng.h"
#include <cmath>

// ---- Polya-Gamma PG(1, z) ---------------------------------------------------------------------
// PG(1, z) = J*(1, z / 2) / 4, where J* is sampled by Devroye's method: a proposal mixing a
// truncated exponential (x > t) and a truncated inverse Gaussian (x < t), accepted by the
// alternating series of the J* density (Polson, Scott & Windle, 2013, Section 4).
static const double PG_T = 0.64;

static double pg_coef(int n, double x) {            // n-th term of the series for the J* density
  double k = (n + 0.5) * M_PI;
  if (x > PG_T) return k * std::exp(-0.5 * k * k * x);
  if (x <= 0) return 0.0;
  return std::exp(-1.5 * (std::log(0.5 * M_PI) + std::log(x)) + std::log(k) - 2.0 * (n + 0.5) * (n + 0.5) / x);
}

static double pg_mass_exp(double z) {               // probability of the exponential piece
  double fz = 0.125 * M_PI * M_PI + 0.5 * z * z;
  double b = std::sqrt(1.0 / PG_T) * (PG_T * z - 1.0);
  double a = -std::sqrt(1.0 / PG_T) * (PG_T * z + 1.0);
  double x0 = std::log(fz) + fz * PG_T;
  double xb = x0 - z + R::pnorm(b, 0.0, 1.0, 1, 1);
  double xa = x0 + z + R::pnorm(a, 0.0, 1.0, 1, 1);
  return 1.0 / (1.0 + 4.0 / M_PI * (std::exp(xb) + std::exp(xa)));
}

static double pg_tinvgauss(double z) {              // inverse Gaussian(1 / z, 1) truncated to (0, t)
  double mu = 1.0 / z, x = PG_T + 1.0;
  if (mu > PG_T) {
    double alpha = 0.0;
    while (R::unif_rand() > alpha) {
      double e1 = R::exp_rand(), e2 = R::exp_rand();
      while (e1 * e1 > 2.0 * e2 / PG_T) { e1 = R::exp_rand(); e2 = R::exp_rand(); }
      x = PG_T / ((1.0 + PG_T * e1) * (1.0 + PG_T * e1));
      alpha = std::exp(-0.5 * z * z * x);
    }
  } else {
    while (x > PG_T) {
      double y = R::norm_rand(); y *= y;
      x = mu + 0.5 * mu * mu * y - 0.5 * mu * std::sqrt(4.0 * mu * y + (mu * y) * (mu * y));
      if (R::unif_rand() > mu / (mu + x)) x = mu * mu / x;
    }
  }
  return x;
}

double rpg1(double z) {
  if (!std::isfinite(z))
    Rcpp::stop("non-finite linear predictor in the Polya-Gamma step: a parameter diverged (e.g. a discrimination, "
               "a variance going to 0 or a person value running off); check the data (items without variation, "
               "extreme times) or simplify the model");
  z = std::fabs(z) * 0.5;
  double fz = 0.125 * M_PI * M_PI + 0.5 * z * z;
  for (int tries = 0; tries < 100000; ++tries) {
    double x = (R::unif_rand() < pg_mass_exp(z)) ? PG_T + R::exp_rand() / fz : pg_tinvgauss(z);
    double s = pg_coef(0, x), y = R::unif_rand() * s;
    // the alternating series decides within a few terms; the cap is a safety net
    for (int n = 1; n < 1000; ++n) {
      if (n % 2 == 1) { s -= pg_coef(n, x); if (y <= s) return 0.25 * x; }
      else            { s += pg_coef(n, x); if (y > s) break; }
    }
  }
  Rcpp::stop("Polya-Gamma sampler did not accept a draw (z = %f)", 2.0 * z);
}

// ---- truncated normal by inversion, on the log scale in the tails ------------------------------
double rtnorm(double mu, double sd, double lo, double hi) {
  double a = (lo - mu) / sd, b = (hi - mu) / sd, u = R::unif_rand(), x;
  if (a >= 0) {                                      // both bounds in the upper tail: invert Q = 1 - Phi
    double la = R::pnorm(a, 0.0, 1.0, 0, 1), lb = R::pnorm(b, 0.0, 1.0, 0, 1);
    x = R::qnorm(la + std::log1p(-u * (1.0 - std::exp(lb - la))), 0.0, 1.0, 0, 1);
  } else if (b <= 0) {                               // both bounds in the lower tail: invert Phi
    double la = R::pnorm(a, 0.0, 1.0, 1, 1), lb = R::pnorm(b, 0.0, 1.0, 1, 1);
    x = R::qnorm(lb + std::log1p(-u * (1.0 - std::exp(la - lb))), 0.0, 1.0, 1, 1);
  } else {
    double pa = R::pnorm(a, 0.0, 1.0, 1, 0), pb = R::pnorm(b, 0.0, 1.0, 1, 0);
    x = R::qnorm(pa + u * (pb - pa), 0.0, 1.0, 1, 0);
  }
  x = std::min(std::max(x, a), b);                   // guard against rounding at the bounds
  return mu + sd * x;
}

// ---- inverse Gaussian ---------------------------------------------------------------------------
double rinvgauss(double mu, double lambda) {
  double y = R::norm_rand(); y *= y;
  double x = mu + 0.5 * mu * mu * y / lambda - 0.5 * mu / lambda * std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  return (R::unif_rand() <= mu / (mu + x)) ? x : mu * mu / x;
}

double rinvgamma(double shape, double scale) { return 1.0 / R::rgamma(shape, 1.0 / scale); }

arma::vec rmvn_prec(const arma::mat& Q, const arma::vec& b) {
  arma::mat U = arma::chol(arma::symmatu(Q));        // Q = U' U
  arma::vec mean = arma::solve(arma::trimatu(U), arma::solve(arma::trimatl(U.t()), b));
  arma::vec z(b.n_elem);
  for (arma::uword i = 0; i < b.n_elem; ++i) z(i) = R::norm_rand();
  return mean + arma::solve(arma::trimatu(U), z);    // Cov = U^-1 U^-T = Q^-1
}

// exported for testing
// [[Rcpp::export(.rpg1)]]
Rcpp::NumericVector rpg1_r(Rcpp::NumericVector z) {
  Rcpp::NumericVector out(z.size());
  for (R_xlen_t i = 0; i < z.size(); ++i) out[i] = rpg1(z[i]);
  return out;
}
// [[Rcpp::export(.rtnorm)]]
Rcpp::NumericVector rtnorm_r(int n, double mu, double sd, double lo, double hi) {
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = rtnorm(mu, sd, lo, hi);
  return out;
}
// [[Rcpp::export(.rinvgauss)]]
Rcpp::NumericVector rinvgauss_r(int n, double mu, double lambda) {
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = rinvgauss(mu, lambda);
  return out;
}
