// Random number generators used by the Gibbs samplers. All draws use R's RNG, so set.seed()
// makes a run reproducible.
#ifndef BIRTRCPP_RNG_H
#define BIRTRCPP_RNG_H

#include <RcppArmadillo.h>

// Polya-Gamma PG(1, z) (Polson, Scott & Windle, 2013; Devroye's alternating-series sampler)
double rpg1(double z);
// its mean E[omega] = tanh(z / 2) / (2 z) (Polson, Scott & Windle, 2013, Sec. 2.2); the series 1/4 - z^2/48 near 0
inline double pg_mean1(double z) { return std::fabs(z) < 1e-4 ? 0.25 - z * z / 48.0 : std::tanh(z / 2.0) / (2.0 * z); }

// N(mu, sd^2) truncated to (lo, hi); lo / hi may be -Inf / Inf
double rtnorm(double mu, double sd, double lo, double hi);

// inverse Gaussian with mean mu and shape lambda (Michael, Schucany & Haas, 1976)
double rinvgauss(double mu, double lambda);

// inverse gamma with shape and scale: 1 / Gamma(shape, rate = scale)
double rinvgamma(double shape, double scale);

// multivariate normal from a precision matrix Q and a vector b: N(Q^-1 b, Q^-1)
arma::vec rmvn_prec(const arma::mat& Q, const arma::vec& b);

#endif
