// Full conditionals shared by the models. Notation (persons j = 1..N, items k = 1..K):
//   accuracy  logit P(y_jk = 1) = a_k (theta_j - b_k), Polya-Gamma variables omega_jk, kappa = y - 1/2
//   log RT    L_jk = lambda_k + m_jk + k1 nu_jk + e_jk,  W_jk = 1 / Var(e_jk)
//             (m_jk = -zeta_j, or -zeta_j - rho_k theta_j; nu, k1 only for the ALD model)
// Priors: see struct Prior (all full conditionals are conjugate; the half-t through an inverse-gamma
//         auxiliary variable).
#ifndef BIRTRCPP_BLOCKS_H
#define BIRTRCPP_BLOCKS_H

#include <RcppArmadillo.h>
#include <functional>

struct Ald { double k1, k2; };

// Priors (rtirt_priors() in R, passed as a numeric vector in this order):
//   a ~ N(a_mu, a_sd^2) on (0, Inf)
//   b ~ N(b_mu, b_sd^2) on [b_lo, b_hi]; lambda ~ N(l_mu, l_sd^2) on (l_lo, Inf) (l_mu NaN: mean and sd of log T)
//   sigma2t and the speed (residual) variance: fam 0 = IG(p1 = shape, p2 = scale) on the variance,
//     1 = half-t(p1 = df, p2 = scale) on its square root
//   regression coefficients ~ N(0, beta_sd^2), rho ~ N(0, rho_sd^2), covariance c ~ N(0, c_sd^2)
struct Prior { double a_mu, a_sd, b_mu, b_sd, b_lo, b_hi, l_mu, l_sd, l_lo, s2_fam, s2_p1, s2_p2,
               beta_sd, rho_sd, vs_fam, vs_p1, vs_p2, c_sd; };
Prior make_prior(const Rcpp::NumericVector& v);
Prior julia_prior();                                   // the priors of ExtendedRtIrtModeling.jl (tests)

// a variance with likelihood var^-alpha exp(-beta / var) and the prior (fam, p1, p2) above: conjugate
// inverse gamma, or for the half-t two inverse-gamma Gibbs steps through the auxiliary variable of
// Huang & Wand (2013), drawn from the current value of the variance
double draw_var(double alpha, double beta, double fam, double p1, double p2, double current);
Ald ald_constants(double q);                         // k1 = (1-2q)/(q(1-q)), k2 = 2/(q(1-q))

// ---- 2PL accuracy model -------------------------------------------------------------------
arma::mat draw_omega(const arma::vec& a, const arma::vec& b, const arma::vec& theta);
arma::vec draw_a(const arma::mat& kappa, const arma::mat& omega, const arma::vec& theta, const arma::vec& b, const Prior& P);
arma::vec draw_b(const arma::mat& kappa, const arma::mat& omega, const arma::vec& theta, const arma::vec& a, const Prior& P);
// accuracy part of theta's full conditional: precision and precision-weighted mean
void irt_part(const arma::mat& kappa, const arma::mat& omega, const arma::vec& a, const arma::vec& b,
              arma::vec& prec, arma::vec& num);

// ---- log RT model ---------------------------------------------------------------------------
arma::vec draw_lambda(const arma::mat& L, const arma::mat& m, const arma::mat& k1nu, const arma::mat& W,
                      double mu_l, double sd_l, double lo_l = 0.0);
arma::vec draw_sigma2_normal(const arma::mat& L, const arma::mat& m, const arma::vec& lambda, const Prior& P, const arma::vec& current);
arma::vec draw_sigma2_ald(const arma::mat& L, const arma::mat& m, const arma::vec& lambda,
                          const arma::mat& nu, const Ald& c, const Prior& P, const arma::vec& current);
// Location move: log T depends on lambda_k - zeta_j only, so (lambda + c, zeta + c) leaves the
// likelihood unchanged and the common shift c is identified by the priors alone. Plain Gibbs
// moves along this direction slowly; this step draws c from its full conditional (a Gibbs step
// on the translation group, Liu & Sabatti, 2000), which leaves the posterior invariant.
// zeta_j has prior N(mz_j, 1 / wz_j). Returns c (0 if a shifted lambda would leave (0, Inf)).
double draw_location_shift(arma::vec& lambda, arma::vec& zeta, const arma::vec& mz, const arma::vec& wz,
                           double mu_l, double sd_l, double lo_l = 0.0);

// ALD mixing variable for residual r with scale s: GIG(1/2, r^2/(k2 s), (k1^2 + 2 k2)/(k2 s))
double draw_nu(double r, double s, const Ald& c);

// ---- collapsed (partially collapsed) steps, used when gibbs(collapse = TRUE) -----------------
// Given the Polya-Gamma (and ALD) variables, the model is linear-Gaussian in the person and RT item
// parameters, so speed zeta can be integrated out of their updates in closed form. The sampler then
// draws theta and the RT item parameters from these marginals and draws zeta last, from its full
// conditional (a partially collapsed Gibbs sampler, van Dyk & Park, 2008).
// RT part of person j: y_jk = zeta_j + rho_k theta_j + error, weights W_jk (y = lambda + k1 nu - log T),
// and zeta_j | theta_j ~ N(m0_j + g_j theta_j, tau_j). prec / num: the rest of theta's conditional
// (accuracy part and theta's own prior).
arma::vec draw_theta_collapsed(const arma::vec& prec, const arma::vec& num, const arma::mat& W, const arma::mat& y,
                               const arma::vec& rho, const arma::vec& m0, const arma::vec& g, const arma::vec& tau);
// (lambda, rho) jointly with zeta integrated out: z_j = log T_j - k1 nu_j + mz_j = lambda - rho theta_j + u_j,
// Cov(u_j) = diag(1 / W_j) + tau_j 1 1', mz_j = E(zeta_j | theta_j). Priors lambda_k ~ N+(mu_l, sd_l^2),
// rho_k ~ N(0, rho_sd^2). with_rho = false: lambda only (rho stays 0).
void draw_rt_items_collapsed(const arma::mat& L, const arma::mat& k1nu, const arma::mat& W, const arma::vec& theta,
                             const arma::vec& mz, const arma::vec& tau, bool with_rho, double mu_l, double sd_l,
                             double rho_sd, arma::vec& lambda, arma::vec& rho, double lo_l = 0.0);

// Fisher-scoring Metropolis-Hastings for the two parameters of one item (blocks.cpp): the prior adds its
// log density, gradient and metric (terms) and gives its support (proposals outside are rejected).
struct ItemPrior {
  std::function<bool(const arma::vec&)> support;
  std::function<void(const arma::vec&, double&, arma::vec&, arma::mat&)> terms;
};
bool item_mh_step(const arma::vec& y, const arma::vec& theta, arma::vec& x, int link, int form, const ItemPrior& prior);

// Joint (a_k, b_k) update with the Polya-Gamma variables integrated out: item_mh_step() on the exact
// logistic posterior given theta (priors a ~ N+, b normal, see Prior), in the coordinates (a, d = -a b),
// where the posterior is near-normal. For a near-normal posterior the Fisher-scoring proposal draws from
// the posterior itself and nearly always accepts. Far in the tails the reverse proposal density can be
// tiny and the step may stall, so the samplers follow it with the Polya-Gamma steps. Returns the number
// of accepted items.
int draw_ab_joint(const arma::mat& Y, const arma::vec& theta, arma::vec& a, arma::vec& b, const Prior& P);

#endif
