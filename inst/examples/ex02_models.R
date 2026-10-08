# Example 2. The five models on simulated data, fitted by ECM (maximum marginal likelihood)
library(birtRcpp)
set.seed(2)
cond <- set_cond(n_subj = 400, n_item = 8, n_feat = 2)
models <- c(mlirt = "mlirt", null = "null", latreg = "latreg", latent = "latent", cross = "cross")
sims <- lapply(models, function(m) sim_data(cond, sim_para(cond, m), m))
fits <- list(
  mlirt  = ecm(mlirt(sims$mlirt)),
  null   = ecm(rtirt_null(sims$null)),
  latreg = ecm(rtirt_latreg(sims$latreg)),
  latent = ecm(rtirt_latent(sims$latent)),
  cross  = ecm(rtirt_cross(sims$cross)))
sapply(fits, function(f) c(iterations = f$ecm$iterations, npar = f$ecm$npar, secs = round(f$ecm$secs, 2)))
lapply(fits, reliability)
estimates(fits$latreg, pars = "beta|cor|var_speed")     # structural parameters of the joint model
estimates(fits$latent, pars = "beta|b_ability")         # speed regressed on ability and covariates
rho <- estimates(fits$cross, pars = "^rho")             # recovery of the cross-relations
cbind(true = attr(sims$cross, "true_para")$rho, est = rho$est, se = rho$se)
# options: fixed speed variance, Rasch-type items
f2 <- ecm(rtirt_latreg(sims$latreg, speed_var = "fixed", itemtype = "1pl"))
f2
equations(f2)
# the same cross model by Gibbs, started at the ECM solution; compare_para() works on Gibbs fits
g <- gibbs(fits$cross, n_iter = 1500, n_chain = 2, init = "ecm", seed = 1)
compare_para(g, "rho")
