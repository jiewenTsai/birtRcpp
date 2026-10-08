# Example 3. Quantile models: one model per quantile level (ALD working likelihood).
# Gibbs only: ecm() fits the normal (mean) models, which serve as the reference.
library(birtRcpp)
set.seed(3)
cond <- set_cond(n_subj = 500, n_item = 8, n_feat = 2)
d_cross <- sim_data(cond, sim_para(cond, "cross"), "cross", type = "skew")
mean_fit <- ecm(rtirt_cross(d_cross))                    # the mean model by ECM
try(ecm(rtirt_cross(d_cross, quantile = 0.5)))           # refused: quantile models use gibbs()
qs <- gibbs(rtirt_cross(d_cross, quantile = c(0.1, 0.25, 0.5, 0.75, 0.9)), n_iter = 2000, n_chain = 2, seed = 1)
qs                                                       # one line per quantile level
# q = 0.9 on this right-skewed data set: the ALD working posterior is multimodal (R-hat warning). Ability becomes
# a second speed factor (a about 0.15, |rho| about 0.55 with item-specific signs) and each chain settles on its own
# sign pattern, about equally likely; chains started at the ECM fit of the normal model do the same. Not a sampler
# fault: check convergence() and do not read the extreme levels (rho above all) on skewed data.
rho_q <- coef(qs)[grepl("^rho", rownames(coef(qs))), ]   # rho_i(q): items x quantiles
round(cbind(mean_ecm = estimates(mean_fit, pars = "^rho")$est, rho_q), 3)
plot(qs)                                                 # rho_i against q
d_lat <- sim_data(cond, sim_para(cond, "latent"), "latent")
ql <- gibbs(rtirt_latent(d_lat, quantile = 0.25), n_iter = 2000, n_chain = 2, seed = 1)
estimates(ql, pars = "beta|b_ability")                   # the 0.25-quantile regression of speed
# collapse = TRUE: speed integrated out of the regression, ability and RT item updates; same posterior,
# several times the effective sample size of the regression in quantile models
qc <- gibbs(rtirt_latent(d_lat, quantile = 0.25), n_iter = 2000, n_chain = 2, seed = 1, collapse = TRUE)
cbind(original = estimates(ql, pars = "beta|b_ability")$ess, collapsed = estimates(qc, pars = "beta|b_ability")$ess)
equations(ql)
