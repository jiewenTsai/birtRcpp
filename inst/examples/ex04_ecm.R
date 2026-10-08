# Example 4. ECM: maximum marginal likelihood (Meng & Rubin, 1993), MAP for small samples,
# and Gibbs chains started from the ECM solution
library(birtRcpp)
set.seed(4)
cond <- set_cond(n_subj = 800, n_item = 10, n_feat = 2)
dat <- sim_data(cond, sim_para(cond, "latreg", true_corr = -0.3), "latreg")
e <- ecm(rtirt_latreg(dat))          # about a second, with standard errors (observed information)
e
summary(e)                           # the same tables as for Gibbs: est, SE, z, Wald intervals
logLik(e); AIC(e); BIC(e)
c(iterations = e$ecm$iterations, converged = e$ecm$converged)
plot(e$ecm$trace, type = "b", xlab = "iteration", ylab = "log-likelihood")   # never decreases
ecm(rtirt_latreg(dat), nodes = 61, se = FALSE)$ecm$logLik - e$ecm$logLik      # quadrature check: ~0
ecm(rtirt_latreg(dat), accelerate = FALSE, se = FALSE)$ecm$iterations         # without SQUAREM
# Gibbs from the ECM solution: a short burn-in suffices
g <- gibbs(e, n_iter = 1000, n_chain = 2, init = "ecm", seed = 1)
max(estimates(g)$rhat, na.rm = TRUE)
cbind(ecm = estimates(e, pars = "beta|cor")$est, se = estimates(e, pars = "beta|cor")$se,
      gibbs = estimates(g, pars = "beta|cor")$est, post_sd = estimates(g, pars = "beta|cor")$sd)
# small samples: prior = "map" adds the item priors of gibbs() (posterior mode in the item parameters)
cs <- set_cond(n_subj = 100, n_item = 10)
small <- sim_data(cs, sim_para(cs, "null"), "null")
ml <- ecm(rtirt_null(small), se = FALSE)
map <- ecm(rtirt_null(small), prior = "map", se = FALSE)
map_tight <- ecm(rtirt_null(small), prior = "map", priors = rtirt_priors(a = c(mean = 1, sd = 0.3)), se = FALSE)
round(rbind(true = attr(small, "true_para")$a, ml = estimates(ml, pars = "^a\\[")$est,
            map = estimates(map, pars = "^a\\[")$est, map_tight = estimates(map_tight, pars = "^a\\[")$est), 2)
summary(ecm(rtirt_null(small), prior = "map"))
