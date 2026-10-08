# Example 1. From a data frame to the tables of a cross-relation model (dissertation Chapter 5).
# ECM first (estimates, standard errors, AIC / BIC in about a second); then Gibbs from the ECM solution.
library(birtRcpp)
set.seed(1)
cond <- set_cond(n_subj = 500, n_item = 10)
sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
raw <- data.frame(id = sprintf("P%03d", 1:500), sim$Y, sim$T)             # a data frame as it comes
names(raw) <- c("id", paste0("q", 1:10), paste0("q", 1:10, "_sec"))
dat <- input_data(resp = paste0("q", 1:10), time = paste0("q", 1:10, "_sec"), id = "id", data = raw)

# 1. ECM: maximum marginal likelihood
fit <- ecm(rtirt_cross(dat))
fit                                  # logLik, npar, AIC, BIC, reliability
s <- summary(fit)                    # all tables, computed once
s                                    # $fit, $items, $se, $parameters, $reliability
s$items                              # one row per item: a, b, lambda, sigma2t, rho
s$se                                 # their standard errors
print(s, tables = c("fit", "parameters"))   # estimate, SE, z, p and Wald 95% interval
head(estimates(fit, pars = "^rho"))  # long table: est, se, z, q025, q975
head(scores(fit))                    # EAP and posterior SD of ability and speed, with the id
c(logLik = logLik(fit), AIC = AIC(fit), BIC = BIC(fit))
c(iterations = fit$ecm$iterations, converged = fit$ecm$converged)
equations(fit)                       # the model as equations (Unicode)

# 2. Gibbs started at the ECM solution: posterior intervals, R-hat, DIC / WAIC / LOOIC
g <- gibbs(fit, n_iter = 2000, n_chain = 3, init = "ecm", seed = 1)
g
head(as.data.frame(summary(g)))      # the long parameter table, with the priors
convergence(g)
cbind(ecm = estimates(fit, pars = "^rho")$est, se = estimates(fit, pars = "^rho")$se,
      gibbs = estimates(g, pars = "^rho")$est, post_sd = estimates(g, pars = "^rho")$sd)
