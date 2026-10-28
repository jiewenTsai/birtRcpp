# Example 4. ECM: maximum marginal likelihood (Meng & Rubin, 1993), MAP for small samples,
# and Gibbs chains started from the ECM solution
library(birtRcpp)
set.seed(4)
cond <- set_cond(n_subj = 800, n_item = 10, n_feat = 2)
dat <- sim_data(cond, sim_para(cond, "latreg", true_corr = -0.3), "latreg")
e <- ecm(dat, model = "latreg") # about a second, with standard errors (observed information)
e
summary(e) # the same tables as for Gibbs: est, SE, z, Wald intervals
logLik(e)
AIC(e)
BIC(e)
c(iterations = e$ecm$iterations, converged = e$ecm$converged)
plot(e$ecm$trace, type = "b", xlab = "iteration", ylab = "log-likelihood") # never decreases
# quadrature check: about 0
e61 <- ecm(dat, model = "latreg", se = FALSE, control = ecm_control(nodes = 61))
e61$ecm$logLik - e$ecm$logLik
# without SQUAREM
e_plain <- ecm(dat, model = "latreg", se = FALSE, control = ecm_control(accelerate = FALSE))
e_plain$ecm$iterations
# Gibbs from the ECM solution: a short burn-in suffices
g <- gibbs(dat, model = "latreg", init = e, n_iter = 1000, n_chain = 2, seed = 1)
max(estimates(g)$rhat, na.rm = TRUE)
est_ecm <- estimates(e, pars = "beta|cor")
est_gibbs <- estimates(g, pars = "beta|cor")
tibble::tibble(ecm = est_ecm$est, se = est_ecm$se, gibbs = est_gibbs$est, post_sd = est_gibbs$sd)
# small samples: prior = "map" adds the item priors of gibbs() (posterior mode in the item parameters)
cs <- set_cond(n_subj = 100, n_item = 10)
small <- sim_data(cs, sim_para(cs, "null"), "null")
tight <- rtirt_priors(a = c(mean = 1, sd = 0.3))
ml <- ecm(small, model = "null", se = FALSE)
map <- ecm(small, model = "null", prior = "map", se = FALSE)
map_tight <- ecm(small, model = "null", prior = "map", priors = tight, se = FALSE)
discrimination <- function(fit) estimates(fit, pars = "^a\\[")$est
a_table <- rbind(
  true = attr(small, "true_para")$a, ml = discrimination(ml),
  map = discrimination(map), map_tight = discrimination(map_tight)
)
round(a_table, 2)
summary(ecm(small, model = "null", prior = "map"))
