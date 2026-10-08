# Example 5. Comparing latent structures: AIC / BIC / LR test (ECM) and marginal LOOIC (Gibbs)
library(birtRcpp)
set.seed(5)
cond <- set_cond(n_subj = 600, n_item = 10, n_feat = 2)
dat <- sim_data(cond, sim_para(cond, "latent"), "latent")     # speed depends on ability and x
m <- list(null = rtirt_null(dat), latreg = rtirt_latreg(dat), latent = rtirt_latent(dat))
# maximum likelihood by ECM
e <- lapply(m, ecm, se = FALSE)
data.frame(model = names(e), logLik = sapply(e, logLik), npar = sapply(e, function(f) f$ecm$npar),
           AIC = sapply(e, AIC), BIC = sapply(e, BIC))
lr <- 2 * (logLik(e$latreg) - logLik(e$null))                 # null is nested in latreg
df <- e$latreg$ecm$npar - e$null$ecm$npar
c(LR = lr, df = df, p = pchisq(lr, df, lower.tail = FALSE))
# Bayesian marginal criteria: the persons are integrated out, so different latent structures compare;
# the chains start at the ECM solutions
g <- lapply(e, gibbs, n_iter = 2000, n_chain = 2, init = "ecm", seed = 1)
do.call(rbind, lapply(g, function(f) fit_indices(f)[, c("DIC", "WAIC", "LOOIC", "SE_LOOIC")]))
