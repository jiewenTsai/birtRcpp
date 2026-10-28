# Example 6. Extracting the summary tables and writing the model for a paper
library(birtRcpp)
set.seed(6)
cond <- set_cond(n_subj = 400, n_item = 8)
dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
fit <- ecm(dat, model = "cross")
s <- summary(fit)
names(s)
utils::write.csv(as.data.frame(s), file.path(tempdir(), "parameters_ecm.csv"), row.names = FALSE) # long table
items <- dplyr::left_join(s$items, s$se, by = "item", suffix = c("", "_se")) # estimates next to SEs
items
eq <- equations(fit, "latex")
writeLines(as.character(eq), file.path(tempdir(), "model.tex")) # \[ \input{model.tex} \] in the paper
# the Bayesian version of the same tables: posterior SDs, R-hat, ESS and the priors
g <- gibbs(dat, model = "cross", init = fit, n_iter = 1500, n_chain = 2, seed = 1)
sg <- summary(g)
utils::write.csv(as.data.frame(sg), file.path(tempdir(), "parameters_gibbs.csv"), row.names = FALSE)
head(as.data.frame(sg))
equations(g) # with the priors
equations(dat, model = "cross", quantile = 0.25) # also before estimation
