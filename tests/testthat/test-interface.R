# The interface of task 3.1: ecm(data, model = ...), gibbs(data, model = ...), and the
# deprecation warnings of task 3.2.

small_data <- function(model = "cross", n_feat = 0) {
  set.seed(1)
  cond <- set_cond(n_subj = 150, n_item = 5, n_feat = n_feat)
  sim_data(cond, sim_para(cond, model), model)
}

test_that("model names the built-in model and says what is not available", {
  dat <- small_data()
  expect_error(ecm(dat), "name the model")
  expect_error(ecm(dat, model = "crossed"), "unknown model \"crossed\"")
  expect_error(ecm(dat, model = "mixture"), "not available yet")
  expect_error(gibbs(dat, model = "mixture"), "not available yet")
  expect_error(ecm(dat, model = "cross", itemtype = "grm"), "not available yet")
  expect_error(ecm(dat, model = "cross", itemtype = c("2pl", "1pl", "2pl", "2pl", "2pl")), "per item are not available")
  expect_error(ecm(dat, model = "cross", itemtype = c("2pl", "1pl")), "one per item")
  expect_error(ecm(dat, model = "cross", itemtype = "3pl"), "unknown itemtype")
  expect_error(ecm(dat$Y, model = "cross"), "input_data")
  expect_error(gibbs(dat$Y, model = "cross"), "input_data")
})

test_that("model names: any case, did-you-mean, only available models listed, the same names in sim_para()/sim_data()", {
  dat <- small_data()
  fit <- ecm(dat, model = "Cross", se = FALSE)
  expect_equal(estimates(fit), estimates(ecm(dat, model = "cross", se = FALSE)))
  header <- capture.output(print(fit))[1]
  expect_equal(header, "model = \"cross\" (rtirt_cross)")
  expect_error(ecm(dat, model = "crosss"), "did you mean \"cross\"")
  expect_error(ecm(dat, model = "mlirt"), "class of model = \"irt\"; did you mean \"irt\"")
  expect_error(ecm(dat, model = "rtirt_cross"), "did you mean \"cross\"")
  expect_error(ecm(dat, model = "2pl"), "item type \\(itemtype = \"2pl\"\\)")
  msg <- tryCatch(ecm(dat, model = "foo"), error = conditionMessage)
  expect_match(msg, "\"latent\" +speed regressed on ability")       # one line per model
  expect_false(grepl("mixture", msg))                                # only the available ones
  expect_error(ecm(dat, model = "null", itemtype = "rasch"), "Rasch model is itemtype = \"1pl\"")
  expect_false(grepl("grm", tryCatch(ecm(dat, model = "null", itemtype = "3pl"), error = conditionMessage)))
  expect_error(gibbs(dat, "cross", 1000, n_chain = 2), "third argument of ecm\\(\\) and gibbs\\(\\) is itemtype.*n_iter = 1000")
  # sim_para() and sim_data() take "irt" (and the old "mlirt"), with the same draws
  cond <- set_cond(n_subj = 50, n_item = 4, n_feat = 1)
  set.seed(2); para_irt <- sim_para(cond, "irt"); set.seed(2); para_mlirt <- sim_para(cond, "mlirt")
  expect_identical(para_irt, para_mlirt)
  set.seed(3); d_irt <- sim_data(cond, para_irt, "irt"); set.seed(3); d_mlirt <- sim_data(cond, para_irt, "mlirt")
  expect_identical(d_irt, d_mlirt)
  set.seed(4); default <- sim_para(cond); set.seed(4); cross <- sim_para(cond, "cross")
  expect_identical(default, cross)
  expect_error(sim_para(cond, "crosss"), "did you mean \"cross\"")
  expect_error(sim_para(set_cond(n_subj = 10, n_item = 3), "latreg"), "set n_feat >= 1")
  gq <- gibbs(dat, model = "cross", quantile = 0.5, n_iter = 20, n_chain = 1, seed = 1,
              control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_match(capture.output(print(gq))[1], "^model = \"cross\", quantile = 0.5 \\(rtirt_cross\\)$")
})

test_that("itemtype may be given once per item when all are equal", {
  dat <- small_data()
  fit_one <- ecm(dat, model = "cross", itemtype = "1pl", se = FALSE)
  fit_each <- ecm(dat, model = "cross", itemtype = rep("1pl", 5), se = FALSE)
  expect_equal(estimates(fit_one), estimates(fit_each))
})

test_that("cov picks covariates; models without covariates drop them quietly", {
  dat <- small_data("latent", n_feat = 2)
  fit_all <- ecm(dat, model = "latent", se = FALSE)
  fit_x1 <- ecm(dat, model = "latent", cov = "x1", se = FALSE)
  expect_true(any(estimates(fit_all)$parameter == "beta[x2]"))
  expect_false(any(estimates(fit_x1)$parameter == "beta[x2]"))
  expect_error(ecm(dat, model = "latent", cov = "x9"), "not in the data")
  expect_no_warning(ecm(dat, model = "cross", se = FALSE))
  expect_warning(ecm(dat, model = "cross", cov = "x1", se = FALSE), "cov is ignored")
})

test_that("model = \"irt\" fits a 2PL without response times or covariates", {
  set.seed(2)
  Y <- matrix(rbinom(150 * 5, 1, 0.6), 150, 5)
  dat <- input_data(resp = Y)
  fit <- ecm(dat, model = "irt", se = FALSE)
  expect_equal(nrow(estimates(fit)), 10)
  expect_s3_class(scores(fit), "tbl_df")
})

test_that("gibbs() starts from an ecm() fit of the same data and model only", {
  dat <- small_data()
  fit_ecm <- ecm(dat, model = "cross", se = FALSE)
  fit <- gibbs(dat, model = "cross", init = fit_ecm, n_iter = 60, n_chain = 1, seed = 1,
               control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_false(is.null(fit$post))
  expect_error(gibbs(dat, model = "null", init = fit_ecm, n_iter = 60),
               "init is an ecm\\(\\) fit of model = \"cross\", but gibbs\\(\\) was asked for model = \"null\"")
  expect_error(gibbs(dat, model = "cross", speed_var = "fixed", init = fit_ecm, n_iter = 60),
               "init was fitted with speed_var = \"free\", gibbs\\(\\) was asked for speed_var = \"fixed\".*init = ecm\\(data, model = \"cross\", speed_var = \"fixed\"\\)")
  other <- small_data()
  other$Y[1, 1] <- 1 - other$Y[1, 1]
  expect_error(gibbs(other, model = "cross", init = fit_ecm, n_iter = 60), "init was fitted to other responses")
  fewer <- input_data(dat$Y[1:100, ], dat$T[1:100, ])
  expect_error(gibbs(fewer, model = "cross", init = fit_ecm, n_iter = 60), "init was fitted to 150 persons x 5 items, gibbs\\(\\) has 100 x 5")
  # init with quantile models (a fit of the mean model, a fit for several levels), init as text,
  # and a fit given as the data
  expect_error(gibbs(dat, model = "cross", quantile = 0.5, init = fit_ecm, n_iter = 60),
               "init was fitted with the mean model, gibbs\\(\\) was asked for quantile = 0.5.*init = ecm\\(data, model = \"cross\", quantile = 0.5\\)")
  expect_error(gibbs(dat, model = "cross", quantile = c(0.2, 0.8), init = fit_ecm, n_iter = 60),
               "starts one model; for the quantiles 0.2, 0.8 use init = \"ecm\"")
  expect_error(gibbs(dat, model = "cross", init = TRUE, n_iter = 60), "init must be \"random\", \"ecm\" or an ecm\\(\\) fit.*not TRUE")
  fit_upper <- gibbs(dat, model = "cross", init = "ECM", n_iter = 60, n_chain = 1, seed = 1,
                     control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_false(is.null(fit_upper$post))
  expect_error(gibbs(fit_ecm, model = "cross", n_iter = 60),
               "first argument of gibbs\\(\\) is the data from input_data\\(\\), not an ecm\\(\\) fit.*gibbs\\(data, model = \"cross\", init = fit\\)")
  expect_error(ecm(fit_ecm, model = "cross"), "first argument of ecm\\(\\) is the data")
})

test_that("tables are tibbles with the printed layout kept", {
  dat <- small_data()
  fit <- ecm(dat, model = "cross")
  expect_s3_class(estimates(fit), "tbl_df")
  expect_s3_class(fit_indices(fit), "tbl_df")
  tables <- summary(fit)
  expect_s3_class(tables$items, c("rtirt_table", "tbl_df"))
  expect_s3_class(tables$parameters, c("rtirt_partable", "rtirt_table", "tbl_df"))
  expect_output(print(tables), "Discrimination \\(a\\)")
  expect_output(print(tables$items), "item parameters")
})

test_that("deprecated functions warn, once per session unless asked for always", {
  withr::local_options(birtRcpp.deprecation = "always")
  dat <- small_data()
  expect_warning(rtirt_cross(dat), class = "deprecatedWarning")
  expect_warning(rtirt_cross(dat), "ecm\\(data, model = \"cross\"\\)")
  expect_warning(step_gibbs(a, b), class = "deprecatedWarning")
  expect_warning(step_mala(theta), class = "deprecatedWarning")
  expect_warning(step_shift(theta = 1, b = 1), class = "deprecatedWarning")
  Y <- dat$Y
  expect_warning(block_model(Y = Y, theta[person] ~ normal(0, 1), Y ~ bernoulli_logit(theta)),
                 class = "deprecatedWarning")
  fit_ecm <- suppressWarnings(ecm(rtirt_cross(dat), se = FALSE))
  expect_warning(gibbs(fit_ecm, init = "ecm", n_iter = 40, n_chain = 1, verbose = FALSE, fit_indices = FALSE),
                 "init = fit")
  expect_no_warning(GibbsRtIrtCross(set_cond(), dat))
})

test_that("without the option a deprecation warning is shown once", {
  withr::local_options(birtRcpp.deprecation = NULL)
  dat <- small_data()
  shown <- birtRcpp:::deprecation_shown
  rm(list = ls(shown), envir = shown)
  expect_warning(rtirt_null(dat), class = "deprecatedWarning")
  expect_no_warning(rtirt_null(dat))
})

test_that("equations() and algorithm() of data and a model name match the old constructors", {
  dat <- small_data()
  old_eq <- suppressWarnings(equations(rtirt_cross(dat, quantile = 0.25)))
  new_eq <- equations(dat, model = "cross", quantile = 0.25)
  expect_identical(as.character(new_eq), as.character(old_eq))
  old_alg <- suppressWarnings(algorithm(rtirt_null(dat), method = "ecm", format = "latex"))
  new_alg <- algorithm(dat, model = "null", method = "ecm", format = "latex")
  expect_identical(as.character(new_alg), as.character(old_alg))
  expect_error(equations(dat), "model")
})

test_that("ecm() fits quantile models only where gibbs() has them", {
  dat <- small_data()
  expect_error(ecm(dat, model = "null", quantile = 0.9), "quantile = is available for model = \"latent\" or \"cross\" only")
  expect_error(ecm(dat, model = "cross", quantile = 1.2), "quantile must be NULL \\(mean model\\) or a level in \\(0, 1\\)")
})

test_that("unknown arguments of ecm() and gibbs() are errors that show the right call", {
  dat <- small_data()
  # engine settings: the control = form, with the user's own data and model
  expect_error(ecm(dat, model = "cross", nodes = 61),
               "use ecm(dat, model = \"cross\", control = ecm_control(nodes = 61))", fixed = TRUE)
  expect_error(ecm(dat, model = "cross", max_iter = 5, verbose = TRUE), "are settings of ecm_control()", fixed = TRUE)
  expect_error(ecm(dat, model = "cross", nodes = 61, control = ecm_control(tol = 1e-6)),
               "control = ecm_control(tol = 1e-06, nodes = 61)", fixed = TRUE)
  expect_error(gibbs(dat, model = "cross", n_iter = 50, n_burnin = 10),
               "control = gibbs_control(n_burnin = 10)", fixed = TRUE)
  for (name in c("n_thin", "collapse", "fit_indices", "verbose")) {
    args <- c(list(dat, model = "cross"), stats::setNames(list(TRUE), name))
    expect_error(do.call(gibbs, args), "gibbs_control(", fixed = TRUE)
  }
  # misspelled names: did you mean
  expect_error(gibbs(dat, model = "cross", niter = 200), "unknown argument niter; did you mean n_iter?", fixed = TRUE)
  expect_error(gibbs(dat, model = "cross", chains = 1), "did you mean n_chain?", fixed = TRUE)
  expect_error(ecm(dat, model = "cross", covs = "x1"), "did you mean cov?", fixed = TRUE)
  expect_error(gibbs(dat, model = "cross", burnin = 50), "control = gibbs_control(n_burnin = 50)", fixed = TRUE)
  expect_error(ecm(dat, model = "cross", foo = 1), "unknown argument foo; the arguments of ecm() are", fixed = TRUE)
  # the other engine's arguments
  expect_error(ecm(dat, model = "cross", n_iter = 100), "n_iter is an argument of gibbs(), not of ecm()", fixed = TRUE)
  expect_error(gibbs(dat, model = "cross", nodes = 61), "setting of ecm_control() for ecm(), not of gibbs()", fixed = TRUE)
  # an error, not a warning after the run
  expect_error(gibbs(dat, model = "cross", niter = 200, n_chain = 1), "niter")
  # the deprecated interface keeps its arguments
  fit_old <- suppressWarnings(ecm(rtirt_cross(dat), nodes = 21, se = FALSE))
  fit_new <- ecm(dat, model = "cross", se = FALSE, control = ecm_control(nodes = 21))
  expect_equal(estimates(fit_old), estimates(fit_new))
})

test_that("control = takes the settings of its own engine; unknown settings name the closest one", {
  dat <- small_data()
  expect_s3_class(ecm_control(), "ecm_control")
  expect_s3_class(gibbs_control(), "gibbs_control")
  expect_error(gibbs(dat, model = "cross", n_iter = 20, n_chain = 1, control = ecm_control()),
               "control = ecm_control\\(\\) holds the settings of ecm\\(\\); gibbs\\(\\) takes control = gibbs_control")
  expect_error(ecm(dat, model = "cross", control = gibbs_control()),
               "holds the settings of gibbs\\(\\); ecm\\(\\) takes control = ecm_control")
  expect_error(gibbs_control(burnin = 50), "unknown setting burnin of gibbs_control\\(\\); did you mean n_burnin")
  expect_error(gibbs_control(n_iter = 200), "n_iter is an argument of gibbs\\(\\), not a setting of gibbs_control")
  expect_error(gibbs_control(nodes = 61), "nodes is a setting of ecm_control\\(\\)")
  expect_error(ecm(dat, model = "cross", control = list(nodez = 61)), "did you mean nodes")
  expect_error(ecm(dat, model = "cross", control = 61), "control must come from ecm_control")
  # a plain list of settings still works, as before
  fit_list <- ecm(dat, model = "cross", se = FALSE, control = list(nodes = 21))
  fit_ctrl <- ecm(dat, model = "cross", se = FALSE, control = ecm_control(nodes = 21))
  expect_equal(estimates(fit_list), estimates(fit_ctrl))
  expect_output(print(gibbs_control(n_burnin = 500)), "n_burnin = 500")
})

test_that("gibbs() refuses quantile = for models without a quantile version, as ecm() does", {
  dat <- small_data("latent", n_feat = 1)
  for (model in c("null", "latreg", "irt")) {
    expect_error(gibbs(dat, model = model, quantile = 0.9, n_iter = 20),
                 sprintf("quantile = is available for model = \"latent\" or \"cross\" only; model = \"%s\"", model), fixed = TRUE)
  }
  expect_error(gibbs(dat, model = "null", quantile = c(0.25, 0.75)), "\"latent\" or \"cross\" only")
  expect_error(equations(dat, model = "null", quantile = 0.5), "\"latent\" or \"cross\" only")
  fit <- gibbs(dat, model = "latent", quantile = 0.9, n_iter = 20, n_chain = 1, seed = 1,
               control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_true(fit$qr)
})

test_that("print() of a Gibbs fit warns with the criteria of convergence(); convergence() prints once and covers quantile fits", {
  dat <- small_data()
  quiet <- gibbs_control(verbose = FALSE, fit_indices = FALSE)
  one_chain <- gibbs(dat, model = "null", n_iter = 60, n_chain = 1, seed = 1, control = quiet)
  printed <- paste(capture.output(print(one_chain)), collapse = "\n")
  expect_match(printed, "1 chain, so R-hat cannot be computed (it needs n_chain >= 2)", fixed = TRUE)
  expect_match(printed, "parameters have ESS > 400 and R-hat < 1.05", fixed = TRUE)
  expect_output(print(convergence(one_chain)), "only the ESS is checked")

  two_chains <- gibbs(dat, model = "null", n_iter = 60, n_chain = 2, seed = 1, control = quiet)
  expect_silent(check <- convergence(two_chains))                   # assigned: not printed
  expect_s3_class(check, "rtirt_convergence")
  expect_lt(check$share, 1)
  printed <- paste(capture.output(print(two_chains)), collapse = "\n")
  expect_match(printed, sprintf("only %.1f%% of %d parameters have ESS > 400", 100 * check$share, check$n), fixed = TRUE)
  expect_false(grepl("1 chain", printed))
  shown <- capture.output(print(check))
  expect_equal(sum(grepl("parameters with ESS", shown)), 1)
  expect_false(any(grepl("^\\$share", shown)))                     # no second print of the raw list

  levels <- gibbs(dat, model = "cross", quantile = c(0.25, 0.75), n_iter = 60, n_chain = 2, seed = 1, control = quiet)
  check_q <- convergence(levels)
  expect_length(check_q$levels, 2)
  expect_equal(check_q$levels[[1]]$share, convergence(levels$models[[1]])$share)
  expect_output(print(check_q), "q = 0.75: ")
  expect_match(paste(capture.output(print(levels)), collapse = "\n"), "Warning (q = 0.25): only", fixed = TRUE)
})

test_that("sig marks only parameters where 0 is a meaningful value", {
  dat <- small_data("latreg", n_feat = 1)
  fit <- ecm(dat, model = "latreg")
  e <- estimates(fit)
  tested <- grepl("^(b|beta_ability|beta_speed)\\[", e$parameter) | e$parameter == "cor_ability_speed"
  expect_true(all(is.na(e$sig[!tested])))                          # a, lambda, sigma2t, var_speed
  expect_true(all(e$sig[tested] %in% c("", "*")))
  expect_true(all(is.na(e$sig[grepl("^a\\[", e$parameter)])))
  g <- gibbs(dat, model = "latreg", n_iter = 60, n_chain = 1, seed = 1,
             control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  eg <- estimates(g)
  expect_identical(is.na(eg$sig), !(grepl("^(b|beta_ability|beta_speed)\\[", eg$parameter) | eg$parameter == "cor_ability_speed"))
  cross <- estimates(ecm(small_data(), model = "cross", se = TRUE))
  expect_true(all(cross$sig[startsWith(cross$parameter, "rho[")] %in% c("", "*")))
  expect_identical(is.na(e$z), !tested)                              # no Wald z where sig is NA
  parameters <- summary(fit)$parameters
  expect_true(all(is.na(parameters$pvalue[!tested])))
  expect_false(anyNA(parameters$pvalue[tested]))
})

test_that("gibbs(prior = \"map\") says that prior = is a setting of ecm()", {
  dat <- small_data("null")
  expect_error(gibbs(dat, model = "null", prior = "map", n_iter = 20),   # partial match of priors =
               "prior = \"ml\" or \"map\" is a setting of ecm()", fixed = TRUE)
})

test_that("ECM fits with extreme estimates or standard errors warn in print() and summary()", {
  set.seed(9)
  cond <- set_cond(n_subj = 40, n_item = 10)
  small <- sim_data(cond, sim_para(cond, "cross"), "cross")
  fit <- ecm(small, model = "cross")                                 # b[item02] = -34, SE 4105
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "Warning: extreme estimates or standard errors (", fixed = TRUE)
  expect_match(printed, "b[item02] = -34.1", fixed = TRUE)
  expect_match(printed, "try ecm(..., prior = \"map\") or gibbs()", fixed = TRUE)
  expect_true(any(grepl("extreme estimates", summary(fit)$header)))
  expect_true(any(grepl("ECM, 41 Gauss-Hermite nodes", summary(fit)$header)))   # the method line of the summary
  map <- ecm(small, model = "cross", prior = "map")
  expect_false(grepl("extreme", paste(capture.output(print(map)), collapse = "\n")))
  expect_false(grepl("extreme", paste(capture.output(print(ecm(small_data(), model = "cross"))), collapse = "\n")))
})

test_that("messages and checks point to the current interface, and ecm() fits are accepted where models were", {
  dat <- small_data()
  fit <- ecm(dat, model = "cross", se = FALSE)
  expect_error(check_sampler(dat), "checks the sampler of a fit, not data: check_sampler(ecm(dat", fixed = TRUE)
  expect_error(check_sampler(list()), "model must be an ecm() or gibbs() fit, e.g. check_sampler(ecm(dat", fixed = TRUE)
  expect_error(check_ecm(steps = 1), "check_ecm(): give an ecm() fit, e.g. check_ecm(ecm(dat, model = \"cross\"))", fixed = TRUE)
  check <- check_ecm(fit, prior = "ml", n_start = 2, seed = 1)                   # an ecm() fit, as documented
  expect_true(attr(check, "passed"))
  expect_error(preknowledge_test(fit, compromised = 1:2), "needs a fit of model = \"null\", \"latreg\" or \"latent\"", fixed = TRUE)
  # score_test(): the impact message suggests the latent regression of ecm()
  skip_if_not_installed("strucchange")
  set.seed(3)
  z <- as.numeric(attr(dat, "true_para")$theta > 0)
  expect_message(score_test(fit, z), "ecm(input_data(..., cov = z), model = \"latreg\")", fixed = TRUE)
  with_cov <- input_data(resp = dat$Y, time = dat$T, cov = cbind(z = z - mean(z)))
  expect_no_message(score_test(ecm(with_cov, model = "latreg", se = FALSE), z))
})

test_that("precis() and summary() of quantile fits print once; table and parameter names are checked; plots", {
  dat <- small_data("latreg", n_feat = 1)
  fit <- ecm(dat, model = "latreg")
  g <- gibbs(dat, model = "latreg", init = fit, n_iter = 60, n_chain = 2, seed = 1,
             control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_silent(p <- precis(g))
  expect_s3_class(p, "rtirt_precis")
  expect_true(withVisible(precis(g))$visible)
  expect_output(print(p), "cor_ability_speed")
  qs <- gibbs(small_data("latent", n_feat = 1), model = "latent", quantile = c(0.3, 0.7), n_iter = 40, n_chain = 1,
              seed = 1, control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_silent(sq <- summary(qs))
  expect_output(print(sq), "quantiles 0.3, 0.7")
  # summary tables and estimates(pars = ) that select nothing
  s <- summary(fit)
  expect_error(print(s, tables = "item"), "unknown table \"item\" \\(did you mean \"items\"\\?\\)")
  expect_output(print(s, tables = "items"), "item parameters")
  expect_warning(none <- estimates(fit, pars = "rho"), "no parameters match \"rho\"; this fit has a, b, lambda")
  expect_equal(nrow(none), 0)
  expect_warning(none_g <- estimates(g, pars = "^rho"), "no parameters match")
  expect_equal(nrow(none_g), 0)
  expect_equal(names(none_g), names(estimates(g)))
  # plots: ECM fits have no traces; structural parameters; as.data.frame()
  expect_error(plot(fit, "trace"), "this is an ecm\\(\\) fit: it has no chains or draws; sample with gibbs\\(data, model = \"latreg\", init = fit\\)")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  drawn <- plot(g, "structural")
  expect_equal(drawn$parameter, c("beta_ability[x1]", "beta_speed[x1]", "cor_ability_speed", "var_speed"))
  expect_equal(nrow(plot(fit, "structural")), 4)
  expect_false(withVisible(plot(qs))$visible)
  expect_equal(as.data.frame(fit), as.data.frame(estimates(fit)))
  expect_match(attr(suppressMessages(summary(g))$fit, "title"), "compare models by WAIC or LOOIC")
})

test_that("anova() compares ecm() fits of the same data by likelihood ratio, AIC and BIC; gibbs() fits by WAIC and LOOIC", {
  dat <- small_data("latreg", n_feat = 1)
  f0 <- ecm(dat, model = "null", se = FALSE)
  f1 <- ecm(dat, model = "latreg", se = FALSE)
  table <- anova(f1, f0)                                   # ordered by the number of parameters
  expect_s3_class(table, "anova")
  expect_equal(rownames(table), c("f0", "f1"))
  expect_equal(table$Chisq[2], as.numeric(2 * (logLik(f1) - logLik(f0))))
  expect_equal(table$Df[2], f1$ecm$npar - f0$ecm$npar)
  expect_equal(table[["Pr(>Chisq)"]][2], stats::pchisq(table$Chisq[2], table$Df[2], lower.tail = FALSE))
  expect_equal(table$AIC, c(AIC(f0), AIC(f1)))
  expect_output(print(table), "f1: model = \"latreg\" \\(rtirt_latreg\\)")
  expect_error(anova(f1), "two or more fits")
  expect_error(anova(f1, ecm(dat, model = "irt", se = FALSE)), "one has response times and the other not")
  expect_error(anova(f1, ecm(dat, model = "latreg", prior = "map", se = FALSE)), "prior = \"map\"")
  expect_error(anova(f1, ecm(small_data("cross"), model = "cross", se = FALSE)), "fits of different data")
  f_cross <- ecm(dat, model = "cross", se = FALSE)
  expect_warning(anova(f0, f1, f_cross), "more parameters than f1 but a lower log-likelihood")
  g0 <- gibbs(dat, model = "null", init = f0, n_iter = 60, n_chain = 2, seed = 1, control = gibbs_control(verbose = FALSE))
  g1 <- gibbs(dat, model = "latreg", init = f1, n_iter = 60, n_chain = 2, seed = 1, control = gibbs_control(verbose = FALSE))
  criteria <- anova(g0, g1)
  expect_equal(criteria$WAIC, c(fit_indices(g0)$WAIC, fit_indices(g1)$WAIC))
  expect_error(anova(g0, f1), "not both")
})

test_that("gibbs() reports progress: a line at the start, one per chain with the time left; set_cond() says where its MCMC settings go", {
  dat <- small_data()
  messages <- character()
  fit <- withCallingHandlers(
    gibbs(dat, model = "cross", n_iter = 40, n_chain = 2, seed = 1),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_match(messages[1], "^model = \"cross\": sampling 2 chains x 40 iterations, one after the other; chain 1 of 2")
  expect_match(messages[2], "^model = \"cross\": chain 1 of 2 done \\([0-9.]+ s\\); (about .*|under 1 s) left")
  expect_match(messages[3], "^model = \"cross\": chain 2 of 2 done \\([0-9.]+ s\\)\\n$")
  expect_match(messages[4], "computing the marginal log-likelihood")
  quiet <- gibbs(dat, model = "cross", n_iter = 40, n_chain = 2, seed = 1, control = gibbs_control(verbose = FALSE))
  expect_identical(quiet$post$draws, fit$post$draws)                 # the messages do not touch the draws
  expect_silent(gibbs(dat, model = "cross", n_iter = 40, n_chain = 1, seed = 1,
                      control = gibbs_control(verbose = FALSE, fit_indices = FALSE)))
  expect_output(print(set_cond(n_subj = 10, n_item = 3)), "gibbs\\(data, model = \\.\\.\\.\\) takes n_iter and n_chain")
})

test_that("check_sampler(), check_ecm() and the Julia constructors given data say what they take", {
  dat <- small_data()
  expect_error(check_sampler(dat, model = "cross"), "check_sampler\\(\\) checks the sampler of a fit, not data: check_sampler\\(ecm\\(dat, model = \"cross\", se = FALSE\\)\\)")
  expect_error(check_sampler(dat), "convergence\\(\\) and plot\\(fit, \"trace\"\\)")
  expect_error(check_ecm(dat, model = "latreg"), "check_ecm\\(ecm\\(dat, model = \"latreg\"\\)\\)")
  expect_error(GibbsRtIrtCross(dat), "takes the conditions first.*GibbsRtIrtCross\\(set_cond\\(\\), dat\\). In R: gibbs\\(dat, model = \"cross\"\\)")
  expect_error(GibbsMlIrt(dat), "gibbs\\(dat, model = \"irt\"\\)")
  expect_s3_class(GibbsRtIrtCross(set_cond(), dat), "rtirt_cross")
})

test_that("algorithm(detail = FALSE) lists the steps without the formulas", {
  dat <- small_data()
  full <- capture.output(print(algorithm(dat, model = "cross")))
  brief <- capture.output(print(algorithm(dat, model = "cross", detail = FALSE)))
  expect_lt(length(brief), length(full) / 2)
  steps <- grep("^Step [0-9]+\\.", full, value = TRUE)
  expect_true(length(steps) > 3)
  expect_true(all(sub(":$", "", steps) %in% brief))
  expect_match(brief[length(brief)], "detail = TRUE")
  expect_equal(as.character(algorithm(dat, model = "cross", detail = FALSE)), as.character(algorithm(dat, model = "cross")))
})
