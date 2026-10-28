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
  expect_error(ecm(dat, model = "crossed"), "model must be one of")
  expect_error(ecm(dat, model = "mixture"), "not available yet")
  expect_error(gibbs(dat, model = "mixture"), "not available yet")
  expect_error(ecm(dat, model = "cross", itemtype = "grm"), "not available yet")
  expect_error(ecm(dat, model = "cross", itemtype = c("2pl", "1pl", "2pl", "2pl", "2pl")), "per item are not available")
  expect_error(ecm(dat, model = "cross", itemtype = c("2pl", "1pl")), "one per item")
  expect_error(ecm(dat, model = "cross", itemtype = "3pl"), "unknown itemtype")
  expect_error(ecm(dat$Y, model = "cross"), "input_data")
  expect_error(gibbs(dat$Y, model = "cross"), "input_data")
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
  expect_error(gibbs(dat, model = "null", init = fit_ecm, n_iter = 60), "same data, model")
  expect_error(gibbs(dat, model = "cross", speed_var = "fixed", init = fit_ecm, n_iter = 60), "same data, model")
  other <- small_data()
  other$Y[1, 1] <- 1 - other$Y[1, 1]
  expect_error(gibbs(other, model = "cross", init = fit_ecm, n_iter = 60), "same data, model")
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

test_that("ecm() refuses quantile models", {
  dat <- small_data()
  expect_error(ecm(dat, model = "cross", quantile = 0.9), "quantile models use gibbs")
})
