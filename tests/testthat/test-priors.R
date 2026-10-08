test_that("rtirt_priors() validates its input and prints", {
  p <- rtirt_priors()
  expect_s3_class(p, "rtirt_priors"); expect_length(birtRcpp:::prior_vector(p), 18)
  expect_equal(p$a, c(mean = 1, sd = 1)); expect_equal(p$sigma2t_family, "half-t")
  expect_equal(rtirt_priors(a = c(sd = 0.5, mean = 1))$a, c(mean = 1, sd = 0.5))
  expect_equal(rtirt_priors(sigma2t = c(scale = 1, shape = 2))$sigma2t_family, "inverse-gamma")
  expect_error(rtirt_priors(a = c(mean = 1, sd = -1)), "positive")
  expect_error(rtirt_priors(a = c(m = 1, s = 1)), "needs the elements")
  j <- rtirt_priors_julia(); expect_null(j$lambda); expect_equal(unname(j$b["upper"]), 4)
  expect_error(rtirt_priors(b = c(0, 1, 2, -2)), "lower < mean < upper")
  expect_error(rtirt_priors(sigma2t = c(0, 1)), "positive")
  expect_error(rtirt_priors(rho = c(1, 2)), "rho must be c\\(sd\\)")
  expect_output(print(p), "half-t\\(3, 1\\)"); expect_output(print(j), "N\\+\\(1, 1\\)")
})

test_that("default priors reproduce the sampler and tighter priors shrink", {
  set.seed(3)
  cond <- set_cond(n_subj = 200, n_item = 5)
  m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
  run <- function(...) as.matrix(gibbs(m, n_iter = 300, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE, ...)$post$draws)
  expect_identical(run(), run(priors = rtirt_priors()))
  g <- gibbs(m, n_iter = 800, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE,
             priors = rtirt_priors(a = c(mean = 1, sd = 0.05), rho = c(sd = 0.02)))
  cf <- coef(g)
  expect_lt(max(abs(cf[grep("^a\\[", names(cf))] - 1)), 0.15)
  expect_lt(max(abs(cf[grep("^rho\\[", names(cf))])), 0.06)
  expect_true(any(grepl("0.05^2", summary(g)$parameters$prior, fixed = TRUE)))
  e <- ecm(m, prior = "map", se = FALSE, priors = rtirt_priors(a = c(mean = 1, sd = 0.05)))
  expect_lt(max(abs(e$ecm$par$a - 1)), 0.15)
  expect_error(gibbs(m, n_iter = 100, priors = list(a = 1), verbose = FALSE), "rtirt_priors")
})

test_that("the dissertation priors reproduce the earlier sampler and the defaults use the Gibbs steps", {
  set.seed(4)
  cond <- set_cond(n_subj = 200, n_item = 5)
  m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
  g <- gibbs(m, n_iter = 600, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
  cf <- coef(g); expect_true(all(cf[grep("^a\\[", names(cf))] > 0)); expect_true(all(is.finite(cf)))
  gj <- gibbs(m, n_iter = 600, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE, priors = rtirt_priors_julia())
  expect_true(any(grepl("N+(1, 1)", summary(gj)$parameters$prior, fixed = TRUE)))
  expect_true(any(grepl("half-t", summary(g)$parameters$prior, fixed = TRUE)))
  e <- ecm(m, prior = "map", se = FALSE); expect_true(all(e$ecm$par$a > 0 & e$ecm$par$sigma2t > 0))
})

test_that("the half-t steps (Huang & Wand auxiliary variable) leave the posterior invariant", {
  set.seed(5)
  chain <- function(alpha, beta, n = 40000) { v <- numeric(n); x <- 1
    for (i in seq_len(n)) v[i] <- x <- draw_var(alpha, beta, x, c(df = 3, scale = 1)); sqrt(v) }
  s0 <- chain(0, 0)                                   # no data: the prior, |t_3|
  expect_lt(abs(stats::median(s0) - stats::qt(0.75, 3)), 0.03)
  expect_lt(abs(mean(s0 < 2) - (2 * stats::pt(2, 3) - 1)), 0.01)
  alpha <- 10; beta <- 3                              # 20 residuals with sum of squares 6
  lp <- function(s) -2 * alpha * log(s) - beta / s^2 - 2 * log1p(s^2 / 3)
  den <- stats::integrate(function(s) exp(lp(s)), 0, Inf)$value
  m <- stats::integrate(function(s) s * exp(lp(s)), 0, Inf)$value / den
  expect_lt(abs(mean(chain(alpha, beta)) - m), 0.005)
})

test_that("only rtirt_priors_julia() centres lambda on the data", {
  expect_error(rtirt_priors(lambda = NULL), "only used by rtirt_priors_julia")
  expect_equal(unname(rtirt_priors()$lambda), c(0, 10, -Inf))
  L <- matrix(rnorm(40, 4), 10)
  expect_equal(birtRcpp:::lambda_prior(rtirt_priors(), L), c(0, 10, -Inf))
  expect_equal(birtRcpp:::lambda_prior(rtirt_priors_julia(), L)[1:2], c(mean(L), sd(as.vector(L))))
  expect_true(is.na(birtRcpp:::prior_vector(rtirt_priors_julia())[7]))
  expect_false(is.na(birtRcpp:::prior_vector(rtirt_priors())[7]))
})
