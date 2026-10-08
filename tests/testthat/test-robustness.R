# regressions from the stress test (scripts/41_stress_birtRcpp.R)
mk <- function(N = 80, K = 5, seed = 1) {
  set.seed(seed); th <- stats::rnorm(N)
  Y <- matrix(stats::rbinom(N * K, 1, stats::plogis(outer(th, rep(1, K)))), N)
  T <- exp(matrix(3 + stats::rnorm(N * K, 0, 0.5), N)); colnames(Y) <- colnames(T) <- sprintf("i%d", 1:K)
  list(Y = Y, T = T)
}

test_that("a time column without variation warns and does not break the samplers", {
  r <- mk(); r$T[, 3] <- 30
  expect_warning(d <- input_data(r$Y, r$T), "without variation")
  f <- gibbs(rtirt_cross(d), n_iter = 200, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
  expect_true(all(is.finite(coef(f))))
  f2 <- gibbs(rtirt_cross(d, quantile = 0.3), n_iter = 200, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
  expect_true(all(is.finite(coef(f2))))
})

test_that("too few kept draws is an error, not a coda failure", {
  r <- mk(); d <- input_data(r$Y, r$T)
  expect_error(gibbs(rtirt_cross(d), n_iter = 2, n_chain = 2, verbose = FALSE), "at least 2")
  expect_error(gibbs(rtirt_cross(d), n_iter = 100, n_thin = 60, verbose = FALSE), "at least 2")
  f <- gibbs(rtirt_cross(d), n_iter = 6, n_chain = 2, verbose = FALSE, fit_indices = FALSE)
  expect_s3_class(summary(f), "rtirt_summary")
})

test_that("same seed, same draws; item names and ids are kept", {
  r <- mk(); colnames(r$Y) <- colnames(r$T) <- c("item 1", "item.2", "x-3", "4", "q_5")
  d <- input_data(r$Y, r$T, id = sprintf("P%02d", 1:80))
  a <- gibbs(rtirt_cross(d), n_iter = 60, n_chain = 2, seed = 7, verbose = FALSE, fit_indices = FALSE)
  b <- gibbs(rtirt_cross(d), n_iter = 60, n_chain = 2, seed = 7, verbose = FALSE, fit_indices = FALSE)
  expect_identical(as.matrix(coda::as.mcmc.list(a)), as.matrix(coda::as.mcmc.list(b)))
  expect_identical(summary(a)$items$item, colnames(r$Y))
  expect_identical(scores(a)$id, sprintf("P%02d", 1:80))
  l <- as.character(equations(a, "latex"))
  expect_true(any(grepl("q\\_5", l, fixed = TRUE)))
})

test_that("LaTeX special characters in item names are escaped", {
  r <- mk(); colnames(r$Y) <- colnames(r$T) <- c("a&b", "50%", "x#1", "y$", "q_5")
  l <- as.character(equations(rtirt_cross(input_data(r$Y, r$T)), "latex"))
  set <- l[grepl("mathrm{a", l, fixed = TRUE)][1]
  for (p in c("\\mathrm{a\\&b}", "\\mathrm{50\\%}", "\\mathrm{q\\_5}")) expect_true(grepl(p, set, fixed = TRUE), info = p)
})
