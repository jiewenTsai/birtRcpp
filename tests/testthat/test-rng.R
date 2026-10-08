test_that("Polya-Gamma draws have the right mean and variance", {
  set.seed(1); n <- 5e4
  for (z in c(0, 1, 4, 10)) {
    x <- birtRcpp:::.rpg1(rep(z, n))
    m <- if (z == 0) 0.25 else tanh(z / 2) / (2 * z)
    v <- if (z == 0) 1 / 24 else (sinh(z) - z) / (4 * z^3 * cosh(z / 2)^2)
    expect_lt(abs(mean(x) - m), 4 * sqrt(v / n))
    expect_equal(var(x), v, tolerance = 0.05)
  }
})

test_that("truncated normal stays in its bounds and is exact in the far tail", {
  set.seed(2)
  x <- birtRcpp:::.rtnorm(2e4, 0.5, 1, 8, 9)
  expect_true(all(x >= 8 & x <= 9))
  Q <- function(v) pnorm(v, 0.5, lower.tail = FALSE)
  expect_gt(suppressWarnings(ks.test((Q(8) - Q(x)) / (Q(8) - Q(9)), "punif")$p.value), 0.001)
  y <- birtRcpp:::.rtnorm(2e4, 0, 1, -4, 4)
  expect_equal(mean(y), 0, tolerance = 0.03)
})

test_that("inverse Gaussian moments", {
  set.seed(3); x <- birtRcpp:::.rinvgauss(1e5, 2, 3)
  expect_equal(mean(x), 2, tolerance = 0.02); expect_equal(var(x), 8 / 3, tolerance = 0.05)
})

test_that("both item samplers target the exact (a, b) posterior given theta", {
  set.seed(11)
  theta <- rnorm(150); y <- rbinom(150, 1, plogis(1.8 * (theta + 1.2)))     # an easy, discriminating item
  # exact posterior on a grid: a ~ N+(1, 1), b ~ N(0, 1) on [-4, 4]
  ga <- seq(0.01, 6, length.out = 400); gb <- seq(-4, 4, length.out = 400)
  lp <- outer(ga, gb, Vectorize(function(a, b) { eta <- a * (theta - b)
    sum(y * eta - log1p(exp(eta))) - (a - 1)^2 / 2 - b^2 / 2 }))
  w <- exp(lp - max(lp)); w <- w / sum(w)
  exact <- c(a = sum(w * ga), b = sum(t(t(w) * gb)))
  for (joint in c(FALSE, TRUE)) {
    d <- birtRcpp:::.item_chain(y, theta, 20000, joint)[-(1:1000), ]
    se <- apply(d, 2, sd) / sqrt(coda::effectiveSize(coda::mcmc(d)))
    expect_true(all(abs(colMeans(d) - exact) < 4.5 * se + 1e-3), info = paste("joint =", joint))
  }
})
