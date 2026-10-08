test_that("calc_normal_part() matches the hand computation and skips missing values", {
  set.seed(1); N <- 30; K <- 4
  x <- matrix(rnorm(N * K), N); z <- matrix(rnorm(N * K), N); w <- matrix(runif(N * K), N); rest <- matrix(rnorm(N * K), N)
  g <- calc_normal_part(x, z, w, rest)
  expect_equal(g$prec, colSums(w * x^2)); expect_equal(g$num, colSums(w * x * (z - rest)))
  gp <- calc_normal_part(x, z, w, rest, by = "person"); expect_equal(gp$num, rowSums(w * x * (z - rest)))
  z2 <- z; z2[1, 1] <- NA; g2 <- calc_normal_part(x, z2, w, rest)
  expect_equal(g2$prec[1], sum((w * x^2)[-1, 1])); expect_equal(g2$prec[-1], g$prec[-1])
  expect_equal(calc_normal_part(2, z, 1)$prec, rep(4 * N, K))
  expect_error(calc_normal_part(x[, 1:2], z), "matrices like z")
})

test_that("draw_mvn() draws from the precision form, batched and truncated", {
  set.seed(2); Q <- matrix(c(2, .5, .5, 1), 2); h <- c(1, -1)
  x <- replicate(20000, draw_mvn(Q, h))
  expect_equal(rowMeans(x), drop(solve(Q, h)), tolerance = 0.02)
  expect_equal(cov(t(x)), solve(Q), tolerance = 0.03)
  y <- draw_mvn(array(rep(diag(2), 50), c(2, 2, 50)), matrix(0, 2, 50), lower = 0)
  expect_equal(dim(y), c(2L, 50L)); expect_true(all(y >= 0))
})

test_that("draw_items_mh() leaves the posterior of (a, b) invariant (logit and probit)", {
  skip_on_cran()
  set.seed(3); N <- 150; th <- rnorm(N)
  for (lk in c("logit", "probit")) {
    fm <- if (lk == "logit") "difficulty" else "intercept"
    y <- rbinom(N, 1, if (lk == "logit") plogis(1.5 * (th - 0.5)) else pnorm(1.5 * th - 0.8))
    a <- 1; b <- 0; out <- matrix(NA, 8000, 2)
    for (t in 1:8000) { r <- draw_items_mh(matrix(y), th, a, b, c(1, 0), diag(c(1, 1 / 9)), lk, fm); a <- r$a; b <- r$b; out[t, ] <- c(a, b) }
    ga <- seq(0.05, 5, length.out = 200); gb <- seq(-4, 4, length.out = 200)
    lp <- outer(ga, gb, Vectorize(function(A, B) { eta <- if (fm == "difficulty") A * (th - B) else A * th - B
      sum(if (lk == "logit") y * eta - log1p(exp(eta)) else y * pnorm(eta, log.p = TRUE) + (1 - y) * pnorm(-eta, log.p = TRUE)) +
        dnorm(A, 1, 1, log = TRUE) + dnorm(B, 0, 3, log = TRUE) }))
    w <- exp(lp - max(lp)); w <- w / sum(w)
    expect_equal(colMeans(out[-(1:500), ]), c(sum(w * ga), sum(t(w) * gb)), tolerance = 0.03)
  }
})

test_that("draw_invwishart(), draw_mean_mvn(), draw_cor() and calc_person_fit() work", {
  set.seed(4)
  S <- replicate(4000, draw_invwishart(10, diag(2) * 7)); expect_equal(apply(S, 1:2, mean), diag(2), tolerance = 0.1)  # mean scale / (df - 3)
  X <- matrix(rnorm(200, 3), 100); m <- replicate(2000, draw_mean_mvn(X, diag(2), c(0, 0), diag(2) * 100))
  expect_equal(rowMeans(m), colMeans(X), tolerance = 0.02)
  x1 <- rnorm(400); x2 <- 0.5 * x1 + sqrt(0.75) * rnorm(400); r <- 0; rr <- numeric(3000)
  for (t in 1:3000) rr[t] <- r <- draw_cor(x1, x2, r)
  expect_lt(abs(mean(rr[-(1:300)]) - cor(x1, x2)), 0.05)
  L <- matrix(c(1, 2, NA, 4), 2); pf <- calc_person_fit(L, matrix(0, 2, 2), c(1, 4))
  expect_equal(pf$lt, c(1, 2^2 / 1 + 4^2 / 4)); expect_equal(pf$df, c(1, 2))
})

test_that("draw_mala() leaves the target invariant and draw_cor(method = 'mala') mixes faster", {
  set.seed(5); n <- 400; mu <- rnorm(n); s <- runif(n, 0.5, 2)
  lp <- function(x) list(lp = -(x - mu)^2 / (2 * s^2), grad = -(x - mu) / s^2)
  x <- mu; eps <- rep(1, n); m1 <- m2 <- 0
  for (t in 1:3000) { x <- draw_mala(x, lp, eps); if (t <= 500) eps <- mcmc_adapt(eps, attr(x, "accepted"), t) else { m1 <- m1 + x; m2 <- m2 + x^2 } }
  m1 <- m1 / 2500; v <- m2 / 2500 - m1^2
  expect_lt(abs(mean((m1 - mu) / s)), 0.05); expect_equal(mean(v / s^2), 1, tolerance = 0.05)
  # matrix input: n independent bivariate normals
  X <- matrix(0, 50, 2); lp2 <- function(X) list(lp = -rowSums(X^2) / 2, grad = -X)
  X <- draw_mala(X, lp2, 0.5); expect_equal(dim(X), c(50L, 2L)); expect_length(attr(X, "accepted"), 50)
  # correlation: same posterior as the random walk, larger ESS
  x1 <- rnorm(500); x2 <- 0.6 * x1 + 0.8 * rnorm(500)
  run <- function(m) { r <- 0; st <- 0.05; out <- numeric(4000)
    for (t in 1:5000) { r <- draw_cor(x1, x2, r, steps = 1, step = st, method = m)
      if (t <= 1000) st <- mcmc_adapt(st, attr(r, "accepted"), t, if (m == "rw") 0.44 else 0.574) else out[t - 1000] <- r }
    out }
  rw <- run("rw"); ml <- run("mala")
  expect_lt(abs(mean(rw) - mean(ml)), 0.005); expect_equal(sd(ml), sd(rw), tolerance = 0.15)
  ess <- function(z) { a <- acf(z, plot = FALSE, lag.max = 50)$acf[-1]; length(z) / (1 + 2 * sum(a[seq_len(match(TRUE, a < 0.05, 50))])) }
  expect_gt(ess(ml), 1.5 * ess(rw))
})
