sim2pl <- function(N = 300, K = 8, seed = 1) {
  set.seed(seed); a <- runif(K, 0.8, 2); b <- rnorm(K)
  Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(rnorm(N), a), 2, a * b))), N)
  Y[sample(length(Y), 50)] <- NA
  list(Y = Y, a = a, b = b, N = N, K = K)
}

test_that("common and person-specific nodes give the same E-step quantities", {
  d <- sim2pl(); nd <- calc_nodes(21); TH <- matrix(nd$z, d$N, 21, byrow = TRUE)
  L1 <- calc_loglik_nodes(d$Y, d$a, d$b, nd$z); L2 <- calc_loglik_nodes(d$Y, d$a, d$b, TH)
  expect_equal(L1, L2, tolerance = 1e-10)
  # direct check of one cell
  y <- d$Y[3, ]; eta <- d$a * (nd$z[5] - d$b); ok <- !is.na(y)
  expect_equal(L1[3, 5], sum(y[ok] * eta[ok] - log1p(exp(eta[ok]))), tolerance = 1e-10)
  W <- calc_node_weights(L1, nd$logw)$W
  expect_equal(rowSums(W), rep(1, d$N), tolerance = 1e-12)
  for (sc in c(FALSE, TRUE)) {
    S1 <- calc_item_sums(d$Y, d$a, d$b, nd$z, W, score = sc); S2 <- calc_item_sums(d$Y, d$a, d$b, TH, W, score = sc)
    expect_equal(S1$S0, S2$S0, tolerance = 1e-10); expect_equal(S1$S1, S2$S1, tolerance = 1e-10)
  }
})

test_that("the component EM reaches the maximum (score zero) and SQUAREM agrees with plain EM", {
  d <- sim2pl(N = 500, K = 6); Y <- d$Y; K <- d$K; nd <- calc_nodes(31); kap <- Y - 0.5
  kap0 <- kap; kap0[is.na(kap0)] <- 0
  Fmap <- function(x) {
    a <- exp(x[1:K]); b <- x[K + 1:K]
    nw <- calc_node_weights(calc_loglik_nodes(Y, a, b, nd$z), nd$logw)
    S <- calc_item_sums(Y, a, b, nd$z, nw$W)
    new <- update_items_logistic(S, colSums(kap0), colSums(kap0 * drop(nw$W %*% nd$z)), a, b)
    list(x = c(log(new$a), new$b), obj = sum(nw$ll))
  }
  f1 <- em_squarem(rep(0, 2 * K), Fmap, tol = 1e-12)
  f0 <- em_squarem(rep(0, 2 * K), Fmap, tol = 1e-12, accelerate = FALSE, max_iter = 5000)
  expect_true(f1$converged); expect_lt(f1$iterations, f0$iterations)
  expect_equal(f1$x, f0$x, tolerance = 1e-4)
  expect_true(all(diff(f0$trace) > -1e-8))                     # plain EM is monotone
  score <- function(x) {                                       # gradient in (log a, b) by the Fisher identity
    a <- exp(x[1:K]); b <- x[K + 1:K]
    nw <- calc_node_weights(calc_loglik_nodes(Y, a, b, nd$z), nd$logw)
    S <- calc_item_sums(Y, a, b, nd$z, nw$W, score = TRUE)
    c(a * (S$S1 - b * S$S0), -a * S$S0)
  }
  expect_lt(max(abs(score(f1$x))), 1e-3)
  H <- calc_hessian(score, f1$x)
  expect_equal(H, t(H)); expect_true(all(eigen(H, symmetric = TRUE)$values < 0))
})

test_that("calc_normal_given_nodes integrates a normal latent variable exactly", {
  set.seed(2); K <- 4; o <- matrix(rnorm(3 * K), 3); w <- c(1, 2, 0.5, 3); m <- c(-1, 0, 1); v <- 0.7
  r <- calc_normal_given_nodes(o, w, m, v)
  f <- function(u, j, q) exp(sum(dnorm(o[j, ], u, 1 / sqrt(w), log = TRUE)) + dnorm(u, m[q], sqrt(v), log = TRUE))
  I <- integrate(Vectorize(function(u) f(u, 2, 3)), -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(r$logI[2, 3], log(I), tolerance = 1e-7)
  M <- integrate(Vectorize(function(u) u * f(u, 2, 3)), -Inf, Inf, rel.tol = 1e-10)$value / I
  expect_equal(r$mean[2, 3], M, tolerance = 1e-7)
})

test_that("pg_mean() is E[omega] of PG(1, eta), elementwise, and keeps the dimensions", {
  eta <- matrix(c(-40, -2, -1e-5, 0, 1e-5, 0.3, 5, 700), 2)
  ref <- ifelse(abs(eta) < 1e-4, 0.25 - eta^2 / 48, tanh(eta / 2) / (2 * eta))
  expect_equal(pg_mean(eta), ref, tolerance = 1e-15)
  expect_equal(dim(pg_mean(eta)), c(2L, 4L))
  set.seed(4); z <- 1.3; expect_equal(mean(draw_pg(rep(z, 4e4))), pg_mean(z), tolerance = 0.01)   # against Polya-Gamma draws
})

test_that("the 2PL log-likelihood stays finite for very large eta (softplus), in every place it is computed", {
  set.seed(6); N <- 30; K <- 4; th <- rnorm(N); Y <- matrix(rbinom(N * K, 1, 0.5), N)
  a <- c(400, 1, 1, 1); b <- c(-6, 0, 0.5, -0.5)                  # eta of item 1 about 2400: exp(eta) overflows
  sp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  eta <- outer(th, a) - matrix(a * b, N, K, byrow = TRUE)
  ref <- rowSums(Y * eta - sp(eta))
  expect_equal(drop(calc_loglik_nodes(Y, a, b, matrix(th))), ref, tolerance = 1e-12)
  expect_false(is.finite(sum(Y * eta - log1p(exp(eta)))))         # the old form
  cond <- set_cond(n_subj = N, n_item = K, n_feat = 1)
  d <- sim_data(cond, sim_para(cond, "mlirt"), "mlirt"); d$Y[] <- Y
  f <- gibbs(mlirt(d), n_iter = 20, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
  f$post$mean$a <- a; f$post$mean$b <- b; f$post$mean$theta <- th
  expect_true(is.finite(birtRcpp:::complete_loglik(f)))           # dic()
  p <- list(a = a, b = b, beta = f$post$mean$beta)
  expect_true(all(is.finite(birtRcpp:::marg_ll_draw(f, p))))       # marginal_loglik()
})

test_that("one formula for the normal integral: calc_normal_given_nodes() = the sums form used by ecm() and marginal_loglik()", {
  set.seed(3); N <- 5; K <- 4; Q <- 3; o <- matrix(rnorm(N * K), N); w <- c(1, 2, 0.5, 3); m <- matrix(rnorm(N * Q), N); v <- 0.4
  r <- calc_normal_given_nodes(o, w, m, v)
  # the Sherman-Morrison form of the old marginal_loglik(): o ~ N(m 1, diag(1/w) + v 1 1')
  P <- 1 / v + sum(w)
  sm <- sapply(seq_len(Q), function(q) sapply(seq_len(N), function(j) {
    r0 <- o[j, ] - m[j, q]; s1 <- sum(w * r0); s2 <- sum(w * r0^2)
    -0.5 * K * log(2 * pi) + 0.5 * sum(log(w)) - 0.5 * log(1 + v * sum(w)) - 0.5 * (s2 - v * s1^2 / (1 + v * sum(w))) }))
  expect_equal(r$logI, sm, tolerance = 1e-12)
  expect_equal(r$var, 1 / P)
})
