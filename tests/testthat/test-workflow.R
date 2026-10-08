sim2pl_w <- function(N = 300, K = 6, seed = 1) {
  set.seed(seed); a <- runif(K, 0.8, 2); b <- rnorm(K)
  list(Y = matrix(rbinom(N * K, 1, plogis(sweep(outer(rnorm(N), a), 2, a * b))), N), a = a, b = b, N = N, K = K)
}

test_that("normal parts add and draw like the precision form", {
  p <- normal_prior(c(1, 2), 2) + normal_part(c(1, 1), c(0.5, 0))
  expect_equal(p$prec, c(1.25, 1.25)); expect_equal(p$num, c(0.75, 0.5))
  set.seed(1); x <- draw(p); set.seed(1); expect_equal(x, draw_gauss(p$prec, p$num))
  expect_true(all(draw(normal_prior(rep(0, 50), 1), lower = 0) >= 0))
  expect_error(draw(normal_part(0, 1)), "precision must be positive")
  expect_error(normal_prior(0, 1) + list(prec = 1, num = 1), "only normal parts")
  expect_output(print(p), "normal part of 2")
})

test_that("run_sampler recovers a 2PL and keeps person means; errors are readable", {
  d <- sim2pl_w()
  fit <- run_sampler(list(Y = d$Y), init = function(x) list(theta = rnorm(d$N), a = rep(1, d$K), b = rep(0, d$K)),
    step = function(s, x) { s$omega <- draw_omega(s$a, s$b, s$theta); s[c("a", "b")] <- draw_items_pg(x$Y, s$omega, s$theta, s$a)
      s$theta <- draw(irt_part(x$Y, s$omega, s$a, s$b) + normal_prior(0, 1)); s },
    monitor = c("a", "b"), keep_mean = "theta", n_iter = 1500, seed = 3)
  sm <- summary(fit)
  expect_equal(sm$parameter, c(sprintf("a[%d]", 1:6), sprintf("b[%d]", 1:6)))
  expect_lt(sqrt(mean((sm$mean[7:12] - d$b)^2)), 0.3)
  expect_equal(nrow(fit$mean$theta), d$N); expect_equal(coda::nchain(fit$draws), 2)
  expect_error(run_sampler(list(), init = function(x) list(a = 1), step = function(s, x) 1, monitor = "a", n_iter = 10), "must return the state")
  expect_error(run_sampler(list(), init = function(x) list(a = 1), step = function(s, x) s, monitor = "zz", n_iter = 10), "not in the starting state")
})

test_that("input checks and missing responses in the blocks", {
  d <- sim2pl_w(N = 100, K = 4); th <- rnorm(100)
  expect_error(draw_omega(th, d$a, d$b), "items first, persons last")
  expect_error(irt_part(d$Y, matrix(1, 100, 4), d$a[-1], d$b), "a has 3 values")
  expect_error(draw_items_pg(d$Y, matrix(1, 4, 100), th, d$a), "omega is 4 x 100")
  Y <- d$Y; Y[1, 1] <- NA
  om <- draw_omega(d$a, d$b, th)
  p1 <- irt_part(Y, om, d$a, d$b); p0 <- irt_part(d$Y[, -1], om[, -1], d$a[-1], d$b[-1])
  expect_equal(p1$prec[1], p0$prec[1]); expect_equal(p1$num[1], p0$num[1]); expect_false(anyNA(p1$num))
  expect_false(anyNA(unlist(draw_items_pg(Y, om, th, d$a)))); expect_false(anyNA(unlist(draw_items_joint(Y, th, d$a, d$b))))
})

test_that("draw_mala: explicit joint target and checked logpost", {
  lp <- function(x) list(lp = -sum(x^2) / 2, grad = -x)
  expect_error(draw_mala(c(0, 0), lp, 1), "joint = TRUE")
  expect_error(draw_mala(c(0, 0), function(x) -sum(x^2), 1), "must return list")
  x <- draw_mala(matrix(0, 2, 2), lp, 0.5, joint = TRUE); expect_equal(dim(x), c(2L, 2L)); expect_length(attr(x, "accepted"), 1)
})

test_that("draw_shift draws the shift from the priors of the shifted parameters", {
  set.seed(2); cs <- replicate(4000, attr(draw_shift(list(b = c(1, 2), g0 = 0, theta = 1:3), prior = list(b = c(0, 1), g0 = c(0, 1))), "shift"))
  # c ~ N(h / Q, 1 / Q): Q = 3, h = -(1 + 2) - 0 = -3
  expect_equal(mean(cs), -1, tolerance = 0.03); expect_equal(var(cs), 1 / 3, tolerance = 0.05)
  s <- draw_shift(list(b = c(1, 2), theta = 1:3), prior = list(b = c(0, 1)))
  expect_equal(s$theta - 1:3, rep(attr(s, "shift"), 3))
  expect_error(draw_shift(list(b = 1), prior = list(zz = c(0, 1))), "not in values")
})

test_that("run_em: the maximum, standard errors from the Hessian, MAP and missing data", {
  d <- sim2pl_w(N = 600, K = 5, seed = 4)
  st <- function(s, x) { e <- estep_irt(x$Y, s$a, s$b); s[c("a", "b")] <- mstep_items(e, s$a, s$b); s$.loglik <- e$loglik; s }
  ll <- function(s, x) estep_irt(x$Y, s$a, s$b)$loglik
  fit <- run_em(list(Y = d$Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), st, loglik = ll, positive = "a")
  expect_true(fit$converged); expect_false(anyNA(fit$estimates$se))
  f <- function(z) estep_irt(d$Y, exp(z[1:5]), z[6:10])$loglik; z <- c(log(fit$state$a), fit$state$b)
  g <- sapply(1:10, function(k) { e <- replace(numeric(10), k, 1e-5); (f(z + e) - f(z - e)) / 2e-5 })
  expect_lt(max(abs(g)), 1e-2)
  lpost <- function(s, x) estep_irt(x$Y, s$a, s$b)$loglik + sum(dnorm(s$a, 1, 0.2, log = TRUE)) + sum(dnorm(s$b, 0, 0.5, log = TRUE))
  stm <- function(s, x) { e <- estep_irt(x$Y, s$a, s$b)
    s$.loglik <- e$loglik + sum(dnorm(s$a, 1, 0.2, log = TRUE)) + sum(dnorm(s$b, 0, 0.5, log = TRUE))     # at the input state
    s[c("a", "b")] <- mstep_items(e, s$a, s$b, priors = rtirt_priors(a = c(1, 0.2), b = c(0, 0.5))); s }
  fm <- expect_no_warning(run_em(list(Y = d$Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), stm, loglik = lpost, positive = "a"))
  expect_lt(sd(fm$state$a), sd(fit$state$a))                       # shrunk towards 1
  # per-item prior means of b: the mode moves towards them
  bm <- c(-2, -1, 0, 1, 2)
  stb <- function(s, x) { e <- estep_irt(x$Y, s$a, s$b); s$.loglik <- e$loglik + sum(dnorm(s$b, bm, 0.3, log = TRUE))
    s[c("a", "b")] <- mstep_items(e, s$a, s$b, priors = rtirt_priors(a = c(1, 100), b = c(0, 0.3)), b_mean = bm); s }
  fb <- run_em(list(Y = d$Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), stb, positive = "a")
  expect_lt(sum((fb$state$b - bm)^2), sum((fit$state$b - bm)^2))
  expect_error(mstep_items(estep_irt(d$Y, rep(1, 5), rep(0, 5)), rep(1, 5), rep(0, 5), b_mean = bm), "priors")
  Y <- d$Y; Y[sample(length(Y), 200)] <- NA
  expect_true(run_em(list(Y = Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), st, positive = "a")$converged)
  expect_error(run_em(list(Y = d$Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), function(s, x) s), "s\\$.loglik")
  e1 <- estep_irt(d$Y, fit$state$a, fit$state$b, mean = rep(0.3, 600)); e2 <- estep_irt(d$Y, fit$state$a, fit$state$b, mean = 0.3)
  expect_equal(e1$loglik, e2$loglik, tolerance = 1e-8)
})

test_that("check_sampler(step = ) passes a correct sampler and rejects a wrong prior", {
  skip_on_cran()
  N <- 12; K <- 3; pr <- rtirt_priors(a = c(1, 0.5), b = c(0, 1))
  prior <- function() list(a = draw_tnorm(rep(1, K), 0.5, 0, Inf), b = rnorm(K), theta = rnorm(N))
  simulate <- function(s) list(Y = matrix(rbinom(N * K, 1, plogis(sweep(outer(s$theta, s$a), 2, s$a * s$b))), N))
  mk <- function(sd) function(s, d) { s$omega <- draw_omega(s$a, s$b, s$theta); s[c("a", "b")] <- draw_items_pg(d$Y, s$omega, s$theta, s$a, pr)
    s$theta <- draw(irt_part(d$Y, s$omega, s$a, s$b) + normal_prior(0, sd)); s }
  expect_true(attr(check_sampler(step = mk(1), prior = prior, simulate = simulate, monitor = c("a", "b", "theta"), n_iter = 2e4, seed = 2), "passed"))
  expect_false(attr(check_sampler(step = mk(1.5), prior = prior, simulate = simulate, monitor = c("a", "b", "theta"), n_iter = 2e4, seed = 2), "passed"))
})

test_that("reg_part: the conditional of regression coefficients", {
  set.seed(5); X <- cbind(1, rnorm(50)); y <- drop(X %*% c(0.5, -1)) + rnorm(50); w <- runif(50, 0.5, 2)
  p <- reg_part(X, y, w) + normal_prior(0, 10)
  Q <- crossprod(X, w * X) + diag(0.01, 2); h <- crossprod(X, w * y)
  expect_equal(p$prec, Q, ignore_attr = TRUE); expect_equal(p$num, drop(h))
  G <- replicate(4000, draw(p)); expect_lt(max(abs(rowMeans(G) - solve(Q, h))), 0.02)
  y2 <- y; y2[1:5] <- NA; expect_equal(reg_part(X, y2)$prec, crossprod(X[-(1:5), ]), ignore_attr = TRUE)
  expect_error(reg_part(X, y) + normal_prior(rep(0, 3), 1), "cannot be added")
  expect_error(draw(reg_part(X[, c(1, 1)], y)), "positive definite")
  expect_output(print(p), "regression part of 2")
})

test_that("b_mean: per-item prior means of b in draw_items_pg() and draw_items_joint()", {
  d <- sim2pl_w(N = 60, K = 3, seed = 6); th <- rnorm(60); om <- draw_omega(d$a, d$b, th); pr <- rtirt_priors(b = c(0.4, 0.5))
  set.seed(1); r1 <- draw_items_pg(d$Y, om, th, d$a, pr); set.seed(1); r2 <- draw_items_pg(d$Y, om, th, d$a, pr, b_mean = rep(0.4, 3))
  expect_equal(r1, r2)
  set.seed(1); j1 <- draw_items_joint(d$Y, th, d$a, d$b, pr); set.seed(1); j2 <- draw_items_joint(d$Y, th, d$a, d$b, pr, b_mean = rep(0.4, 3))
  expect_equal(j1, j2)
  bm <- c(-1, 0, 2); kap <- d$Y - 0.5
  Q <- 1 / 0.25 + d$a^2 * colSums(om); h <- bm / 0.25 + d$a * colSums(sweep(om, 2, d$a, `*`) * th - kap)
  B <- replicate(4000, draw_items_pg(d$Y, om, th, d$a, pr, one_pl = TRUE, b_mean = bm)$b)
  expect_lt(max(abs(rowMeans(B) - h / Q) / sqrt(1 / Q)), 0.06)
  expect_error(draw_items_pg(d$Y, om, th, d$a, pr, b_mean = 1:2), "one prior mean per item")
})

test_that("print() of a sampler that has not mixed says so", {
  fit <- run_sampler(list(), function(d) list(x = 0), function(s, d) { s$x <- s$x + rnorm(1, 0, 0.01); s }, "x", n_iter = 400)
  expect_output(print(fit), "have not mixed")
})

test_that("run_em() warns about a wrong M-step and a .loglik at the wrong state", {
  set.seed(7); N <- 800; K <- 8; a0 <- runif(K, 0.8, 2); b0 <- rnorm(K)
  Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(rnorm(N, 0, 1.4), a0), 2, a0 * b0))), N)
  ini <- function(d) list(b = rep(0, K), sigma = 1)
  ll <- function(s, d) estep_irt(d$Y, rep(1, K), s$b, sd = s$sigma)$loglik
  mk <- function(sig) function(s, d) { e <- estep_irt(d$Y, rep(1, K), s$b, sd = s$sigma); s$.loglik <- e$loglik
    s$b <- mstep_items(e, rep(1, K), s$b, one_pl = TRUE)$b; s$sigma <- sig(e); s }
  ok <- expect_no_warning(run_em(list(Y = Y), ini, mk(function(e) sqrt(mean(e$eap^2 + e$psd^2))), ll, positive = "sigma"))
  w <- capture_warnings(run_em(list(Y = Y), ini, mk(function(e) 0.9 * sqrt(mean(e$eap^2 + e$psd^2))), ll, positive = "sigma"))
  expect_match(w, "not a maximum", all = FALSE)
  w <- capture_warnings(run_em(list(Y = Y), ini, mk(function(e) sd(e$eap)), ll, positive = "sigma"))
  expect_match(w, "decreased", all = FALSE); expect_match(w, "non-finite", all = FALSE)
  late <- function(s, d) { s <- mk(function(e) sqrt(mean(e$eap^2 + e$psd^2)))(s, d); s$.loglik <- ll(s, d); s }
  expect_warning(run_em(list(Y = Y), ini, late, ll, positive = "sigma"), "differs from loglik")
  expect_output(print(run_em(list(Y = Y), ini, mk(function(e) sqrt(mean(e$eap^2 + e$psd^2))), positive = "sigma")), "give run_em")
  e_only <- function(s, d) { s$.loglik <- ll(s, d); s }               # no M-step: the start is not a maximum
  w <- capture_warnings(run_em(list(Y = Y), ini, e_only, ll, positive = "sigma"))
  expect_match(w, "did not change any parameter", all = FALSE); expect_match(w, "not a maximum|not negative definite", all = FALSE)
})

test_that("runSimulation() is run_simulation() with the Julia defaults (same draws)", {
  c2 <- set_cond(n_subj = 40, n_item = 4, n_feat = 2, n_iter = 60, n_chain = 1, n_rep = 2)
  set.seed(1); p <- sim_para(c2, "cross")
  set.seed(2); a <- runSimulation(c2, p, setDataRtIrtCross, GibbsRtIrtCross, fit_indices = FALSE)
  set.seed(2); b <- run_simulation(c2, p, "cross", speed_var = "fixed", pars = c("a", "b", "theta"), fit_indices = FALSE)
  expect_identical(a, b)
  set.seed(1); p <- sim_para(c2, "latreg")
  set.seed(3); a <- runSimulation(c2, p, fit_indices = FALSE)                  # Julia default: setDataRtIrt + GibbsRtIrt
  set.seed(3); b <- run_simulation(c2, p, "latreg", speed_var = "fixed", pars = c("a", "b", "theta"), fit_indices = FALSE)
  expect_identical(a, b)
  expect_error(runSimulation(c2, p, setDataRtIrtCross, GibbsRtIrt), "matching pair")
})
