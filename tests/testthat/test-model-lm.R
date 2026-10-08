sim_dif <- function(N = 800, K = 6, dif = 0, impact = 0, seed = 1) {
  set.seed(seed); z <- rbinom(N, 1, 0.5); th <- rnorm(N, impact * z); a0 <- runif(K, 0.8, 1.8); b0 <- seq(-1, 1, length.out = K)
  eta <- sweep(outer(th, a0), 2, a0 * b0); eta[, 2] <- eta[, 2] + dif * z
  list(Y = matrix(rbinom(N * K, 1, plogis(eta)), N), z = z, K = K)
}

test_that("score_test(add = ) equals the likelihood ratio statistic asymptotically (Glas, 1998, 1999)", {
  d <- sim_dif(dif = 0.6, seed = 2); Y <- d$Y; z <- d$z; K <- d$K
  m0 <- block_model(Y = Y, z = z, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
                    Y ~ bernoulli_logit(a * (theta - b)))
  f0 <- ecm(m0, prior = "ml", se = FALSE)
  st <- score_test(f0, add = Y ~ delta[item] * z)
  expect_s3_class(st, "block_lm_test"); expect_equal(st$test, c("all", sprintf("delta[%d]", 1:K))); expect_equal(st$df, rep(1:0, c(1, K)) * (K - 1) + 1)
  m1 <- block_model(Y = Y, z = z, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), delta[item] ~ normal(0, 1),
                    Y ~ bernoulli_logit(a * (theta - b) + delta * z))
  lr <- 2 * (ecm(m1, prior = "ml", se = FALSE)$loglik - f0$loglik)
  expect_lt(abs(st$statistic[1] - lr) / lr, 0.03)
  zi <- outer(z, seq_len(K) == 2) * 1                                         # item 2 alone: a scalar new parameter and a data matrix in data =
  m2 <- block_model(Y = Y, zi = zi, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), d ~ normal(0, 1),
                    Y ~ bernoulli_logit(a * (theta - b) + d * zi))
  f2 <- ecm(m2, prior = "ml", se = FALSE)
  s2 <- score_test(f0, add = Y ~ d * zi, data = list(zi = zi))
  expect_equal(s2$statistic[1], st$statistic[3], tolerance = 1e-8)            # the per-item test is the test of that item alone
  expect_lt(abs(s2$statistic[1] - 2 * (f2$loglik - f0$loglik)), 0.05)
  expect_lt(abs(st$estimate[3] - f2$state$d), 0.02)                           # one-step estimate (expected parameter change)
  expect_lt(st$p.value[3], 0.01)
  expect_output(print(st), "Lagrange multiplier")
  expect_equal(score_test(mml(m0, prior = "ml", se = FALSE), add = Y ~ delta[item] * z)$statistic, st$statistic, tolerance = 1e-4)
})

test_that("score_test(add = ) refuses what it cannot test, and handles absorbed directions", {
  d <- sim_dif(N = 400, impact = 0.5, seed = 3); Y <- d$Y; z <- d$z; K <- d$K
  m0 <- block_model(Y = Y, z = z, theta[person] ~ normal(g * z, 1), g ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
                    Y ~ bernoulli_logit(a * (theta - b)))
  f0 <- ecm(m0, prior = "ml", se = FALSE)
  st <- score_test(f0, add = Y ~ delta[item] * z)
  expect_equal(st$df[1], K - 1)                                               # a shift of every item along a is absorbed by g
  expect_output(print(st), "absorbed by the null parameters")
  zz <- outer(z, c(0, rep(1, K - 1)))                                         # item 1 without the term: no information, and the rest identified
  s2 <- score_test(f0, add = Y ~ delta[item] * zz, data = list(zz = zz))
  expect_true(is.na(s2$statistic[2])); expect_equal(s2$df[1], K - 1); expect_lt(abs(s2$statistic[1] - st$statistic[1]), 0.5)      # two 7-df tests of the same hypothesis
  expect_error(score_test(f0), "give the term")
  expect_error(score_test(f0, add = X ~ delta[item] * z), "not observed data")
  expect_error(score_test(f0, add = Y ~ a[item] * z), "already a name")
  expect_error(score_test(f0, add = Y ~ d * w), "multiplies the new parameters")
  expect_error(score_test(f0, add = Y ~ delta[item] * w, data = list(w = c(NA, rep(1, 399)))), "missing values")
  expect_error(score_test(f0, add = Y ~ exp(delta) * z), "must be 0")
  expect_error(score_test(f0, add = Y ~ delta[group] * z), "dimension")
})

test_that("score_test(add = ) reproduces ci_test() of the built-in model (a normal part)", {
  skip_on_cran()
  set.seed(4); cond <- set_cond(n_subj = 500, n_item = 6); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  ref <- ecm(rtirt_cross(sim), se = FALSE)
  m <- block_model(Y = unname(sim$Y), logT = unname(sim$log_t),
    theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, sqrt(v)), v ~ half_t(3, 1),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
    sigma2[item] ~ half_t(3, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  f <- ecm(m, prior = "ml", se = FALSE)
  ci <- ci_test(ref); st <- score_test(f, add = logT ~ delta[item] * Y)        # log times of correct responses shifted
  expect_lt(max(abs(st$statistic[-1] - ci$statistic)), 0.01)
  expect_lt(max(abs(st$estimate[-1] - ci$delta)), 1e-3)
  expect_lt(abs(st$statistic[1] - attr(ci, "overall")[["statistic"]]), 0.01)
})

test_that("Louis' information with person-level parameters that are not integrated out equals the numerical derivative of the score", {
  set.seed(5); cond <- set_cond(n_subj = 150, n_item = 5); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  LT <- unname(sim$log_t); LT[sample(length(LT), 15)] <- NA
  m <- block_model(Y = unname(sim$Y), logT = LT, theta[person] ~ normal(c * zeta, 1), zeta[person] ~ normal(0, 1), c ~ normal(0, 1),
    u[person] ~ normal(0, 0.5), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
    sigma2[item] ~ inv_gamma(2, 1), Y ~ bernoulli_logit(a * (theta - b) + 0.3 * zeta + u), logT ~ normal(lambda - zeta - rho * theta - 0.5 * u, sqrt(sigma2)))
  for (lat in list("theta", c("theta", "u"))) {                               # zeta (and u) fixed per person; with u integrated exactly
    fit <- ecm(m, latent = lat, se = FALSE)
    E <- bm_setup(m, NULL, lat, NULL, "map", NULL); s <- fit$state[E$fixed]; C <- fit$ecm$C
    sc <- bm_score(E, C); g <- function(x) unlist(sc(utils::relist(x, s), NULL)[E$fixed])
    Hn <- calc_hessian(g, unlist(s), h = 1e-5); Hl <- bm_information(E, C)(s, NULL)
    expect_lt(max(abs(Hl - Hn)) / max(abs(Hn)), 1e-5)
  }
  fit <- ecm(m, latent = "theta")
  expect_true(all(is.finite(fit$estimates$se)))
})

test_that("algorithm() writes ECM as the objective, the E-step and the conditional modes; Newton as one maximization", {
  set.seed(6); cond <- set_cond(n_subj = 100, n_item = 4); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  m <- block_model(Y = unname(sim$Y), logT = unname(sim$log_t), theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, 1),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
    sigma2[item] ~ inv_gamma(2, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  f <- ecm(m, se = FALSE)
  u <- as.character(algorithm(f)); l <- paste(as.character(algorithm(f, format = "latex")), collapse = "\n")
  expect_true(all(c("Model p:", "Objective:", "E-step:", "M-step: conditional maximizations, in this order:", "Then:") %in% u))
  expect_match(l, "\\lambda_{i} &\\leftarrow \\frac{m_{i}}{p_{i}}", fixed = TRUE)
  expect_match(l, "m_{i} &= \\sum_{j} \\frac{\\mathbb{E}\\left[\\mathrm{logT}_{ji} + \\zeta_{j} + \\rho_{i}\\,\\theta_{j}\\right]}{\\sigma_{i}^2}", fixed = TRUE)
  expect_match(l, "\\sigma_{i}^2 &\\leftarrow \\frac{1 + \\tfrac{1}{2}", fixed = TRUE)
  expect_match(l, "\\zeta_{j} \\mid \\theta_{j}, y_{j} &\\sim \\mathcal{N}", fixed = TRUE)
  expect_match(paste(as.character(algorithm(ecm(m, prior = "ml", se = FALSE, quadrature = "ba81"), format = "latex")), collapse = "\n"),
               "\\frac{\\left(\\sum_{j} \\mathbb{E}", fixed = TRUE)                # the ML mode of a variance: the mean squared residual
  n <- as.character(algorithm(mml(m, se = FALSE)))
  expect_true("Maximization:" %in% n); expect_false(any(grepl("conditional maximizations", n)))
})

test_that("vi(quadrature = \"ba81\") and the 2L points for an exact latent: the same fixed point, the ELBO never decreases", {
  set.seed(8); cond <- set_cond(n_subj = 200, n_item = 5); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  m <- block_model(Y = unname(sim$Y), logT = unname(sim$log_t), theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, sqrt(v)), v ~ half_t(3, 1),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
    sigma2[item] ~ half_t(3, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  E <- vb_setup(bm_setup(m, NULL, NULL, NULL, "map", NULL, TRUE, from_data = FALSE)); expect_true(E$sigma)
  E2 <- vb_setup(bm_setup(m, NULL, NULL, NULL, "map", NULL, TRUE, sigma = FALSE, from_data = FALSE))   # points m -/+ sd at every node
  E2$affine <- FALSE                                                                                # and the exact step node by node
  set.seed(1); Q <- vb_init(E); C <- bm_centres(E, E$s, NULL); Q2 <- Q
  for (k in 1:3) {
    e <- vb_estep(E, Q, C); e2 <- vb_estep(E2, Q2, C)
    expect_equal(e$ll, e2$ll, tolerance = 1e-10)
    for (st in E$steps) for (v in st$vars) { Q <- vb_update(E, Q, e, v); Q2 <- vb_update(E2, Q2, e2, v) }
  }
  expect_equal(unlist(Q$norm), unlist(Q2$norm), tolerance = 1e-10); expect_equal(unlist(Q$var), unlist(Q2$var), tolerance = 1e-10)
  va <- vi(m); vb <- vi(m, quadrature = "ba81")
  expect_equal(unname(vb$elbo_drop[["worst"]]), "0")
  expect_lt(max(abs(vb$summary$mean - va$summary$mean)), 0.01)
  expect_lt(abs(vb$elbo[length(vb$elbo)] - va$elbo[length(va$elbo)]), 0.05)
  expect_true(attr(check_vi(m, n_start = 2, sweeps = 6, seed = 1, quadrature = "ba81"), "passed"))
  al <- paste(as.character(algorithm(vb)), collapse = "\n")
  expect_match(al, "common node", fixed = TRUE); expect_match(al, "summed over persons at each node", fixed = TRUE)
  expect_match(paste(as.character(algorithm(va)), collapse = "\n"), "through 4 points", fixed = TRUE)
})
