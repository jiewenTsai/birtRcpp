# ECM of the quantile (ALD) models: designs in R/em_quantile.R

quantile_data <- function(model = "cross", n = 150, k = 5, seed = 1) {
  set.seed(seed)
  cond <- set_cond(n_subj = n, n_item = k, n_feat = 2)
  sim_data(cond, sim_para(cond, model), model)
}

# context of the internal functions, at starting values moved away from the fit
quantile_context <- function(model, itemtype = "2pl", q = 0.3, speed_var = "free") {
  dat <- quantile_data(model)
  object <- build_model(dat, model, itemtype, speed_var, NULL, q)
  E <- em_setup(object, 21, FALSE)
  E$q <- q
  x <- em_pack(em_start(E), E)
  set.seed(2)
  list(E = E, x = x + stats::rnorm(length(x), 0, 0.2))
}

numerical_gradient <- function(f, x) {
  vapply(seq_along(x), function(k) {
    h <- 1e-5 * max(1, abs(x[k]))
    (f(replace(x, k, x[k] + h)) - f(replace(x, k, x[k] - h))) / (2 * h)
  }, 0)
}

test_that("the speed integral of \"cross\" has the derivatives and moments of its closed form", {
  set.seed(1)
  q <- 0.3
  B <- matrix(stats::rnorm(4 * 5, 0, 1.5), 4)
  s <- stats::runif(5, 0.3, 1)
  v <- 0.7
  g <- .ald_rt_grad(B, s, v, q)
  expect_equal(g$ll, .log_ald_rt(B, s, v, q))
  h <- 1e-6
  for (k in 1:5) {
    e <- matrix(0, nrow(B), 5)
    e[, k] <- h
    expect_equal(g$gb[, k], (.log_ald_rt(B + e, s, v, q) - .log_ald_rt(B - e, s, v, q)) / (2 * h), tolerance = 1e-6)
    d <- replace(numeric(5), k, h)
    expect_equal(g$gs[, k], (.log_ald_rt(B, s + d, v, q) - .log_ald_rt(B, s - d, v, q)) / (2 * h), tolerance = 1e-6)
  }
  expect_equal(g$gv, (.log_ald_rt(B, s, v + h, q) - .log_ald_rt(B, s, v - h, q)) / (2 * h), tolerance = 1e-6)
  # posterior moments of speed for row 1, by numerical integration between the kinks
  log_ald <- function(u) log(q * (1 - q) / s) - u * (q - (u < 0)) / s
  dens <- function(z) vapply(z, function(zz) exp(sum(log_ald(B[1, ] + zz)) + stats::dnorm(zz, 0, sqrt(v), log = TRUE)), 0)
  cuts <- c(-Inf, sort(-B[1, ]), Inf)
  piece <- function(f) sum(vapply(seq_len(length(cuts) - 1), function(j) stats::integrate(f, cuts[j], cuts[j + 1], rel.tol = 1e-10)$value, 0))
  total <- piece(dens)
  expect_equal(log(total), g$ll[1], tolerance = 1e-8)
  expect_equal(piece(function(z) z * dens(z)) / total, g$ez[1], tolerance = 1e-6)
  expect_equal(piece(function(z) z^2 * dens(z)) / total, g$ez2[1], tolerance = 1e-6)
})

test_that("the analytic score is the gradient of the working log-likelihood", {
  for (model in c("cross", "latent")) {
    for (itemtype in c("2pl", "1pl")) {
      ctx <- quantile_context(model, itemtype)
      loglik <- function(x) sum(em_qr_estep(em_unpack(x, ctx$E), ctx$E)$ll)
      numerical <- numerical_gradient(loglik, ctx$x)
      expect_equal(em_qr_score(ctx$x, ctx$E), numerical, tolerance = 1e-6, info = paste(model, itemtype))
    }
  }
  ctx <- quantile_context("cross", speed_var = "fixed")
  loglik <- function(x) sum(em_qr_estep(em_unpack(x, ctx$E), ctx$E)$ll)
  expect_equal(em_qr_score(ctx$x, ctx$E), numerical_gradient(loglik, ctx$x), tolerance = 1e-6)
})

test_that("speed of \"latent\" a posteriori: the moments of the two truncated normals", {
  ctx <- quantile_context("latent", q = 0.75)
  E <- ctx$E
  p <- em_unpack(ctx$x, E)
  TH <- matrix(E$z, E$N, length(E$z), byrow = TRUE)
  lat <- em_qr_latent(p, E, TH)
  centre <- lat$zh[1]
  mean_speed <- (em_xb(E$X, p$beta) + p$gamma * TH)[1, 4]
  q <- E$q
  dens <- function(z) stats::dnorm(z, centre, sqrt(1 / lat$A)) * q * (1 - q) / p$s * exp(-(z - mean_speed) * (q - (z < mean_speed)) / p$s)
  part <- function(f) stats::integrate(f, -Inf, mean_speed, rel.tol = 1e-12)$value + stats::integrate(f, mean_speed, Inf, rel.tol = 1e-12)$value
  total <- part(dens)
  expect_equal(lat$F[1, 4], log(total), tolerance = 1e-8)
  expect_equal(lat$Ez[1, 4], part(function(z) z * dens(z)) / total, tolerance = 1e-8)
  expect_equal(lat$Ez2[1, 4], part(function(z) z^2 * dens(z)) / total, tolerance = 1e-8)
})

test_that("ecm() fits quantile models: monotone, a maximum, sandwich standard errors", {
  dat <- quantile_data("cross", n = 200)
  fit <- ecm(dat, model = "cross", quantile = 0.7)
  expect_true(fit$ecm$converged)
  expect_true(all(diff(fit$ecm$trace) > -1e-8 * abs(fit$ecm$trace[-1])))
  expect_identical(fit$ecm$se_type, "sandwich")
  e <- estimates(fit)
  expect_true(all(is.finite(e$se) & e$se > 0))
  expect_true(isSymmetric(unname(fit$ecm$vcov)))
  checks <- check_ecm(fit, n_start = 2, seed = 1)
  expect_true(all(checks$passed))
  expect_output(print(fit), "maximum marginal ALD working likelihood")
  expect_output(print(fit), "Standard errors: sandwich")
  expect_output(print(summary(fit)), "ALD scale of log RT \\(sigma2t\\)")
  expect_output(print(summary(fit)), "ALD working-likelihood estimate, sandwich SE")
  # the q = 0.5 fit and the mean fit on normal residuals: the same cross-relations, about
  rho_median <- estimates(ecm(dat, model = "cross", quantile = 0.5, se = FALSE), pars = "^rho")$est
  rho_mean <- estimates(ecm(dat, model = "cross", se = FALSE), pars = "^rho")$est
  expect_lt(max(abs(rho_median - rho_mean)), 0.15)
})

test_that("quantile ECM fits: several levels, latent covariates, gibbs() starts, anova() refuses", {
  dat <- quantile_data("latent", n = 200, k = 4)
  qs <- ecm(dat, model = "latent", quantile = c(0.25, 0.75), cov = "x1", se = FALSE)
  expect_s3_class(qs, "rtirt_qset")
  expect_output(print(qs), "q=0.25")
  e <- estimates(qs)
  expect_true(all(c(0.25, 0.75) %in% e$q))
  expect_true(all(vapply(qs$models, has_em, TRUE)))
  fit <- qs$models[[2]]
  g <- gibbs(dat, model = "latent", quantile = 0.75, cov = "x1", init = fit, n_iter = 60, n_chain = 1, seed = 1,
             control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_false(is.null(g$post))
  g2 <- gibbs(dat, model = "latent", quantile = c(0.25, 0.75), cov = "x1", init = "ecm", n_iter = 60, n_chain = 1, seed = 1,
              control = gibbs_control(verbose = FALSE, fit_indices = FALSE))
  expect_s3_class(g2, "rtirt_qset")
  expect_error(gibbs(dat, model = "latent", quantile = 0.25, cov = "x1", init = fit, n_iter = 60),
               "init was fitted with quantile = 0.75, gibbs\\(\\) was asked for quantile = 0.25")
  expect_error(anova(fit, ecm(dat, model = "latent", cov = "x1", se = FALSE)), "ALD working likelihood, not comparable")
  expect_output(print(algorithm(fit)), "sandwich")
})
