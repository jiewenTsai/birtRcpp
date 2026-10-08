fit_em <- function(type, n = 300, k = 6, ...) {
  cond <- set_cond(n_subj = n, n_item = k, n_feat = 2)
  dat <- sim_data(cond, sim_para(cond, type), type)
  ctor <- switch(type, mlirt = mlirt, null = rtirt_null, latreg = rtirt_latreg, latent = rtirt_latent, cross = rtirt_cross)
  ctor(dat, ...)
}

test_that("EM is monotone, its logLik is the marginal likelihood, and the score is the gradient", {
  set.seed(11)
  for (type in c("mlirt", "null", "latreg", "latent", "cross")) {
    m <- fit_em(type)
    f0 <- ecm(m, accelerate = FALSE, se = FALSE, max_iter = 300)
    expect_gte(min(diff(f0$ecm$trace)), -1e-8)
    f <- ecm(m, se = FALSE)
    expect_true(f$ecm$converged)
    expect_gte(f$ecm$logLik, f0$ecm$logLik - 1e-6)
    # logLik equals the marginal_loglik() code at the ECM estimates
    e <- estimates(f); ff <- f; ff$post <- list(settings = list(X = f$ecm$X, q = NA))
    p <- birtRcpp:::draw_pars(stats::setNames(e$est, e$parameter), f$data$items, colnames(f$ecm$X) %||% character())
    expect_equal(sum(birtRcpp:::marg_ll_draw(ff, p, 41)), f$ecm$logLik, tolerance = 1e-8)
    # analytic score vs numerical gradient, away from the optimum
    E <- birtRcpp:::em_setup(m, 41, FALSE); x <- f$ecm$x + stats::rnorm(length(f$ecm$x), 0, 0.05)
    obj <- function(x) sum(birtRcpp:::em_estep(birtRcpp:::em_unpack(x, E), E)$ll)
    ng <- vapply(seq_along(x), function(j) { h <- 1e-5; d <- replace(numeric(length(x)), j, h); (obj(x + d) - obj(x - d)) / (2 * h) }, 0)
    expect_lt(max(abs(ng - birtRcpp:::em_score(x, E))), 1e-4)
  }
})

test_that("EM recovers the parameters and gives standard errors", {
  set.seed(12)
  m <- fit_em("cross", n = 1500, k = 8)
  f <- ecm(m); e <- estimates(f); tp <- attr(m$data, "true_para")
  expect_lt(max(abs(e$est[startsWith(e$parameter, "a[")] - tp$a)), 0.35)
  expect_lt(max(abs(e$est[startsWith(e$parameter, "rho[")] - tp$rho)), 0.2)
  expect_true(all(e$se > 0))
  z <- (e$est - c(tp$a, tp$b, tp$lambda, tp$sigma2t, tp$rho, tp$sigma_p[2, 2])) / e$se
  expect_lt(mean(abs(z) > 2.58), 0.1)
  expect_equal(dim(vcov(f)), c(nrow(e), nrow(e)))
  expect_s3_class(logLik(f), "logLik"); expect_equal(AIC(f), f$ecm$AIC)
  expect_named(fit_indices(f), c("method", "logLik", "npar", "N", "AIC", "BIC"))
  expect_true(all(c("ability", "ability_psd", "speed", "speed_psd") %in% names(scores(f))))
  s <- summary(f); expect_output(print(s), "Structural:"); expect_true(all(c("z", "pvalue") %in% names(as.data.frame(s))))
})

test_that("EM options: fixed speed variance, 1PL, item priors; quantile models refused", {
  set.seed(13)
  f <- ecm(fit_em("latreg", speed_var = "fixed"), se = FALSE)
  e <- estimates(f); expect_equal(e$est[e$parameter == "var_speed"], 1)
  f <- ecm(fit_em("null", itemtype = "1pl"))
  e <- estimates(f); expect_true(all(e$est[startsWith(e$parameter, "a[")] == 1)); expect_true(all(is.na(e$se[startsWith(e$parameter, "a[")])))
  m <- fit_em("cross", n = 150); fm <- ecm(m, prior = "map", se = FALSE); fl <- ecm(m, se = FALSE)
  expect_lt(fm$ecm$logLik, fl$ecm$logLik + 1e-8)                # MAP is not the likelihood maximum
  expect_error(ecm(fit_em("cross", quantile = 0.3)), "quantile models use gibbs")
  expect_error(ecm(fit_em("latent", quantile = c(0.3, 0.5))), "quantile models use gibbs")
})

test_that("Gibbs chains can start from the ECM fit", {
  set.seed(14)
  m <- fit_em("latent")
  g <- gibbs(m, n_iter = 200, n_chain = 2, init = "ecm", seed = 1, verbose = FALSE, fit_indices = FALSE)
  expect_false(is.null(g$ecm)); expect_false(is.null(g$post))
  expect_true("rhat" %in% names(estimates(g))); expect_true("se" %in% names(estimates(g, method = "ecm")))
})

test_that("ECM fits: labels, MAP without AIC/BIC, convergence(), compare_para() and plot() work", {
  cond <- set_cond(n_subj = 200, n_item = 5)
  d <- sim_data(cond, sim_para(cond, "cross"), "cross")
  f <- ecm(rtirt_cross(d))
  expect_match(paste(capture.output(print(f)), collapse = "\n"), "ECM, 41 Gauss-Hermite nodes")
  expect_match(paste(capture.output(print(f)), collapse = "\n"), 'init = "ecm"', fixed = TRUE)
  expect_equal(fit_indices(f)$method, "ECM"); expect_equal(f$ecm$nodes, 41)
  expect_output(cv <- convergence(f), "ECM: converged"); expect_true(cv$converged)
  cp <- compare_para(f, "a"); expect_equal(cp$n, 5); expect_lt(abs(cp$bias), 1)
  pdf(NULL); expect_silent(plot(f)); dev.off()
  expect_error(plot(f, type = "trace"), "gibbs")
  fm <- ecm(rtirt_cross(d), prior = "map", se = FALSE)
  expect_true(is.na(fit_indices(fm)$AIC)); expect_match(paste(capture.output(print(fm)), collapse = "\n"), "no AIC / BIC")
})
