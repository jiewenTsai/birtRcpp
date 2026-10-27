# Regression snapshot (step 0 of the refactor).
# fixtures/make-regression.R stores what every case below returns; test-regression.R refits
# the cases and compares. Both files use the functions here, so the two sides cannot drift.

regression_fixture <- function() {
  test_path("fixtures", "regression.rds")
}

# ---- data --------------------------------------------------------------------

# One small simulated data set per built-in model, each with its own seed.
regression_data <- function(model) {
  seeds <- c(cross = 101, latent = 102, latreg = 103, null = 104, mlirt = 105)
  set.seed(seeds[[model]])
  needs_cov <- model %in% c("latent", "latreg", "mlirt")
  n_feat <- if (needs_cov) 2 else 0
  cond <- set_cond(n_subj = 200, n_item = 6, n_feat = n_feat)
  true_para <- sim_para(cond, model)
  sim_data(cond, true_para, model)
}

regression_2pl_matrix <- function() {
  set.seed(201)
  n_person <- 200
  n_item <- 6
  theta <- rnorm(n_person)
  a <- runif(n_item, 0.8, 2)
  b <- seq(-1.5, 1.5, length.out = n_item)
  eta <- sweep(outer(theta, a), 2, a * b)
  matrix(rbinom(n_person * n_item, 1, plogis(eta)), n_person)
}

regression_block_2pl <- function() {
  Y <- regression_2pl_matrix()
  block_model(Y = Y,
    theta[person] ~ normal(0, 1),
    a[item] ~ normal(1, 1, lower = 0),
    b[item] ~ normal(0, 3),
    Y ~ bernoulli_logit(a * (theta - b)))
}

regression_block_cross <- function() {
  sim <- regression_data("cross")
  block_model(Y = unname(sim$Y), logT = unname(sim$log_t),
    theta[person] ~ normal(0, 1),
    zeta[person] ~ normal(0, sqrt(v)),
    v ~ half_t(3, 1),
    a[item] ~ normal(1, 1, lower = 0),
    b[item] ~ normal(0, 3),
    lambda[item] ~ normal(0, 10),
    rho[item] ~ normal(0, 1),
    sigma2[item] ~ half_t(3, 1),
    Y ~ bernoulli_logit(a * (theta - b)),
    logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
}

# ---- cases -------------------------------------------------------------------

regression_constructors <- list(
  cross = rtirt_cross,
  latent = rtirt_latent,
  latreg = rtirt_latreg,
  null = rtirt_null,
  mlirt = mlirt)

regression_model <- function(model, ...) {
  dat <- regression_data(model)
  constructor <- regression_constructors[[model]]
  constructor(dat, ...)
}

regression_gibbs <- function(object, ...) {
  gibbs(object, n_iter = 400, n_chain = 2, seed = 11, verbose = FALSE, ...)
}

# the same fits through the interface of task 3.1: ecm(data, model = ...), gibbs(data, model = ...)
regression_names <- c(cross = "cross", latent = "latent", latreg = "latreg", null = "null", mlirt = "irt")

regression_gibbs_new <- function(model, collapse = FALSE, ...) {
  dat <- regression_data(model)
  settings <- gibbs_control(collapse = collapse, verbose = FALSE)
  gibbs(dat, model = regression_names[[model]], n_iter = 400, n_chain = 2, seed = 11, control = settings, ...)
}

# Each case is a function that returns a fitted object; `kind` says how to snapshot it.
regression_cases <- function() {
  cases <- list()
  for (model in names(regression_constructors)) {
    cases[[paste0("ecm_", model)]] <- local({
      m <- model
      list(kind = "ecm",
           fit = function() ecm(regression_model(m)),
           fit_new = function() ecm(regression_data(m), model = regression_names[[m]]))
    })
    cases[[paste0("gibbs_", model)]] <- local({
      m <- model
      list(kind = "gibbs",
           fit = function() regression_gibbs(regression_model(m)),
           fit_new = function() regression_gibbs_new(m))
    })
    cases[[paste0("gibbs_collapse_", model)]] <- local({
      m <- model
      list(kind = "gibbs",
           fit = function() regression_gibbs(regression_model(m), collapse = TRUE),
           fit_new = function() regression_gibbs_new(m, collapse = TRUE))
    })
  }
  cases$gibbs_quantile_latent <- list(kind = "gibbs", fit = function() {
    regression_gibbs(regression_model("latent", quantile = 0.25))
  }, fit_new = function() {
    regression_gibbs_new("latent", quantile = 0.25)
  })
  cases$gibbs_quantile_cross <- list(kind = "gibbs", fit = function() {
    regression_gibbs(regression_model("cross", quantile = 0.25))
  }, fit_new = function() {
    regression_gibbs_new("cross", quantile = 0.25)
  })
  cases$gibbs_init_ecm_cross <- list(kind = "gibbs", fit = function() {
    fit_ecm <- ecm(regression_model("cross"))
    regression_gibbs(fit_ecm, init = "ecm")
  }, fit_new = function() {
    fit_ecm <- ecm(regression_data("cross"), model = "cross")
    regression_gibbs_new("cross", init = fit_ecm)
  })
  cases$block_ecm_2pl <- list(kind = "block_ecm", fit = function() {
    ecm(regression_block_2pl())
  })
  cases$block_gibbs_2pl <- list(kind = "block_gibbs", fit = function() {
    model <- regression_block_2pl()
    steps <- list(step_mala(theta), step_gibbs(a, b))
    gibbs(model, steps, n_iter = 400, seed = 11)
  })
  cases$block_vi_2pl <- list(kind = "block_vi", fit = function() {
    vi(regression_block_2pl())
  })
  cases$block_ecm_cross <- list(kind = "block_ecm", fit = function() {
    model <- regression_block_cross()
    start <- list(lambda = colMeans(regression_data("cross")$log_t))
    ecm(model, prior = "ml", init = start)
  })
  cases
}

# ---- snapshots ---------------------------------------------------------------

# the hint lines ("Use summary(), ...") name the interface, which task 3.1 changed on purpose
drop_hints <- function(lines) {
  lines[!startsWith(lines, "Use ")]
}

# printed text, with run times removed
regression_text <- function(x) {
  old_locale <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", old_locale))
  Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
  lines <- utils::capture.output(print(x))
  gsub("[0-9]+(\\.[0-9]+)? s\\b", "<time> s", lines, perl = TRUE)
}

# the tables of summary(), as plain data frames
regression_tables <- function(summary_list) {
  is_table <- vapply(summary_list, is.data.frame, logical(1))
  tables <- summary_list[is_table]
  lapply(tables, as.data.frame)
}

snapshot_ecm <- function(fit) {
  list(
    estimates = as.data.frame(estimates(fit)),
    logLik = as.numeric(logLik(fit)),
    scores = as.data.frame(utils::head(scores(fit), 10)),
    tables = regression_tables(summary(fit)),
    equations = regression_text(equations(fit)),
    algorithm = regression_text(algorithm(fit)),
    print = regression_text(fit))
}

snapshot_gibbs <- function(fit) {
  draws <- as.matrix(fit$post$draws)
  est <- as.data.frame(estimates(fit))
  list(
    draws = draws,
    estimates = est,
    mcse = stats::setNames(est$sd / sqrt(est$ess), est$parameter),
    scores = as.data.frame(utils::head(scores(fit), 10)),
    tables = regression_tables(summary(fit)),
    equations = regression_text(equations(fit)),
    algorithm = regression_text(algorithm(fit)))
}

snapshot_block_ecm <- function(fit) {
  list(
    estimates = as.data.frame(fit$estimates),
    logLik = fit$loglik,
    algorithm = regression_text(algorithm(fit)),
    equations = regression_text(equations(fit$model)))
}

snapshot_block_gibbs <- function(fit) {
  draws <- as.matrix(fit$draws)
  est <- as.data.frame(summary(fit))
  list(
    draws = draws,
    estimates = est,
    mcse = stats::setNames(est$sd / sqrt(est$ess), est$parameter),
    algorithm = regression_text(algorithm(fit)))
}

snapshot_block_vi <- function(fit) {
  list(
    estimates = as.data.frame(fit$summary),
    elbo = utils::tail(fit$elbo, 1),
    algorithm = regression_text(algorithm(fit)))
}

# `interface = "new"` fits the case through ecm(data, model = ...) / gibbs(data, model = ...)
regression_snapshot <- function(case, interface = c("old", "new")) {
  interface <- match.arg(interface)
  fit <- if (interface == "new") case$fit_new() else case$fit()
  snapshot <- switch(case$kind,
    ecm = snapshot_ecm(fit),
    gibbs = snapshot_gibbs(fit),
    block_ecm = snapshot_block_ecm(fit),
    block_gibbs = snapshot_block_gibbs(fit),
    block_vi = snapshot_block_vi(fit))
  c(list(kind = case$kind), snapshot)
}

# ---- comparisons -------------------------------------------------------------

# largest absolute difference between the numeric columns of two tables
max_table_diff <- function(new, old) {
  is_num <- vapply(old, is.numeric, logical(1))
  diffs <- mapply(function(x, y) max(abs(x - y), 0, na.rm = TRUE), new[is_num], old[is_num])
  max(diffs, 0)
}

# the character columns must agree exactly
expect_same_labels <- function(new, old, label) {
  is_chr <- !vapply(old, is.numeric, logical(1))
  expect_identical(names(new), names(old), label = paste(label, "column names"))
  expect_identical(new[is_chr], old[is_chr], label = paste(label, "labels"))
}

expect_table_close <- function(new, old, tol, label) {
  expect_same_labels(new, old, label)
  expect_lte(max_table_diff(new, old), tol, label = paste(label, "max abs difference"))
}

# Same draws (same random numbers in the same order), or else posterior means within
# 3 Monte Carlo SEs (the MC error of both runs combined). With many parameters a few land
# outside 3 by chance, so up to 2% of them (at least one) may.
expect_draws_match <- function(new, old, label) {
  same_shape <- identical(dim(new$draws), dim(old$draws)) &&
    identical(colnames(new$draws), colnames(old$draws))
  if (same_shape && max(abs(new$draws - old$draws)) < 1e-10) {
    expect_table_close(new$estimates, old$estimates, 1e-10, paste(label, "estimates"))
    return(invisible("identical"))
  }
  common <- intersect(names(old$mcse), names(new$mcse))
  expect_setequal(names(new$mcse), names(old$mcse))
  mean_new <- stats::setNames(new$estimates[[2]], names(new$mcse))[common]
  mean_old <- stats::setNames(old$estimates[[2]], names(old$mcse))[common]
  mcse_both <- sqrt(new$mcse[common]^2 + old$mcse[common]^2)
  z <- abs(mean_new - mean_old) / mcse_both
  allowed <- max(1, floor(0.02 * length(z)))
  expect_lte(sum(z > 3), allowed, label = paste(label, "means beyond 3 MCSE"))
  invisible("monte carlo")
}
