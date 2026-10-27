# Simulation tools mirroring SimTools.jl. The generators draw from exactly the models the
# samplers fit (the Julia generators truncated log T at 0 and ignored the true sigma2t in the
# cross-relation data).
#
# The random numbers are drawn in a fixed order (item parameters a, b, sigma2t first, then the
# model-specific ones; covariates, persons, responses, log times), so a seed gives the same data
# as before the code was rewritten.

#' True parameters and simulated data
#'
#' `sim_para()` draws true parameters for a condition; `sim_data()` simulates data from them
#' (the person parameters are stored in `attr(data, "true_para")`). Item parameters:
#' a ~ N+(1, 0.2^2), b ~ N(0, 0.5^2), sigma2t ~ logN(log 0.3, 0.2^2).
#'
#' @param cond Conditions from [set_cond()].
#' @param model `"cross"` ([rtirt_cross()]), `"latent"` ([rtirt_latent()]), `"latreg"`
#'   ([rtirt_latreg()]), `"null"` ([rtirt_null()]) or `"mlirt"` ([mlirt()]).
#' @param para True parameters from `sim_para()` with the same `model`.
#' @param true_corr Correlation of ability and speed (`"latreg"`, `"null"`).
#' @param type Residual distribution: `"norm"`, `"tail"` (t with 5 df, scaled to the same
#'   variance) or `"skew"` (centred Gamma(0.5, 1), scaled to the same variance); `"cross"`
#'   (log RT residuals) and `"latent"` (speed residuals).
#' @return `sim_para`: an [input_para()]; `sim_data`: an [input_data()] (items `item01`, ...;
#'   covariates `x1`, ...) with the true parameters (item and person) in `attr(, "true_para")`,
#'   which the model constructors pick up.
#' @name simulation
NULL

sim_models <- c("cross", "latent", "latreg", "null", "mlirt")
models_with_cov <- c("mlirt", "latreg", "latent")

# ---- true parameters ---------------------------------------------------------------------

#' @rdname simulation
#' @export
sim_para <- function(cond, model = sim_models, true_corr = 0.3) {
  model <- match.arg(model)
  n_item <- cond$n_item
  n_feat <- cond$n_feat %||% 0L
  check_sim_cond(cond, model, n_feat)

  a <- abs(stats::rnorm(n_item, 1, 0.2))
  b <- stats::rnorm(n_item, 0, 0.5)
  sigma2t <- stats::rlnorm(n_item, log(0.3), 0.2)

  switch(model,
    mlirt = sim_para_mlirt(a, b, n_feat),
    latreg = sim_para_bivariate(a, b, sigma2t, n_feat, true_corr, with_cov = TRUE),
    null = sim_para_bivariate(a, b, sigma2t, n_feat, true_corr, with_cov = FALSE),
    latent = sim_para_latent(a, b, sigma2t, n_feat),
    cross = sim_para_cross(a, b, sigma2t))
}

check_sim_cond <- function(cond, model, n_feat) {
  if (is.null(cond$n_item) || is.null(cond$n_subj)) {
    stop("set n_subj and n_item in set_cond() for simulations")
  }
  if (model %in% models_with_cov && n_feat < 1) {
    stop("model '", model, "' needs n_feat >= 1")
  }
}

sim_para_mlirt <- function(a, b, n_feat) {
  beta <- stats::rnorm(n_feat, 0, 0.5)
  input_para(a = a, b = b, beta = beta)
}

# latreg and null: (ability, speed) bivariate normal given the covariates (latreg only)
sim_para_bivariate <- function(a, b, sigma2t, n_feat, true_corr, with_cov) {
  n_item <- length(a)
  lambda <- stats::rnorm(n_item, 4, 0.2)
  beta <- NULL
  if (with_cov) {
    beta_draws <- stats::rnorm(2 * n_feat, 0, 0.3)
    beta <- matrix(beta_draws, n_feat, 2, dimnames = list(NULL, c("theta", "speed")))
  }
  sd_speed <- sqrt(0.3)
  cov_ability_speed <- true_corr * sd_speed
  sigma_p <- matrix(c(1, cov_ability_speed, cov_ability_speed, sd_speed^2), 2)
  input_para(a = a, b = b, lambda = lambda, sigma2t = sigma2t, beta = beta, sigma_p = sigma_p)
}

sim_para_latent <- function(a, b, sigma2t, n_feat) {
  lambda <- stats::rnorm(length(a), 3, 0.2)
  beta <- stats::rnorm(n_feat, 0, 0.3)
  input_para(a = a, b = b, lambda = lambda, sigma2t = sigma2t, beta = beta, b_theta = -0.2,
             sigma_p = diag(c(1, 0.09)))
}

sim_para_cross <- function(a, b, sigma2t) {
  lambda <- stats::rnorm(length(a), 3, 0.5)
  rho <- stats::rnorm(length(a), 0, 0.2)
  input_para(a = a, b = b, lambda = lambda, sigma2t = sigma2t, rho = rho, sigma_p = diag(2))
}

# ---- data ----------------------------------------------------------------------------------

# standardized residuals: mean 0, variance 1
std_resid <- function(n, type) {
  switch(type,
    norm = stats::rnorm(n),
    tail = stats::rt(n, 5) / sqrt(5 / 3),
    skew = (stats::rgamma(n, 0.5, 1) - 0.5) / sqrt(0.5))
}

# 2PL responses: logit P(y = 1) = a (theta - b)
sim_y <- function(theta, a, b) {
  eta <- outer(theta, a) - matrix(a * b, length(theta), length(a), byrow = TRUE)
  prob <- stats::plogis(eta)
  matrix(stats::rbinom(length(prob), 1, prob), nrow(prob))
}

# log T = lambda - zeta - extra + e, e with variance sigma2t
log_rt <- function(lambda, zeta, sigma2t, extra = 0, type = "norm") {
  n_person <- length(zeta)
  n_item <- length(lambda)
  mean_log_t <- matrix(lambda, n_person, n_item, byrow = TRUE) - zeta - extra
  resid <- matrix(std_resid(n_person * n_item, type), n_person)
  resid_sd <- matrix(sqrt(sigma2t), n_person, n_item, byrow = TRUE)
  mean_log_t + resid * resid_sd
}

sim_cov <- function(n_person, n_feat) {
  draws <- stats::rnorm(n_person * n_feat)
  matrix(draws, n_person, dimnames = list(NULL, sprintf("x%d", seq_len(n_feat))))
}

# input_data() with the standard item names; the responses are simulated before the times
as_sim_data <- function(resp, time = NULL, cov = NULL) {
  colnames(resp) <- item_names(ncol(resp))
  suppressWarnings(input_data(resp, time, cov))
}

#' @rdname simulation
#' @export
sim_data <- function(cond, para, model = sim_models, type = c("norm", "tail", "skew")) {
  model <- match.arg(model)
  type <- match.arg(type)
  n_person <- cond$n_subj
  cov <- if (model %in% models_with_cov) sim_cov(n_person, cond$n_feat) else NULL

  sim <- switch(model,
    mlirt = sim_data_mlirt(para, cov),
    latreg = sim_data_bivariate(para, cov, n_person),
    null = sim_data_bivariate(para, cov, n_person),
    latent = sim_data_latent(para, cov, n_person, type),
    cross = sim_data_cross(para, n_person, type))
  data <- sim$data
  attr(data, "true_para") <- sim$para
  data
}

sim_data_mlirt <- function(para, cov) {
  n_person <- nrow(cov)
  para$theta <- drop(cov %*% para$beta) + stats::rnorm(n_person)
  resp <- sim_y(para$theta, para$a, para$b)
  list(data = as_sim_data(resp, cov = cov), para = para)
}

sim_data_bivariate <- function(para, cov, n_person) {
  person_mean <- if (is.null(cov)) matrix(0, n_person, 2) else cov %*% para$beta
  person_resid <- matrix(stats::rnorm(2 * n_person), n_person) %*% chol(para$sigma_p)
  para$theta <- person_mean[, 1] + person_resid[, 1]
  para$zeta <- person_mean[, 2] + person_resid[, 2]
  resp <- sim_y(para$theta, para$a, para$b)
  time <- exp(log_rt(para$lambda, para$zeta, para$sigma2t))
  list(data = as_sim_data(resp, time, cov), para = para)
}

sim_data_latent <- function(para, cov, n_person, type) {
  para$theta <- stats::rnorm(n_person)
  speed_resid <- sqrt(para$sigma_p[2, 2]) * std_resid(n_person, type)
  para$zeta <- drop(cov %*% para$beta) + para$b_theta * para$theta + speed_resid
  resp <- sim_y(para$theta, para$a, para$b)
  time <- exp(log_rt(para$lambda, para$zeta, para$sigma2t))
  list(data = as_sim_data(resp, time, cov), para = para)
}

sim_data_cross <- function(para, n_person, type) {
  para$theta <- stats::rnorm(n_person)
  para$zeta <- sqrt(para$sigma_p[2, 2]) * stats::rnorm(n_person)
  resp <- sim_y(para$theta, para$a, para$b)
  cross_effect <- outer(para$theta, para$rho)
  time <- exp(log_rt(para$lambda, para$zeta, para$sigma2t, cross_effect, type))
  list(data = as_sim_data(resp, time), para = para)
}

# ---- recovery ------------------------------------------------------------------------------

#' Bias, RMSE and correlation
#' @param est,true Estimates and true values.
#' @export
get_bias <- function(est, true) {
  mean(est - true)
}

#' @rdname get_bias
#' @export
get_rmse <- function(est, true) {
  sqrt(mean((est - true)^2))
}

#' @rdname get_bias
#' @export
get_corr <- function(est, true) {
  stats::cor(est, true)
}

#' Simulation study
#'
#' As `runSimulation()` in Julia: `cond$n_rep` replications with fixed item parameters and new
#' persons (and data) in each replication; returns bias, RMSE and correlation per parameter
#' block and replication.
#' @param cond Conditions.
#' @param para True parameters from [sim_para()].
#' @param model Model name as in [sim_para()].
#' @param quantile Quantile level of the ALD version (`"cross"`, `"latent"`), or `NULL`.
#' @param pars Parameter blocks to evaluate.
#' @param type Residual distribution of the simulated data.
#' @param speed_var Variance of speed (`"cross"`, `"latreg"`, `"null"`), as in the constructors.
#' @param ... Passed to [gibbs()].
#' @return A tibble with one row per replication and parameter block.
#' @export
run_simulation <- function(cond, para, model = sim_models, quantile = NULL, pars = c("a", "b", "ability"),
                           type = "norm", speed_var = c("free", "fixed"), ...) {
  model <- match.arg(model)
  speed_var <- match.arg(speed_var)
  replications <- seq_len(cond$n_rep)

  one_replication <- function(rep) {
    data <- sim_data(cond, para, model, type)
    fit_model <- sim_model(model, data, cond, quantile, speed_var)
    fit <- gibbs_rtirt(fit_model, verbose = FALSE, ...)
    recovery <- purrr::map(pars, function(par) compare_para(fit, par))
    dplyr::bind_rows(recovery) |>
      dplyr::mutate(rep = rep, .before = 1)
  }

  purrr::map(replications, one_replication) |>
    dplyr::bind_rows()
}

# the model object of run_simulation() (without the deprecation warnings of the constructors)
sim_model <- function(model, data, cond, quantile, speed_var) {
  switch(model,
    cross = new_model("rtirt_cross", data, cond, quantile = quantile, speed_var = speed_var),
    latent = new_model("rtirt_latent", data, cond, quantile = quantile),
    latreg = new_model("rtirt_latreg", data, cond, speed_var = speed_var),
    null = new_model("rtirt_null", data, cond, speed_var = speed_var),
    mlirt = new_model("mlirt", data, cond))
}
