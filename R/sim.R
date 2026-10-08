# Simulation tools mirroring SimTools.jl. The generators draw from exactly the models the
# samplers fit (the Julia generators truncated log T at 0 and ignored the true sigma2t in the
# cross-relation data).

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

#' @rdname simulation
#' @export
sim_para <- function(cond, model = sim_models, true_corr = 0.3) {
  model <- match.arg(model); K <- cond$n_item; P <- cond$n_feat %||% 0L
  if (is.null(K) || is.null(cond$n_subj)) stop("set n_subj and n_item in set_cond() for simulations")
  if (model %in% c("mlirt", "latreg", "latent") && P < 1) stop("model '", model, "' needs n_feat >= 1")
  a <- abs(stats::rnorm(K, 1, 0.2)); b <- stats::rnorm(K, 0, 0.5); s2 <- stats::rlnorm(K, log(0.3), 0.2)
  switch(model,
    mlirt = input_para(a = a, b = b, beta = stats::rnorm(P, 0, 0.5)),
    latreg = , null = {
      sd2 <- sqrt(0.3)                               # SD of speed given the covariates
      input_para(a = a, b = b, lambda = stats::rnorm(K, 4, 0.2), sigma2t = s2,
                 beta = if (model == "latreg") matrix(stats::rnorm(2 * P, 0, 0.3), P, 2, dimnames = list(NULL, c("theta", "speed"))),
                 sigma_p = matrix(c(1, true_corr * sd2, true_corr * sd2, sd2^2), 2)) },
    latent = input_para(a = a, b = b, lambda = stats::rnorm(K, 3, 0.2), sigma2t = s2,
                        beta = stats::rnorm(P, 0, 0.3), b_theta = -0.2, sigma_p = diag(c(1, 0.09))),
    cross = input_para(a = a, b = b, lambda = stats::rnorm(K, 3, 0.5), sigma2t = s2,
                       rho = stats::rnorm(K, 0, 0.2), sigma_p = diag(2)))
}

std_resid <- function(n, type) switch(type, norm = stats::rnorm(n), tail = stats::rt(n, 5) / sqrt(5 / 3),
                                      skew = (stats::rgamma(n, 0.5, 1) - 0.5) / sqrt(0.5))
sim_y <- function(theta, a, b) {
  P <- stats::plogis(outer(theta, a) - matrix(a * b, length(theta), length(a), byrow = TRUE))
  matrix(stats::rbinom(length(P), 1, P), nrow(P))
}
log_rt <- function(lambda, zeta, sigma2t, extra = 0, type = "norm") {
  N <- length(zeta); K <- length(lambda)
  matrix(lambda, N, K, byrow = TRUE) - zeta - extra +
    matrix(std_resid(N * K, type), N) * matrix(sqrt(sigma2t), N, K, byrow = TRUE)
}

#' @rdname simulation
#' @export
sim_data <- function(cond, para, model = sim_models, type = c("norm", "tail", "skew")) {
  model <- match.arg(model); type <- match.arg(type); N <- cond$n_subj
  X <- if (model %in% c("mlirt", "latreg", "latent")) matrix(stats::rnorm(N * cond$n_feat), N, dimnames = list(NULL, sprintf("x%d", seq_len(cond$n_feat)))) else NULL
  inp <- function(Y, T = NULL, X = NULL) { colnames(Y) <- item_names(ncol(Y)); suppressWarnings(input_data(Y, T, X)) }
  switch(model,
    mlirt = { para$theta <- drop(X %*% para$beta) + stats::rnorm(N); d <- inp(sim_y(para$theta, para$a, para$b), X = X) },
    latreg = , null = {
      mu <- if (is.null(X)) matrix(0, N, 2) else X %*% para$beta
      E <- matrix(stats::rnorm(2 * N), N) %*% chol(para$sigma_p)
      para$theta <- mu[, 1] + E[, 1]; para$zeta <- mu[, 2] + E[, 2]
      d <- inp(sim_y(para$theta, para$a, para$b), exp(log_rt(para$lambda, para$zeta, para$sigma2t)), X) },
    latent = {
      para$theta <- stats::rnorm(N)
      para$zeta <- drop(X %*% para$beta) + para$b_theta * para$theta + sqrt(para$sigma_p[2, 2]) * std_resid(N, type)
      d <- inp(sim_y(para$theta, para$a, para$b), exp(log_rt(para$lambda, para$zeta, para$sigma2t)), X) },
    cross = {
      para$theta <- stats::rnorm(N); para$zeta <- sqrt(para$sigma_p[2, 2]) * stats::rnorm(N)
      d <- inp(sim_y(para$theta, para$a, para$b),
                      exp(log_rt(para$lambda, para$zeta, para$sigma2t, outer(para$theta, para$rho), type))) })
  attr(d, "true_para") <- para
  d
}

#' Bias, RMSE and correlation
#' @param est,true Estimates and true values.
#' @export
get_bias <- function(est, true) mean(est - true)
#' @rdname get_bias
#' @export
get_rmse <- function(est, true) sqrt(mean((est - true)^2))
#' @rdname get_bias
#' @export
get_corr <- function(est, true) stats::cor(est, true)

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
#' @export
run_simulation <- function(cond, para, model = sim_models, quantile = NULL, pars = c("a", "b", "ability"),
                           type = "norm", speed_var = c("free", "fixed"), ...) {
  model <- match.arg(model); speed_var <- match.arg(speed_var)
  ctor <- switch(model, cross = function(d) rtirt_cross(d, quantile, speed_var = speed_var, cond = cond),
                 latent = function(d) rtirt_latent(d, quantile, cond = cond),
                 latreg = function(d) rtirt_latreg(d, speed_var = speed_var, cond = cond),
                 null = function(d) rtirt_null(d, speed_var = speed_var, cond = cond), mlirt = function(d) mlirt(d, cond = cond))
  out <- lapply(seq_len(cond$n_rep), function(r) {
    d <- sim_data(cond, para, model, type)
    fit <- gibbs(ctor(d), verbose = FALSE, ...)
    cbind(rep = r, do.call(rbind, lapply(pars, function(p) compare_para(fit, p))))
  })
  do.call(rbind, out)
}
