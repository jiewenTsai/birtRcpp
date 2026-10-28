# The built-in models by name, as in mirt and TAM:
#   fit_ecm <- ecm(dat, model = "cross")
#   fit_gibbs <- gibbs(dat, model = "cross", init = fit_ecm)
# `model` picks the structure and `itemtype` the item model; the engine settings that are rarely
# changed live in ecm_control() and gibbs_control(). Both build the model object of the old
# constructors (new_model()) and run the same engine, so old and new code give the same fit.

# model names -> internal model types (NA: not available yet)
model_types <- c(irt = "mlirt", null = "rtirt_null", latreg = "rtirt_latreg",
                 latent = "rtirt_latent", cross = "rtirt_cross", mixture = NA)

model_type <- function(model) {
  choices <- paste0("\"", names(model_types), "\"", collapse = ", ")
  if (is.null(model)) {
    stop(sprintf("name the model: model = one of %s", choices), call. = FALSE)
  }
  is_name <- is.character(model) && length(model) == 1 && model %in% names(model_types)
  if (!is_name) {
    stop(sprintf("model must be one of %s", choices), call. = FALSE)
  }
  if (model == "mixture") {
    stop("model = \"mixture\" (Wang & Xu, 2015) is not available yet", call. = FALSE)
  }
  model_types[[model]]
}

check_itemtype <- function(itemtype, n_item) {
  allowed <- c("2pl", "1pl", "grm")
  right_length <- length(itemtype) %in% c(1, n_item)
  if (!is.character(itemtype) || !right_length) {
    stop(sprintf("itemtype must be one value, or one per item (%d)", n_item), call. = FALSE)
  }
  unknown <- setdiff(itemtype, allowed)
  if (length(unknown)) {
    stop(sprintf("unknown itemtype: %s (use \"2pl\", \"1pl\" or \"grm\")", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  if ("grm" %in% itemtype) {
    stop("itemtype = \"grm\" is not available yet", call. = FALSE)
  }
  types <- unique(itemtype)
  if (length(types) > 1) {
    stop("different item types per item are not available yet: give one itemtype for all items", call. = FALSE)
  }
  types
}

# the covariates the model uses: all of `data` (cov = NULL), the named ones, or none
select_cov <- function(data, type, cov) {
  uses_cov <- .model_info[[type]]$X
  if (!uses_cov) {
    if (!is.null(cov)) {
      warning("this model has no covariates; cov is ignored", call. = FALSE)
    }
    data$X <- NULL
    return(data)
  }
  if (!is.null(cov)) {
    missing_cov <- setdiff(cov, colnames(data$X))
    if (length(missing_cov)) {
      stop(sprintf("cov not in the data: %s", paste(missing_cov, collapse = ", ")), call. = FALSE)
    }
    data$X <- data$X[, cov, drop = FALSE]
  }
  no_cov <- is.null(data$X) || ncol(data$X) == 0
  if (type == "mlirt" && no_cov) {
    data$X <- matrix(0, nrow(data$Y), 0)            # "irt" without covariates: plain 2PL
  }
  data
}

build_model <- function(data, model, itemtype, speed_var, cov, quantile = NULL, intercept = FALSE) {
  type <- model_type(model)
  itemtype <- check_itemtype(itemtype, ncol(data$Y))
  data <- select_cov(data, type, cov)
  new_model(type, data, quantile = quantile, speed_var = speed_var, itemtype = itemtype, intercept = intercept)
}

#' Settings of the engines
#'
#' The settings of [ecm()] and [gibbs()] that are rarely changed.
#'
#' @param max_iter Maximum number of ECM (SQUAREM) iterations.
#' @param tol Convergence: relative change of the log-likelihood below `tol` and largest
#'   parameter change (unconstrained scale) below `sqrt(tol)`.
#' @param nodes Gauss-Hermite nodes for ability.
#' @param accelerate Use SQUAREM.
#' @param verbose `ecm_control()`: print the log-likelihood every 10 iterations;
#'   `gibbs_control()`: print a line per chain.
#' @param n_burnin Burn-in per chain (default `n_iter / 2`).
#' @param n_thin Thinning of the stored draws.
#' @param collapse `FALSE` (default): the Gibbs sampler of ExtendedRtIrtModeling.jl, plus the
#'   two group moves. `TRUE`: a partially collapsed Gibbs sampler (van Dyk & Park, 2008).
#'   Given the Polya-Gamma (and ALD) variables the model is linear-Gaussian, so speed is
#'   integrated out in closed form when ability, the RT item parameters (lambda, and rho of
#'   `"cross"`) and, in `"latent"`, the structural regression are drawn, and speed is drawn
#'   after them; in `"irt"` the regression coefficients are drawn with ability integrated
#'   out. The 2PL item parameters (a, b) are drawn jointly with the Polya-Gamma variables
#'   integrated out (a Metropolis-Hastings step with a Fisher-scoring proposal; the
#'   acceptance rate is in `post$ab_accept`). Same posterior, fewer slow directions.
#' @param fit_indices Compute the marginal log-likelihood after sampling (for DIC, WAIC and
#'   PSIS-LOO in [fit_indices()]; see [marginal_loglik()]). `FALSE` skips it (it can take
#'   seconds to a minute for large `"cross"` quantile models); [fit_indices()] then computes
#'   it on request.
#' @return A list of settings.
#' @examples
#' fit_settings <- gibbs_control(n_burnin = 500, collapse = TRUE)
#' fit_settings
#' @export
ecm_control <- function(max_iter = 1000, tol = 1e-10, nodes = 41, accelerate = TRUE, verbose = FALSE) {
  list(max_iter = max_iter, tol = tol, nodes = nodes, accelerate = accelerate, verbose = verbose)
}

#' @rdname ecm_control
#' @export
gibbs_control <- function(n_burnin = NULL, n_thin = 1, collapse = FALSE, fit_indices = TRUE, verbose = TRUE) {
  list(n_burnin = n_burnin, n_thin = n_thin, collapse = collapse, fit_indices = fit_indices, verbose = verbose)
}

#' @rdname ecm
#' @export
ecm.input_data <- function(data, model = NULL, itemtype = "2pl", speed_var = c("free", "fixed"), cov = NULL,
                           prior = c("ml", "map"), priors = NULL, se = TRUE, control = ecm_control(), ...) {
  speed_var <- match.arg(speed_var)
  object <- build_model(data, model, itemtype, speed_var, cov)
  ecm.rtirt(object, max_iter = control$max_iter, tol = control$tol, nodes = control$nodes, se = se,
            prior = prior, priors = priors, accelerate = control$accelerate, verbose = control$verbose, ...)
}

# an ecm() fit as the start of gibbs(): it must be the same data and model
attach_ecm <- function(object, fit) {
  same_model <- identical(class(fit)[1], class(object)[1])
  same_resp <- identical(fit$data$Y, object$data$Y)
  same_times <- identical(fit$data$log_t, object$data$log_t)
  same_cov <- identical(fit$data$X, object$data$X)
  same_data <- same_resp && same_times && same_cov
  key <- c("speed_var", "itemtype")
  same_settings <- identical(fit$settings[key], object$settings[key])
  if (!same_model || !same_data || !same_settings) {
    stop("init must be an ecm() fit of the same data, model, speed_var and itemtype", call. = FALSE)
  }
  object$ecm <- fit$ecm
  object
}

#' @rdname gibbs
#' @export
gibbs.input_data <- function(data, model = NULL, itemtype = "2pl", speed_var = c("free", "fixed"), quantile = NULL,
                             cov = NULL, intercept = FALSE, priors = NULL, n_iter = 5000, n_chain = 4, seed = NULL,
                             init = "random", control = gibbs_control(), ...) {
  speed_var <- match.arg(speed_var)
  object <- build_model(data, model, itemtype, speed_var, cov, quantile, intercept)
  if (inherits(init, "rtirt")) {
    if (!has_em(init)) {
      stop("init must be an ecm() fit, \"random\" or \"ecm\"", call. = FALSE)
    }
    if (inherits(object, "rtirt_qset")) {
      stop("quantile models start from random values", call. = FALSE)
    }
    object <- attach_ecm(object, init)
    init <- "ecm"
  }
  settings <- list(n_iter = n_iter, n_chain = n_chain, n_burnin = control$n_burnin, n_thin = control$n_thin,
                   verbose = control$verbose, fit_indices = control$fit_indices, init = init,
                   priors = priors, collapse = control$collapse)
  args <- c(list(object, seed = seed), settings, list(...))
  sampler <- if (inherits(object, "rtirt_qset")) gibbs.rtirt_qset else gibbs_rtirt
  do.call(sampler, args)
}

#' @rdname equations
#' @param model,itemtype,speed_var,quantile,cov,intercept For `input_data`: the model as in
#'   [ecm()] and [gibbs()], so its equations can be read before estimation, e.g.
#'   `equations(dat, model = "cross", quantile = 0.25)`.
#' @export
equations.input_data <- function(object, model = NULL, itemtype = "2pl", speed_var = c("free", "fixed"), quantile = NULL,
                                 cov = NULL, intercept = FALSE, format = c("unicode", "latex"), priors = TRUE, ...) {
  speed_var <- match.arg(speed_var)
  unfitted <- build_model(object, model, itemtype, speed_var, cov, quantile, intercept)
  equations(unfitted, format = format, priors = priors)
}

#' @rdname algorithm
#' @param model,itemtype,speed_var,quantile,cov,intercept For `input_data`: the model as in
#'   [ecm()] and [gibbs()], e.g. `algorithm(dat, model = "cross", method = "ecm")`.
#' @export
algorithm.input_data <- function(object, model = NULL, itemtype = "2pl", speed_var = c("free", "fixed"), quantile = NULL,
                                 cov = NULL, intercept = FALSE, method = NULL, collapse = NULL,
                                 format = c("unicode", "latex"), ...) {
  speed_var <- match.arg(speed_var)
  unfitted <- build_model(object, model, itemtype, speed_var, cov, quantile, intercept)
  algorithm(unfitted, method = method, collapse = collapse, format = format)
}
