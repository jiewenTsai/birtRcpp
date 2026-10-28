# Model constructors and the sampler. Each constructor returns an S3 object with the fields of
# the Julia structs (cond, data, true_para, para, post); its class is c(<model>, "rtirt"), and
# the mean and quantile (ALD) versions share a class and differ in `qr` / `q`.

.model_info <- list(
  mlirt        = list(engine = "mlirt",  rt = FALSE, X = TRUE,  qr = FALSE),
  rtirt_null   = list(engine = "rtirt",  rt = TRUE,  X = FALSE, qr = FALSE),
  rtirt_latreg = list(engine = "rtirt",  rt = TRUE,  X = TRUE,  qr = FALSE),
  rtirt_latent = list(engine = "latent", rt = TRUE,  X = TRUE,  qr = TRUE),
  rtirt_cross  = list(engine = "cross",  rt = TRUE,  X = FALSE, qr = TRUE))

# conditions for a data set: dimensions from the data, MCMC settings from `cond` (if given)
data_cond <- function(cond, data, type) {
  cond <- cond %||% set_cond()
  if (!inherits(cond, "sim_cond")) {
    stop("cond must come from set_cond()")
  }
  uses_cov <- .model_info[[type]]$X
  n_feat <- if (uses_cov) ncol(data$X) else 0L
  dims <- c(n_subj = nrow(data$Y), n_item = ncol(data$Y), n_feat = n_feat)
  for (name in names(dims)) {
    if (is.null(cond[[name]])) {
      cond[[name]] <- unname(dims[name])
      next
    }
    checked <- name != "n_feat" || uses_cov
    if (checked && cond[[name]] != dims[name]) {
      stop_dim_mismatch(dims, cond, uses_cov)
    }
  }
  cond
}

stop_dim_mismatch <- function(dims, cond, uses_cov) {
  data_cov <- if (uses_cov) sprintf(" x %d covariates", dims[3]) else ""
  cond_cov <- if (uses_cov) sprintf(" x %s", cond$n_feat) else ""
  stop(sprintf("data has %d persons x %d items%s; cond says %s x %s%s", dims[1], dims[2], data_cov,
               cond$n_subj, cond$n_item, cond_cov), call. = FALSE)
}

julia_constructor_names <- c(mlirt = "GibbsMlIrt()", rtirt_null = "GibbsRtIrtNull()", rtirt_latreg = "GibbsRtIrt()",
                             rtirt_latent = "GibbsRtIrtLatent()", rtirt_cross = "GibbsRtIrtCross()")

new_model <- function(type, data, cond = NULL, true_para = NULL, para = NULL, quantile = NULL,
                      speed_var = "free", itemtype = "2pl", intercept = FALSE) {
  check_model_data(type, data)
  if (length(quantile) > 1) {
    return(new_qset(type, data, cond, true_para, para, quantile, speed_var, itemtype, intercept))
  }
  if (!is.null(quantile) && !(is.numeric(quantile) && quantile > 0 && quantile < 1)) {
    stop("quantile must be NULL (mean model) or a level in (0, 1)", call. = FALSE)
  }
  cond <- data_cond(cond, data, type)
  true_para <- true_para %||% attr(data, "true_para")
  settings <- list(speed_var = speed_var, itemtype = itemtype, intercept = intercept)
  model <- list(cond = cond, data = data, true_para = true_para, para = para, qr = !is.null(quantile),
                q = if (is.null(quantile)) NA_real_ else quantile, settings = settings, post = NULL)
  structure(model, class = c(type, "rtirt"))
}

check_model_data <- function(type, data) {
  info <- .model_info[[type]]
  if (inherits(data, "sim_cond")) {
    stop(sprintf("the data come first: %s(data, cond = cond) (the Julia alias %s keeps the order (Cond, Data))", type,
                 julia_constructor_names[[type]]), call. = FALSE)
  }
  if (!inherits(data, "input_data")) {
    stop("data must come from input_data()", call. = FALSE)
  }
  if (info$rt && is.null(data$log_t)) {
    stop(type, " needs response times: input_data(..., time = )", call. = FALSE)
  }
  if (info$X && is.null(data$X)) {
    stop(type, " needs covariates: input_data(..., cov = )", call. = FALSE)
  }
  if (!info$X && !is.null(data$X)) {
    hint <- if (info$rt) " (covariate effects on speed: model = \"latent\" or \"latreg\")" else ""
    warning(sprintf("%s does not use covariates; the %d columns of cov are ignored%s", type, ncol(data$X), hint),
            call. = FALSE)
  }
}

# several quantile levels: one model each
new_qset <- function(type, data, cond, true_para, para, quantile, speed_var, itemtype, intercept) {
  one_level <- function(q) new_model(type, data, cond, true_para, para, q, speed_var, itemtype, intercept)
  models <- lapply(quantile, one_level)
  structure(list(models = models, quantile = quantile, type = type), class = "rtirt_qset")
}

#' Model constructors (deprecated)
#'
#' **Deprecated**: name the model in [ecm()] or [gibbs()] instead, e.g.
#' `rtirt_cross(dat)` becomes `ecm(dat, model = "cross")` or `gibbs(dat, model = "cross")`.
#' The constructors still work (with a warning once per session) and get no new features.
#'
#' ```
#' fit <- gibbs(dat, model = "cross", quantile = 0.25, n_iter = 10000, n_chain = 3, seed = 1)
#' ```
#'
#' * `mlirt`: 2PL with a latent regression of ability (no response times).
#' * `rtirt_null`: joint RT-IRT (van der Linden, 2007), ability and speed correlated.
#' * `rtirt_latreg`: joint RT-IRT with latent regressions of ability and speed on covariates
#'   (dissertation Chapter 3).
#' * `rtirt_latent`: speed regressed on ability and covariates, mean regression or, with
#'   `quantile`, quantile regression (ALD; Chapter 4).
#' * `rtirt_cross`: cross-relations of ability with each log response time, normal or, with
#'   `quantile`, ALD residuals (Chapter 5).
#'
#' The Julia constructors (`GibbsRtIrtCross(Cond, Data)` etc.) are aliases with the Julia
#' argument order and defaults (see [julia_aliases]).
#'
#' @param data Data from [input_data()].
#' @param quantile `NULL` for the mean model, or the quantile level of the ALD model
#'   (`rtirt_latent`: of speed given ability and covariates; `rtirt_cross`: of each log
#'   response time). A vector fits one model per level (an `rtirt_qset`; [gibbs()] samples
#'   all, `plot()` draws the coefficients against the level).
#' @param speed_var `"free"` (default) estimates the variance of speed; `"fixed"` sets it
#'   to 1 (`rtirt_latreg`, `rtirt_null`: ability and speed then have a correlation matrix),
#'   as the Julia package does. For log seconds the speed variance is usually far below 1,
#'   so fixing it over-constrains the model.
#' @param itemtype `"2pl"` or `"1pl"` (discriminations fixed at 1).
#' @param intercept Add an intercept to the latent regressions (default `FALSE`: the
#'   covariates are centred and the item parameters carry the location).
#' @param cond Optional [set_cond()] (simulations): the dimensions are checked against the
#'   data, and its MCMC settings are the defaults of [gibbs()].
#' @param true_para Optional true parameters ([input_para()]) for simulations; by default
#'   those stored by [sim_data()].
#' @param para Optional initial values (not needed; each chain draws its own).
#' @return A model object of class `c(<model>, "rtirt")`; see [gibbs()].
#' @name models
NULL

#' @rdname models
#' @export
mlirt <- function(data, itemtype = c("2pl", "1pl"), intercept = FALSE, cond = NULL, true_para = NULL, para = NULL) {
  deprecate_constructor("mlirt", "irt")
  new_model("mlirt", data, cond, true_para, para, itemtype = match.arg(itemtype), intercept = intercept)
}
#' @rdname models
#' @export
rtirt_null <- function(data, speed_var = c("free", "fixed"), itemtype = c("2pl", "1pl"), cond = NULL, true_para = NULL, para = NULL) {
  deprecate_constructor("rtirt_null", "null")
  new_model("rtirt_null", data, cond, true_para, para, speed_var = match.arg(speed_var), itemtype = match.arg(itemtype))
}
#' @rdname models
#' @export
rtirt_latreg <- function(data, speed_var = c("free", "fixed"), itemtype = c("2pl", "1pl"), intercept = FALSE,
                         cond = NULL, true_para = NULL, para = NULL) {
  deprecate_constructor("rtirt_latreg", "latreg")
  new_model("rtirt_latreg", data, cond, true_para, para, speed_var = match.arg(speed_var), itemtype = match.arg(itemtype),
            intercept = intercept)
}
#' @rdname models
#' @export
rtirt_latent <- function(data, quantile = NULL, itemtype = c("2pl", "1pl"), intercept = FALSE, cond = NULL,
                         true_para = NULL, para = NULL) {
  deprecate_constructor("rtirt_latent", "latent")
  new_model("rtirt_latent", data, cond, true_para, para, quantile, itemtype = match.arg(itemtype), intercept = intercept)
}
#' @rdname models
#' @export
rtirt_cross <- function(data, quantile = NULL, speed_var = c("free", "fixed"), itemtype = c("2pl", "1pl"),
                        cond = NULL, true_para = NULL, para = NULL) {
  deprecate_constructor("rtirt_cross", "cross")
  new_model("rtirt_cross", data, cond, true_para, para, quantile, speed_var = match.arg(speed_var), itemtype = match.arg(itemtype))
}

#' Run the Gibbs sampler
#'
#' Samples a built-in model from its full conditionals: `n_chain` independent chains, each
#' from its own random starting values or from an [ecm()] fit.
#'
#' ```
#' fit_ecm <- ecm(dat, model = "cross")
#' fit_gibbs <- gibbs(dat, model = "cross", init = fit_ecm)
#' ```
#'
#' @param data Data from [input_data()].
#' @param model The model, by name:
#'   * `"irt"`: 2PL (no response times); with covariates, a latent regression of ability.
#'   * `"null"`: joint RT-IRT (van der Linden, 2007), ability and speed correlated.
#'   * `"latreg"`: joint RT-IRT with latent regressions of ability and speed on the covariates.
#'   * `"latent"`: speed regressed on ability and the covariates (mean or, with `quantile`,
#'     quantile regression).
#'   * `"cross"`: cross-relations of ability with each log response time (normal or, with
#'     `quantile`, ALD residuals).
#'   * `"mixture"`: the mixture model of Wang and Xu (2015); not available yet.
#' @param itemtype `"2pl"`, `"1pl"` (discriminations fixed at 1) or `"grm"` (not available
#'   yet); one value, or one per item (all equal for now).
#' @param speed_var `"free"` (default) estimates the variance of speed; `"fixed"` sets it to 1
#'   (`"null"`, `"latreg"`, `"cross"`).
#' @param quantile `NULL` for the mean model, or the quantile level of the ALD model
#'   (`"latent"`: of speed given ability and covariates; `"cross"`: of each log response
#'   time). A vector fits one model per level (an `rtirt_qset`). Only `"latent"` and
#'   `"cross"` have quantile versions; `quantile =` with another model is an error.
#' @param cov Names of the covariates of `data` to use (default all); models without
#'   covariates ignore them.
#' @param intercept Add an intercept to the latent regressions (default `FALSE`: the
#'   covariates are centred and the item parameters carry the location).
#' @param priors Prior distributions, from [rtirt_priors()] (default: `rtirt_priors()`).
#' @param n_iter,n_chain Iterations per chain and chains.
#' @param seed Random seed (`NULL`: the current RNG state).
#' @param init `"random"` (default): each chain draws its own starting values; an [ecm()]
#'   fit of the same data and model: start the chains at its estimates and EAP scores, which
#'   shortens the burn-in; `"ecm"`: run [ecm()] first and start there.
#' @param control Further sampler settings, from [gibbs_control()], e.g.
#'   `gibbs(dat, model = "cross", control = gibbs_control(n_burnin = 500))`.
#' @param ... Not used: any other argument is an error before sampling starts. A setting of
#'   [gibbs_control()] (`n_burnin`, `n_thin`, `collapse`, `fit_indices`, `verbose`) gets the
#'   `control = gibbs_control(...)` form in the message, and a misspelled name (`niter`) the
#'   closest argument (`n_iter`).
#' @return A fitted model with `post`: `draws` (a `coda::mcmc.list` of the item and structural
#'   parameters, all iterations), `mean` (posterior means after burn-in, an [input_para()]),
#'   `person` (posterior mean and SD of ability and speed), `loglik` (complete-data
#'   log-likelihood per iteration and chain) and `secs`. Accessors: `summary()`, `coef()`,
#'   [estimates()], [scores()], [reliability()], [convergence()], [fit_indices()], `plot()`.
#' @references Wang, C., & Xu, G. (2015). A mixture hierarchical model for response times and
#'   response accuracy. *British Journal of Mathematical and Statistical Psychology*, 68,
#'   456-477.
#' @examples
#' \donttest{
#' cond <- set_cond(n_subj = 300, n_item = 8)
#' dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' fit <- gibbs(dat, model = "cross", n_iter = 1000, n_chain = 2, seed = 1)
#' summary(fit)
#' head(scores(fit))
#' }
#' @export
gibbs <- function(data, ...) UseMethod("gibbs")

#' Sample a model object (deprecated)
#'
#' The interface before `gibbs(data, model = ...)`: [gibbs()] on a model object from the
#' deprecated constructors ([models]), or on an [ecm()] fit with `init = "ecm"`. Use
#' `gibbs(data, model = ..., init = fit_ecm)` instead.
#'
#' @param data A model object (see [models]), an [ecm()] fit, or an `rtirt_qset` (several
#'   quantile levels).
#' @param seed Random seed (`NULL`: the current RNG state).
#' @param ... Settings of the sampler: `n_iter`, `n_chain`, `n_burnin` (default `n_iter / 2`) and
#'   `n_thin`, with defaults from the `cond` of the model (5000 iterations, 4 chains when no `cond`
#'   was given); `intercept` and `itemtype`, which override the settings of the constructor;
#'   `verbose`, `fit_indices` and `collapse` (see [gibbs_control()]); `init` (`"random"` or
#'   `"ecm"`); `priors`, from [rtirt_priors()].
#' @return The model with `post`, as [gibbs()].
#' @name gibbs_model
#' @keywords internal
NULL

#' @rdname gibbs_model
#' @export
gibbs.rtirt <- function(data, ...) {
  if (has_em(data)) {
    deprecate_warn("gibbs(fit, init = \"ecm\")", "gibbs(data, model = ..., init = fit)")
  }
  gibbs_rtirt(data, ...)
}

# the sampler of a model object
gibbs_rtirt <- function(model, n_iter = NULL, n_chain = NULL, n_burnin = NULL, n_thin = NULL, seed = NULL,
                        intercept = NULL, itemtype = NULL, verbose = TRUE, fit_indices = TRUE,
                        init = c("random", "ecm"), priors = NULL, collapse = FALSE, ...) {
  init <- match.arg(init)
  if (!is.logical(collapse) || length(collapse) != 1 || is.na(collapse)) {
    stop("collapse must be TRUE or FALSE", call. = FALSE)
  }
  priors <- priors %||% rtirt_priors()
  if (is.character(priors)) {
    # prior = "map" of ecm() reaches here through partial matching of priors =
    stop("priors must come from rtirt_priors(); prior = \"ml\" or \"map\" is a setting of ecm(): ",
         "gibbs() always samples the posterior under the priors of rtirt_priors()", call. = FALSE)
  }
  if (!inherits(priors, "rtirt_priors")) {
    stop("priors must come from rtirt_priors()", call. = FALSE)
  }
  prior_values <- prior_vector(priors)
  if (length(list(...))) {
    warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  type <- class(model)[1]
  info <- .model_info[[type]]
  settings <- model$settings
  settings$intercept <- intercept %||% settings$intercept
  if (!is.null(itemtype)) {
    settings$itemtype <- match.arg(itemtype, c("2pl", "1pl"))
  }
  cond <- mcmc_cond(model$cond, n_iter, n_chain, n_burnin, n_thin)
  model$cond <- cond
  model$settings <- settings

  design <- sampler_design(model, info, settings, collapse)
  if (init == "ecm") {
    model <- ecm_start(model, settings)
  }

  start_time <- Sys.time()
  chains <- lapply(seq_len(cond$n_chain), function(chain) {
    run_gibbs_chain(model, info, design, init, prior_values, chain, verbose)
  })

  model$post <- make_post(chains, info, model$qr, cond, model$data$items)
  model$post$secs <- as.numeric(Sys.time() - start_time, units = "secs")
  if (design$ab_joint && !design$one_pl) {
    model$post$ab_accept <- mean(vapply(chains, function(r) r$ab_accept, 0))
  }
  model$post$settings <- list(intercept = settings$intercept, itemtype = settings$itemtype, cov2one = design$fixed,
                              speed_var = settings$speed_var, collapse = collapse, priors = priors,
                              q = if (model$qr) design$q else NA_real_, X = design$X)
  model$post$estimates <- NULL
  model$post$estimates <- estimates(model)
  if (fit_indices) {
    model$post$fit <- marginal_fit(model)
  }
  model
}

# MCMC settings: the arguments of gibbs() override those of the model's conditions
mcmc_cond <- function(cond, n_iter, n_chain, n_burnin, n_thin) {
  if (!is.null(n_iter)) {
    cond$n_iter <- n_iter
    cond$n_burnin <- n_burnin %||% round(n_iter / 2)
  }
  cond$n_burnin <- n_burnin %||% cond$n_burnin
  cond$n_chain <- n_chain %||% cond$n_chain
  cond$n_thin <- n_thin %||% cond$n_thin
  if (cond$n_burnin >= cond$n_iter) {
    stop("n_burnin must be below n_iter")
  }
  if ((cond$n_iter - cond$n_burnin) / cond$n_thin < 2) {
    stop("n_iter - n_burnin must leave at least 2 (thinned) draws per chain", call. = FALSE)
  }
  cond
}

# what the C++ sampler needs besides the data: design matrix, quantile level and switches
sampler_design <- function(model, info, settings, collapse) {
  n_person <- model$cond$n_subj
  X <- if (info$X) model$data$X else matrix(0, n_person, 0)
  if (info$X && settings$intercept) {
    X <- cbind(`(Intercept)` = 1, X)
  }
  one_pl <- settings$itemtype == "1pl"
  list(X = X,
       fixed = settings$speed_var == "fixed" && info$engine %in% c("rtirt", "cross"),
       q = if (model$qr) model$q else 0.5,
       one_pl = one_pl,
       collapse = collapse,
       ab_joint = collapse && !one_pl)               # (a, b) jointly with omega integrated out
}

# init = "ecm": an ECM fit of the same model (an intercept column starts at 0)
ecm_start <- function(model, settings) {
  if (model$qr) {
    stop("init = \"ecm\" needs a normal model; quantile models start from random values", call. = FALSE)
  }
  key <- c("speed_var", "itemtype")
  if (is.null(model$ecm) || !identical(model$ecm$settings[key], settings[key])) {
    model_ecm <- model
    model_ecm$post <- NULL
    model_ecm$settings$intercept <- FALSE
    model$ecm <- ecm(model_ecm, se = FALSE)$ecm
  }
  model
}

run_gibbs_chain <- function(model, info, design, init, prior_values, chain, verbose) {
  data <- model$data
  cond <- model$cond
  start <- if (init == "ecm") {
    em_inits(model, info$engine, design$X, design$fixed)
  } else {
    initial_values(info$engine, data, design$X, design$fixed)
  }
  chain_start <- Sys.time()
  result <- call_gibbs_cpp(info$engine, data, model$qr, design, cond, start, prior_values)
  colnames(result$draws) <- par_names(info$engine, data$items, colnames(design$X))
  if (verbose) {
    seconds <- as.numeric(Sys.time() - chain_start, units = "secs")
    message(sprintf("%s: chain %d of %d done (%.1f s)", model_label(model), chain, cond$n_chain, seconds))
  }
  result
}

# one chain of the C++ sampler of the engine
call_gibbs_cpp <- function(engine, data, qr, design, cond, start, prior_values) {
  d <- design
  switch(engine,
    mlirt  = .gibbs_mlirt(data$Y, d$X, d$one_pl, d$collapse, d$ab_joint, cond$n_iter, cond$n_burnin, start, prior_values),
    rtirt  = .gibbs_rtirt(data$Y, data$log_t, d$X, d$fixed, d$one_pl, d$collapse, d$ab_joint, cond$n_iter, cond$n_burnin,
                          start, prior_values),
    latent = .gibbs_latent(data$Y, data$log_t, d$X, qr, d$q, d$one_pl, d$collapse, d$ab_joint, cond$n_iter,
                           cond$n_burnin, start, prior_values),
    cross  = .gibbs_cross(data$Y, data$log_t, qr, d$q, d$fixed, d$one_pl, d$collapse, d$ab_joint, cond$n_iter,
                          cond$n_burnin, start, prior_values))
}

#' @rdname gibbs_model
#' @export
gibbs.rtirt_qset <- function(data, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data$models <- lapply(data$models, gibbs_rtirt, ...)
  data
}

# residual variances to start from: a time column without variation has variance 0, and 1 / 0
# breaks the samplers, so floor at 10% of the median (and at 1e-4)
var_floor <- function(v) {
  v <- unname(v)
  pmax(v, 0.1 * stats::median(v), 1e-4)
}

# random starting values for one chain (drawn in this order: theta, b, zeta, rho)
initial_values <- function(engine, D, X, fixed) {
  n_person <- nrow(D$Y)
  n_item <- ncol(D$Y)
  n_cov <- ncol(X)
  theta <- stats::rnorm(n_person)
  b <- stats::rnorm(n_item, 0, 0.1)
  init <- list(theta = theta, a = rep(1, n_item), b = b)
  if (engine != "mlirt") {
    zeta <- 0.3 * stats::rnorm(n_person)
    lambda <- unname(colMeans(D$log_t))
    sigma2 <- var_floor(apply(D$log_t, 2, stats::var))
    init <- c(init, list(zeta = zeta, lambda = lambda, sigma2 = sigma2))
  }
  switch(engine,
    mlirt  = c(init, list(beta = rep(0, n_cov))),
    rtirt  = c(init, list(beta1 = rep(0, n_cov), beta2 = rep(0, n_cov), c = 0, v = if (fixed) 1 else 0.1)),
    latent = c(init, list(beta = rep(0, n_cov + 1), s = 0.1)),
    cross  = c(init, list(rho = stats::rnorm(n_item, 0, 0.1), s = 1)))
}

par_names <- function(engine, items, xn) {
  per_item <- function(p) sprintf("%s[%s]", p, items)
  rt_items <- if (engine != "mlirt") c(per_item("lambda"), per_item("sigma2t"))
  base <- c(per_item("a"), per_item("b"), rt_items)
  switch(engine,
    mlirt  = c(base, sprintf("beta[%s]", xn)),
    rtirt  = c(base, sprintf("beta_ability[%s]", xn), sprintf("beta_speed[%s]", xn), "cor_ability_speed", "var_speed"),
    latent = c(base, sprintf("beta[%s]", xn), "b_ability", "var_speed"),
    cross  = c(base, per_item("rho"), "var_speed"))
}

make_post <- function(chains, info, qr, C, items) {
  keep <- seq(1, C$n_iter, by = C$n_thin)
  as_mcmc <- function(r) coda::mcmc(r$draws[keep, , drop = FALSE], start = 1, thin = C$n_thin)
  draws <- coda::mcmc.list(lapply(chains, as_mcmc))
  after_burnin <- (C$n_burnin + 1):C$n_iter
  pooled_draws <- do.call(rbind, lapply(chains, function(r) r$draws[after_burnin, , drop = FALSE]))
  post_mean <- colMeans(pooled_draws)

  person <- pooled_persons(chains, info)
  mean <- posterior_means(post_mean, person, info)
  if (qr) {
    mean$nu <- Reduce(`+`, lapply(chains, `[[`, "nu_mean")) / length(chains)
  }
  accept <- if (!is.null(chains[[1]]$accept)) mean(sapply(chains, `[[`, "accept")) else NA
  list(draws = draws, mean = mean, person = person,
       loglik = sapply(chains, `[[`, "loglik"), burnin = C$n_burnin, accept = accept)
}

# posterior means and SDs of the person parameters, pooled over the chains
pooled_persons <- function(chains, info) {
  pool <- function(mean_name, sd_name) {
    chain_means <- sapply(chains, `[[`, mean_name)
    chain_sds <- sapply(chains, `[[`, sd_name)
    pooled_mean <- rowMeans(chain_means)
    pooled_var <- rowMeans(chain_sds^2 + chain_means^2) - pooled_mean^2
    list(mean = pooled_mean, sd = sqrt(pmax(pooled_var, 0)))
  }
  person <- list(theta = pool("theta_mean", "theta_sd"))
  if (info$rt) {
    person$zeta <- pool("zeta_mean", "zeta_sd")
  }
  if (!is.null(chains[[1]]$lt_mean)) {                # RT person fit, averaged over the chains
    lt <- rowMeans(sapply(chains, `[[`, "lt_mean"))
    exceed <- rowMeans(sapply(chains, `[[`, "lt_exceed"))
    person$rt_fit <- list(lt = lt, exceed = exceed)
  }
  person
}

# posterior means as an input_para (the Julia names), from the pooled draws
posterior_means <- function(post_mean, person, info) {
  block <- function(p) {
    values <- post_mean[startsWith(names(post_mean), paste0(p, "["))]
    names(values) <- sub("^[^[]*\\[(.*)\\]$", "\\1", names(values))
    values
  }
  person_names <- c("ability", "speed")
  diagonal_cov <- function(var_speed) {
    matrix(c(1, 0, 0, var_speed), 2, dimnames = list(person_names, person_names))
  }

  mean <- input_para(theta = person$theta$mean, a = block("a"), b = block("b"))
  if (info$rt) {
    mean <- input_para(theta = person$theta$mean, a = block("a"), b = block("b"), zeta = person$zeta$mean,
                       lambda = block("lambda"), sigma2t = block("sigma2t"))
  }
  if (info$engine == "mlirt") {
    mean$beta <- block("beta")
  }
  if (info$engine == "rtirt") {
    mean$beta <- cbind(theta = block("beta_ability"), speed = block("beta_speed"))
    cor_ability_speed <- post_mean[["cor_ability_speed"]]
    var_speed <- post_mean[["var_speed"]]
    cov_ability_speed <- cor_ability_speed * sqrt(var_speed)
    mean$sigma_p <- matrix(c(1, cov_ability_speed, cov_ability_speed, var_speed), 2,
                           dimnames = list(person_names, person_names))
  }
  if (info$engine == "latent") {
    mean$beta <- block("beta")
    mean$b_theta <- post_mean[["b_ability"]]
    mean$sigma_p <- diagonal_cov(post_mean[["var_speed"]])
  }
  if (info$engine == "cross") {
    mean$rho <- block("rho")
    mean$sigma_p <- diagonal_cov(post_mean[["var_speed"]])
  }
  mean
}
