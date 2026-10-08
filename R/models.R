# Model constructors and the sampler. Each constructor returns an S3 object with the fields of
# the Julia structs (cond, data, true_para, para, post); its class is c(<model>, "rtirt"), and
# the mean and quantile (ALD) versions share a class and differ in `qr` / `q`.

.model_info <- list(
  mlirt        = list(engine = "mlirt",  rt = FALSE, X = TRUE),
  rtirt_null   = list(engine = "rtirt",  rt = TRUE,  X = FALSE),
  rtirt_latreg = list(engine = "rtirt",  rt = TRUE,  X = TRUE),
  rtirt_latent = list(engine = "latent", rt = TRUE,  X = TRUE),
  rtirt_cross  = list(engine = "cross",  rt = TRUE,  X = FALSE))

# conditions for a data set: dimensions from the data, MCMC settings from `cond` (if given)
data_cond <- function(cond, data, type) {
  if (is.null(cond)) cond <- set_cond()
  if (!inherits(cond, "sim_cond")) stop("cond must come from set_cond()")
  dims <- c(n_subj = nrow(data$Y), n_item = ncol(data$Y), n_feat = if (.model_info[[type]]$X) ncol(data$X) else 0L)
  for (n in names(dims)) {
    if (is.null(cond[[n]])) cond[[n]] <- unname(dims[n])
    else if (n != "n_feat" || .model_info[[type]]$X) if (cond[[n]] != dims[n])
      stop(sprintf("data has %d persons x %d items%s; cond says %s x %s%s", dims[1], dims[2],
                   if (.model_info[[type]]$X) sprintf(" x %d covariates", dims[3]) else "",
                   cond$n_subj, cond$n_item, if (.model_info[[type]]$X) sprintf(" x %s", cond$n_feat) else ""), call. = FALSE)
  }
  cond
}

new_model <- function(type, data, cond = NULL, true_para = NULL, para = NULL, quantile = NULL,
                      speed_var = "free", itemtype = "2pl", intercept = FALSE) {
  info <- .model_info[[type]]
  if (inherits(data, "sim_cond"))
    stop(sprintf("the data come first: %s(data, cond = cond) (the Julia alias %s keeps the order (Cond, Data))", type,
                 c(mlirt = "GibbsMlIrt()", rtirt_null = "GibbsRtIrtNull()", rtirt_latreg = "GibbsRtIrt()",
                   rtirt_latent = "GibbsRtIrtLatent()", rtirt_cross = "GibbsRtIrtCross()")[[type]]), call. = FALSE)
  if (!inherits(data, "input_data")) stop("data must come from input_data()", call. = FALSE)
  if (info$rt && is.null(data$log_t)) stop(type, " needs response times: input_data(..., time = )", call. = FALSE)
  if (info$X && is.null(data$X)) stop(type, " needs covariates: input_data(..., cov = )", call. = FALSE)
  if (!info$X && !is.null(data$X))
    warning(sprintf("%s does not use covariates; the %d columns of cov are ignored%s", type, ncol(data$X),
                    if (info$rt) " (covariate effects on speed: rtirt_latent() or rtirt_latreg())" else ""), call. = FALSE)
  if (length(quantile) > 1) {                       # several quantile levels: one model each
    ms <- lapply(quantile, function(q) new_model(type, data, cond, true_para, para, q, speed_var, itemtype, intercept))
    return(structure(list(models = ms, quantile = quantile, type = type), class = "rtirt_qset"))
  }
  if (!is.null(quantile) && !(is.numeric(quantile) && quantile > 0 && quantile < 1))
    stop("quantile must be NULL (mean model) or a level in (0, 1)", call. = FALSE)
  cond <- data_cond(cond, data, type)
  if (is.null(true_para)) true_para <- attr(data, "true_para")
  structure(list(cond = cond, data = data, true_para = true_para, para = para, qr = !is.null(quantile),
                 q = if (is.null(quantile)) NA_real_ else quantile,
                 settings = list(speed_var = speed_var, itemtype = itemtype, intercept = intercept), post = NULL),
            class = c(type, "rtirt"))
}

#' Model constructors
#'
#' Create a model object from [input_data()], then sample it with [gibbs()]:
#'
#' ```
#' fit <- gibbs(rtirt_cross(dat, quantile = 0.25), n_iter = 10000, n_chain = 3, seed = 1)
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
mlirt <- function(data, itemtype = c("2pl", "1pl"), intercept = FALSE, cond = NULL, true_para = NULL, para = NULL)
  new_model("mlirt", data, cond, true_para, para, itemtype = match.arg(itemtype), intercept = intercept)
#' @rdname models
#' @export
rtirt_null <- function(data, speed_var = c("free", "fixed"), itemtype = c("2pl", "1pl"), cond = NULL, true_para = NULL, para = NULL)
  new_model("rtirt_null", data, cond, true_para, para, speed_var = match.arg(speed_var), itemtype = match.arg(itemtype))
#' @rdname models
#' @export
rtirt_latreg <- function(data, speed_var = c("free", "fixed"), itemtype = c("2pl", "1pl"), intercept = FALSE,
                         cond = NULL, true_para = NULL, para = NULL)
  new_model("rtirt_latreg", data, cond, true_para, para, speed_var = match.arg(speed_var), itemtype = match.arg(itemtype),
            intercept = intercept)
#' @rdname models
#' @export
rtirt_latent <- function(data, quantile = NULL, itemtype = c("2pl", "1pl"), intercept = FALSE, cond = NULL,
                         true_para = NULL, para = NULL)
  new_model("rtirt_latent", data, cond, true_para, para, quantile, itemtype = match.arg(itemtype), intercept = intercept)
#' @rdname models
#' @export
rtirt_cross <- function(data, quantile = NULL, speed_var = c("free", "fixed"), itemtype = c("2pl", "1pl"),
                        cond = NULL, true_para = NULL, para = NULL)
  new_model("rtirt_cross", data, cond, true_para, para, quantile, speed_var = match.arg(speed_var), itemtype = match.arg(itemtype))

#' Run the Gibbs sampler
#'
#' The R counterpart of `sample!(MCMC)`: runs `n_chain` independent chains (each with its
#' own random starting values) and returns the model with `post` filled in.
#'
#' @param model A model object (see [models]), or an `rtirt_qset` (several quantile levels).
#' @param n_iter,n_chain,n_burnin,n_thin Iterations per chain, chains, burn-in (default
#'   `n_iter / 2`) and thinning of the stored draws. Defaults: the `cond` of the model
#'   (5000 iterations, 4 chains when no `cond` was given).
#' @param seed Random seed (`NULL`: the current RNG state).
#' @param intercept,itemtype Override the settings of the constructor.
#' @param verbose Print a line per chain.
#' @param fit_indices Compute the marginal log-likelihood after sampling (for DIC, WAIC and
#'   PSIS-LOO in [fit_indices()]; see [marginal_loglik()]). `FALSE` skips it (it can take
#'   seconds to a minute for large `rtirt_cross` quantile models); [fit_indices()] then
#'   computes it on request.
#' @param init `"random"` (default): each chain draws its own starting values; `"ecm"`:
#'   start the chains near the maximum likelihood estimates and EAP scores of [ecm()] (run
#'   first when the model has no ECM fit with the same settings), which shortens the burn-in.
#' @param priors Prior distributions, from [rtirt_priors()] (default: `rtirt_priors()`).
#' @param collapse `FALSE` (default): the Gibbs sampler of ExtendedRtIrtModeling.jl, plus the
#'   two group moves. `TRUE`: a partially collapsed Gibbs sampler (van Dyk & Park, 2008).
#'   Given the Polya-Gamma (and ALD) variables the model is linear-Gaussian, so speed is
#'   integrated out in closed form when ability, the RT item parameters (lambda, and rho of
#'   `rtirt_cross()`) and, in `rtirt_latent()`, the structural regression are drawn, and speed
#'   is drawn after them; in `mlirt()` the regression coefficients are drawn with ability
#'   integrated out. The 2PL item parameters (a, b) are drawn jointly with the Polya-Gamma
#'   variables integrated out (a Metropolis-Hastings step with a Fisher-scoring proposal;
#'   the acceptance rate is in `post$ab_accept`). Same posterior, fewer slow directions.
#' @param ... Unused.
#' @return The model with `post`: `draws` (a `coda::mcmc.list` of the item and structural
#'   parameters, all iterations), `mean` (posterior means after burn-in, an [input_para()]),
#'   `person` (posterior mean and SD of ability and speed), `loglik` (complete-data
#'   log-likelihood per iteration and chain) and `secs`. Accessors: `summary()`, `coef()`,
#'   [estimates()], [scores()], [reliability()], [convergence()], [fit_indices()], `plot()`.
#' @examples
#' \donttest{
#' cond <- set_cond(n_subj = 300, n_item = 8)
#' dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' fit <- gibbs(rtirt_cross(dat), n_iter = 1000, n_chain = 2, seed = 1)
#' summary(fit)
#' head(scores(fit))
#' }
#' @export
gibbs <- function(model, ...) UseMethod("gibbs")

#' @rdname gibbs
#' @export
gibbs.rtirt <- function(model, n_iter = NULL, n_chain = NULL, n_burnin = NULL, n_thin = NULL, seed = NULL,
                        intercept = NULL, itemtype = NULL, verbose = TRUE, fit_indices = TRUE,
                        init = c("random", "ecm"), priors = NULL, collapse = FALSE, ...) {
  init <- match.arg(init)
  if (!is.logical(collapse) || length(collapse) != 1 || is.na(collapse)) stop("collapse must be TRUE or FALSE", call. = FALSE)
  priors <- priors %||% rtirt_priors()
  if (!inherits(priors, "rtirt_priors")) stop("priors must come from rtirt_priors()", call. = FALSE)
  pv <- prior_vector(priors)
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  type <- class(model)[1]; info <- .model_info[[type]]; S <- model$settings
  if (!is.null(intercept)) S$intercept <- intercept
  if (!is.null(itemtype)) S$itemtype <- match.arg(itemtype, c("2pl", "1pl"))
  C <- model$cond
  if (!is.null(n_iter)) { C$n_iter <- n_iter; if (is.null(n_burnin)) C$n_burnin <- round(n_iter / 2) }
  if (!is.null(n_burnin)) C$n_burnin <- n_burnin
  if (!is.null(n_chain)) C$n_chain <- n_chain
  if (!is.null(n_thin)) C$n_thin <- n_thin
  if (C$n_burnin >= C$n_iter) stop("n_burnin must be below n_iter")
  if ((C$n_iter - C$n_burnin) / C$n_thin < 2) stop("n_iter - n_burnin must leave at least 2 (thinned) draws per chain", call. = FALSE)
  model$cond <- C; model$settings <- S
  D <- model$data; N <- C$n_subj; K <- C$n_item
  fixed <- S$speed_var == "fixed" && info$engine %in% c("rtirt", "cross")
  X <- if (info$X) D$X else matrix(0, N, 0)
  if (info$X && S$intercept) X <- cbind(`(Intercept)` = 1, X)
  q <- if (model$qr) model$q else 0.5
  one_pl <- S$itemtype == "1pl"
  ab_joint <- collapse && !one_pl                    # (a, b) jointly with omega integrated out
  if (init == "ecm") {
    if (model$qr) stop("init = \"ecm\" needs a normal model; quantile models start from random values", call. = FALSE)
    key <- c("speed_var", "itemtype")                 # an intercept column starts at 0
    if (is.null(model$ecm) || !identical(model$ecm$settings[key], S[key])) {
      m0 <- model; m0$post <- NULL; m0$settings$intercept <- FALSE; model$ecm <- ecm(m0, se = FALSE)$ecm
    }
  }
  chains <- vector("list", C$n_chain); t0 <- Sys.time()
  for (ch in seq_len(C$n_chain)) {
    ini <- if (init == "ecm") em_inits(model, info$engine, X, fixed) else initial_values(info$engine, D, X, fixed)
    t1 <- Sys.time()
    r <- switch(info$engine,
      mlirt  = .gibbs_mlirt(D$Y, X, one_pl, collapse, ab_joint, C$n_iter, C$n_burnin, ini, pv),
      rtirt  = .gibbs_rtirt(D$Y, D$log_t, X, fixed, one_pl, collapse, ab_joint, C$n_iter, C$n_burnin, ini, pv),
      latent = .gibbs_latent(D$Y, D$log_t, X, model$qr, q, one_pl, collapse, ab_joint, C$n_iter, C$n_burnin, ini, pv),
      cross  = .gibbs_cross(D$Y, D$log_t, model$qr, q, fixed, one_pl, collapse, ab_joint, C$n_iter, C$n_burnin, ini, pv))
    colnames(r$draws) <- par_names(info$engine, D$items, colnames(X))
    chains[[ch]] <- r
    if (verbose) message(sprintf("%s: chain %d of %d done (%.1f s)", model_label(model), ch, C$n_chain,
                                 as.numeric(Sys.time() - t1, units = "secs")))
  }
  model$post <- make_post(chains, info, model$qr, C, D$items)
  model$post$secs <- as.numeric(Sys.time() - t0, units = "secs")
  if (ab_joint && !one_pl) model$post$ab_accept <- mean(vapply(chains, function(r) r$ab_accept, 0))
  model$post$settings <- list(intercept = S$intercept, itemtype = S$itemtype, cov2one = fixed, speed_var = S$speed_var,
                              collapse = collapse, priors = priors, q = if (model$qr) q else NA_real_, X = X)
  model$post$estimates <- NULL; model$post$estimates <- estimates(model)
  if (fit_indices) model$post$fit <- marginal_fit(model)
  model
}

#' @rdname gibbs
#' @export
gibbs.rtirt_qset <- function(model, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  model$models <- lapply(model$models, gibbs, ...)
  model
}

# residual variances to start from: a time column without variation has variance 0, and 1 / 0
# breaks the samplers, so floor at 10% of the median (and at 1e-4)
var_floor <- function(v) { v <- unname(v); pmax(v, 0.1 * stats::median(v), 1e-4) }

# random starting values for one chain
initial_values <- function(engine, D, X, fixed) {
  N <- nrow(D$Y); K <- ncol(D$Y); P <- ncol(X)
  init <- list(theta = stats::rnorm(N), a = rep(1, K), b = stats::rnorm(K, 0, 0.1))
  if (engine != "mlirt") init <- c(init, list(zeta = 0.3 * stats::rnorm(N), lambda = unname(colMeans(D$log_t)),
                                               sigma2 = var_floor(apply(D$log_t, 2, stats::var))))
  switch(engine,
    mlirt  = c(init, list(beta = rep(0, P))),
    rtirt  = c(init, list(beta1 = rep(0, P), beta2 = rep(0, P), c = 0, v = if (fixed) 1 else 0.1)),
    latent = c(init, list(beta = rep(0, P + 1), s = 0.1)),
    cross  = c(init, list(rho = stats::rnorm(K, 0, 0.1), s = 1)))
}

par_names <- function(engine, items, xn) {
  it <- function(p) sprintf("%s[%s]", p, items)
  base <- c(it("a"), it("b"), if (engine != "mlirt") c(it("lambda"), it("sigma2t")))
  switch(engine,
    mlirt  = c(base, sprintf("beta[%s]", xn)),
    rtirt  = c(base, sprintf("beta_ability[%s]", xn), sprintf("beta_speed[%s]", xn), "cor_ability_speed", "var_speed"),
    latent = c(base, sprintf("beta[%s]", xn), "b_ability", "var_speed"),
    cross  = c(base, it("rho"), "var_speed"))
}

make_post <- function(chains, info, qr, C, items) {
  keep <- seq(1, C$n_iter, by = C$n_thin)
  draws <- coda::mcmc.list(lapply(chains, function(r) coda::mcmc(r$draws[keep, , drop = FALSE], start = 1, thin = C$n_thin)))
  post <- do.call(rbind, lapply(chains, function(r) r$draws[(C$n_burnin + 1):C$n_iter, , drop = FALSE]))
  pm <- colMeans(post)
  get <- function(p) { i <- startsWith(names(pm), paste0(p, "[")); v <- pm[i]
    names(v) <- sub("^[^[]*\\[(.*)\\]$", "\\1", names(v)); v }
  pool <- function(m, s) {                      # pooled posterior mean and SD over chains
    M <- sapply(chains, `[[`, m); S <- sapply(chains, `[[`, s)
    list(mean = rowMeans(M), sd = sqrt(pmax(rowMeans(S^2 + M^2) - rowMeans(M)^2, 0)))
  }
  pn <- c("ability", "speed")
  person <- list(theta = pool("theta_mean", "theta_sd"))
  if (info$rt) person$zeta <- pool("zeta_mean", "zeta_sd")
  if (!is.null(chains[[1]]$lt_mean))                   # RT person fit, averaged over the chains
    person$rt_fit <- list(lt = rowMeans(sapply(chains, `[[`, "lt_mean")), exceed = rowMeans(sapply(chains, `[[`, "lt_exceed")))
  mean <- input_para(theta = person$theta$mean, a = get("a"), b = get("b"))
  if (info$rt) mean <- input_para(theta = person$theta$mean, a = get("a"), b = get("b"), zeta = person$zeta$mean,
                                  lambda = get("lambda"), sigma2t = get("sigma2t"))
  switch(info$engine,
    mlirt = { mean$beta <- get("beta") },
    rtirt = { mean$beta <- cbind(theta = get("beta_ability"), speed = get("beta_speed"))
              r <- pm[["cor_ability_speed"]]; v <- pm[["var_speed"]]
              mean$sigma_p <- matrix(c(1, r * sqrt(v), r * sqrt(v), v), 2, dimnames = list(pn, pn)) },
    latent = { mean$beta <- get("beta"); mean$b_theta <- pm[["b_ability"]]
               mean$sigma_p <- diag(c(1, pm[["var_speed"]])); dimnames(mean$sigma_p) <- list(pn, pn) },
    cross = { mean$rho <- get("rho"); mean$sigma_p <- diag(c(1, pm[["var_speed"]])); dimnames(mean$sigma_p) <- list(pn, pn) })
  if (qr) mean$nu <- Reduce(`+`, lapply(chains, `[[`, "nu_mean")) / length(chains)
  list(draws = draws, mean = mean, person = person,
       loglik = sapply(chains, `[[`, "loglik"), burnin = C$n_burnin,
       accept = if (!is.null(chains[[1]]$accept)) mean(sapply(chains, `[[`, "accept")) else NA)
}
