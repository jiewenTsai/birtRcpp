# Declared models for the building blocks, in the style of PyMC: the model is declared once
# (distributions and dimensions), and samplers are stacked as steps that name the parameters they
# update. Every step derives its full conditional, gradient or shift from the declared model, so the
# steps cannot disagree about the model; the Polya-Gamma variables of logistic parts are drawn
# automatically wherever a step needs them. Persons j are the rows and items i the columns of the
# data matrices.

MODEL_DISTS <- c("normal", "bernoulli_logit", "half_t", "inv_gamma")

#' Declare a model for the building blocks
#'
#' `block_model()` declares a model once, the way PyMC does: the data, and one line
#' `name ~ distribution(...)` per parameter and per observed matrix. Samplers are then stacked as
#' steps that only name the parameters they update ([step_gibbs()], [step_mala()], [step_shift()],
#' run by [gibbs()]). Each step derives what it needs from the model, so
#'
#' * the steps always target the same posterior (a prior cannot be forgotten in one step);
#' * the Pólya–Gamma variables \eqn{\omega} of a logistic part are hidden: they are drawn
#'   whenever a Gibbs step needs them and are fresh (drawn once per iteration and again after a
#'   MALA step has moved a parameter of that part);
#' * positive parameters are moved on the log scale by MALA, with the Jacobian added;
#' * dimensions are declared once and checked.
#'
#' **Data**: named arguments. Matrices are persons x items (all the same size); vectors have one
#' value per person (length N) or per item (length K); missing responses or log times are `NA`.
#'
#' **Parameters**: `name[person] ~ ...`, `name[item] ~ ...` or `name ~ ...` (a single value).
#'
#' **Distributions**
#' * `normal(mean, sd, lower = -Inf, upper = Inf)`: `mean` and `sd` may use parameters and data,
#'   e.g. `zeta[person] ~ normal(beta * theta, sqrt(tau2))`. For an observed matrix the mean is
#'   written with the person and item parameters as they are; they are expanded to the matrix.
#' * `bernoulli_logit(eta)`: logistic responses, e.g. `Y ~ bernoulli_logit(a * (theta - b))`.
#' * `half_t(df, scale)` and `inv_gamma(shape, scale)`: priors of a **variance** `v` that is used as
#'   `sd = sqrt(v)`; `half_t()` is the prior of the SD \eqn{\sqrt v}{sqrt(v)} (drawn through the
#'   inverse-gamma auxiliary variable of Huang & Wand, 2013).
#'
#' Expressions may use `+ - * / ^`, `exp()`, `log()` and `sqrt()`. A line `name <- expression`
#' defines a name for later lines, as in `rethinking::ulam()`, e.g. `eta <- a * (theta - b)` and then
#' `Y ~ bernoulli_logit(eta)`; definitions are substituted (they are not parameters).
#'
#' @param ... Data (named arguments) and declarations (formulas).
#' @param object,x A declared model.
#' @return An object of class `block_model`; `equations()` prints the model.
#' @seealso [gibbs()] and the steps; [run_sampler()] for samplers written by hand.
#' @examples
#' set.seed(1); N <- 300; K <- 6
#' th <- rnorm(N); a0 <- runif(K, 0.8, 2); b0 <- seq(-1.5, 1.5, length.out = K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(th, a0), 2, a0 * b0))), N)
#' m <- block_model(Y = Y,
#'   theta[person] ~ normal(0, 1),
#'   a[item] ~ normal(1, 1, lower = 0),
#'   b[item] ~ normal(0, 3),
#'   Y ~ bernoulli_logit(a * (theta - b)))
#' m
#' fit <- gibbs(m, list(step_mala(theta), step_gibbs(a, b)), n_iter = 600)
#' summary(fit)[1:4, ]
#' algorithm(fit)
#' @export
block_model <- function(...) {
  env <- parent.frame()
  args <- as.list(substitute(list(...)))[-1]
  nm <- names(args) %||% rep("", length(args))
  data <- list(); decl <- list()
  for (k in seq_along(args)) {
    if (nzchar(nm[k])) data[[nm[k]]] <- eval(args[[k]], env)
    else decl[[length(decl) + 1]] <- args[[k]]
  }
  bm_build(data, decl, env)
}

# the model from evaluated data and unevaluated declarations (bounds and hyperparameters are evaluated in env)
bm_build <- function(data, decl, env) {
  err <- function(...) stop(sprintf(...), call. = FALSE)
  mats <- Filter(is.matrix, data)
  if (!length(mats)) err("block_model(): give the data as named persons x items matrices, e.g. Y = Y")
  N <- nrow(mats[[1]]); K <- ncol(mats[[1]])
  for (n in names(mats)) if (!identical(dim(mats[[n]]), c(N, K)))
    err("block_model(): %s is %d x %d but %s is %d x %d (all matrices are persons x items)", n, nrow(mats[[n]]), ncol(mats[[n]]), names(mats)[1], N, K)
  ddim <- vapply(names(data), function(n) {
    x <- data[[n]]
    if (is.matrix(x)) "NK" else if (length(x) == 1) "scalar"
    else if (length(x) == N && N != K) "person" else if (length(x) == K && N != K) "item"
    else if (length(x) == N) err("block_model(): %s has length %d, and N = K = %d: give it as an N x K matrix", n, N, N)
    else err("block_model(): %s has length %d; vectors need one value per person (%d) or per item (%d)", n, length(x), N, K)
  }, "")
  for (n in names(data)) data[[n]] <- if (is.matrix(data[[n]])) { storage.mode(data[[n]]) <- "double"; data[[n]] } else as.numeric(data[[n]])
  data0 <- lapply(data, function(x) { x[is.na(x)] <- 0; x })           # missing values as 0 (masked where used)

  pars <- list(); facs <- list(); defs <- list()
  for (f in decl) {
    if (is.call(f) && identical(f[[1]], as.name("<-"))) {          # a definition, as in rethinking::ulam(): eta <- a * (theta - b)
      if (!is.name(f[[2]])) err("block_model(): the left side of '%s' must be a name", deparse1(f))
      dn <- as.character(f[[2]])
      if (dn %in% c(names(data), names(pars), names(defs))) err("block_model(): %s is defined twice (or is data or a parameter)", dn)
      defs[[dn]] <- do.call(substitute, list(f[[3]], defs))
      next
    }
    if (!is.call(f) || !identical(f[[1]], as.name("~")) || length(f) != 3)
      err("block_model(): '%s' is not a declaration 'name ~ distribution(...)'", deparse1(f))
    lhs <- f[[2]]; rhs <- do.call(substitute, list(f[[3]], defs)); dim <- "scalar"
    if (is.call(lhs) && identical(lhs[[1]], as.name("["))) { dim <- as.character(lhs[[3]]); lhs <- lhs[[2]] }
    if (!is.name(lhs)) err("block_model(): the left side of '%s' must be a name, e.g. theta[person]", deparse1(f))
    y <- as.character(lhs)
    if (!dim %in% c("scalar", "person", "item")) err("block_model(): '%s': the dimension is [person] or [item] (or none for one value)", deparse1(f))
    d <- if (is.call(rhs)) as.character(rhs[[1]]) else ""
    if (!d %in% MODEL_DISTS) err("block_model(): '%s': the distribution must be one of %s", deparse1(f), paste0(MODEL_DISTS, "()", collapse = ", "))
    a <- switch(d, normal = match.call(function(mean, sd, lower = -Inf, upper = Inf) NULL, rhs),
                bernoulli_logit = match.call(function(eta) NULL, rhs),
                half_t = match.call(function(df, scale) NULL, rhs), inv_gamma = match.call(function(shape, scale) NULL, rhs))
    a <- as.list(a)[-1]
    obs <- y %in% names(data)
    if (obs && dim != "scalar") err("block_model(): %s is data; its dimension comes from the data, write '%s ~ ...'", y, y)
    if (!obs && y %in% names(pars)) err("block_model(): %s is declared twice", y)
    fac <- list(y = y, dist = d, observed = obs, shape = if (obs) ddim[[y]] else dim, decl = f)
    if (d == "normal") {
      if (is.null(a$mean) || is.null(a$sd)) err("block_model(): '%s': normal() needs mean and sd", deparse1(f))
      fac$mean <- a$mean; fac$sd <- a$sd
      fac$lower <- eval(a$lower %||% -Inf, env); fac$upper <- eval(a$upper %||% Inf, env)
      if (obs && (is.finite(fac$lower) || is.finite(fac$upper))) err("block_model(): truncation is for parameters, not for data (%s)", y)
    } else if (d == "bernoulli_logit") {
      if (!obs) err("block_model(): bernoulli_logit() is for observed 0/1 data (%s)", y)
      if (is.null(a$eta)) err("block_model(): '%s': bernoulli_logit() needs the linear predictor", deparse1(f))
      fac$eta <- a$eta
      v <- data[[y]][!is.na(data[[y]])]; if (!all(v %in% c(0, 1))) err("block_model(): %s must be 0/1 (or NA) for bernoulli_logit()", y)
    } else {
      if (obs) err("block_model(): %s() is a prior of a variance parameter, not of data (%s)", d, y)
      p <- vapply(a, function(z) as.numeric(eval(z, env)), 0)
      if (length(p) != 2 || any(!(p > 0))) err("block_model(): '%s': %s() needs two positive numbers", deparse1(f), d)
      fac$hyper <- p
    }
    if (!obs) pars[[y]] <- list(dim = dim, dist = d, lower = fac$lower %||% (if (d %in% c("half_t", "inv_gamma")) 0 else -Inf),
                                upper = fac$upper %||% Inf, hyper = fac$hyper, mean = fac$mean, sd = fac$sd)
    facs[[length(facs) + 1]] <- fac
  }
  used <- unique(vapply(Filter(function(z) z$observed, facs), `[[`, "", "y"))
  M <- structure(list(data = data, data0 = data0, ddim = ddim, pars = pars, factors = NULL, N = N, K = K, used = used, defs = defs), class = "block_model")
  pn <- names(pars)
  for (k in seq_along(facs)) {
    f <- facs[[k]]; f$id <- k
    ex <- Filter(Negate(is.null), list(f$mean, f$sd, f$eta))
    sy <- unique(unlist(lapply(ex, all.vars)))
    bad <- setdiff(sy, c(pn, names(data)))
    if (length(bad)) err("block_model(): unknown name(s) in '%s': %s (not data and not declared)", deparse1(f$decl), paste(bad, collapse = ", "))
    for (v in intersect(sy, c(pn, names(data)))) {
      vd <- if (v %in% pn) pars[[v]]$dim else ddim[[v]]
      if (!dim_fits(vd, f$shape))
        err("block_model(): '%s' uses %s (one value per %s) where one value per %s is needed", deparse1(f$decl), v, dim_word(vd), dim_word(f$shape))
    }
    f$vars <- union(intersect(sy, pn), if (!f$observed) f$y)
    if (f$dist == "normal" && (is.finite(f$lower) || is.finite(f$upper)) && length(intersect(c(all.vars(f$mean), all.vars(f$sd)), pn)))
      err("block_model(): '%s': a truncated prior needs a mean and sd without parameters (its normalizing constant would depend on them)", deparse1(f$decl))
    f$syms <- unique(c(sy, f$y))
    f$sd_vars <- if (!is.null(f$sd)) intersect(all.vars(f$sd), pn) else character()
    if (f$dist == "normal") f$r <- call("-", as.name(f$y), f$mean)
    if (f$dist == "bernoulli_logit") { Yd <- data[[f$y]]; f$mask <- !is.na(Yd); f$kappa <- ifelse(f$mask, Yd - 0.5, 0) }
    if (f$observed && f$dist == "normal") f$mask <- !is.na(data[[f$y]])
    expr <- if (f$dist == "normal") f$r else f$eta
    if (!is.null(expr)) {
      f$d1 <- list(); f$dsd <- list()
      for (v in f$vars) {
        f$d1[[v]] <- tryCatch(stats::D(expr, v), error = function(e)
          err("block_model(): '%s' uses a function that cannot be differentiated (use + - * / ^ exp log sqrt)", deparse1(f$decl)))
        if (f$dist == "normal") f$dsd[[v]] <- stats::D(f$sd, v)
      }
    }
    f$ir <- ml_factor(f, M)
    facs[[k]] <- f
  }
  M$factors <- facs
  M
}

# a self-contained declaration of factor f (definitions substituted, bounds and hyperparameters as numbers), to rebuild a model
bm_decl <- function(f, M) {
  lhs <- if (f$observed || M$pars[[f$y]]$dim == "scalar") as.name(f$y) else call("[", as.name(f$y), as.name(M$pars[[f$y]]$dim))
  rhs <- switch(f$dist, bernoulli_logit = call("bernoulli_logit", f$eta),
                normal = as.call(c(list(as.name("normal"), f$mean, f$sd), if (!f$observed) list(lower = f$lower, upper = f$upper))),
                as.call(c(list(as.name(f$dist)), as.list(unname(f$hyper)))))
  call("~", lhs, rhs)
}

dim_fits <- function(vd, shape) vd == "scalar" || vd == shape || shape == "NK"
dim_word <- function(d) switch(d, NK = "person and item", person = "person", item = "item", scalar = "model")
dim_len <- function(M, d) switch(d, person = M$N, item = M$K, scalar = 1L, NK = M$N * M$K)

#' @rdname block_model
#' @export
print.block_model <- function(x, ...) { print(equations(x)); invisible(x) }


# ---- evaluation helpers --------------------------------------------------------------------------

# the parameters and data of a factor, expanded to its shape
fenv <- function(f, s, M) {
  N <- M$N; K <- M$K
  bc <- function(x, d) if (f$shape == "NK") switch(d, person = matrix(x, N, K), item = matrix(x, N, K, byrow = TRUE), x) else x
  vals <- list(); syms <- f$syms %||% c(names(M$pars), names(M$data))       # only the names this part uses
  for (v in intersect(syms, names(M$pars))) if (!is.null(s[[v]])) vals[[v]] <- bc(s[[v]], M$pars[[v]]$dim)
  for (v in intersect(syms, names(M$data))) vals[[v]] <- bc(M$data[[v]], M$ddim[[v]])
  if (f$observed && !is.null(f$mask)) { yv <- vals[[f$y]]; yv[!f$mask] <- 0; vals[[f$y]] <- yv }
  list2env(vals, parent = baseenv())
}
full_of <- function(x, f, M) {
  if (length(x) != 1) return(x)
  switch(f$shape, NK = matrix(x, M$N, M$K), person = rep(x, M$N), item = rep(x, M$K), scalar = x)
}
reduce_to <- function(x, vdim, f) {
  if (f$shape == "NK") return(switch(vdim, person = rowSums(x), item = colSums(x), scalar = sum(x)))
  if (vdim == "scalar") sum(x) else x
}
fmask <- function(f, M) if (is.null(f$mask)) 1 else f$mask
is_zero <- function(e) identical(e, 0) || identical(e, 0L) || (is.numeric(e) && length(e) == 1 && e == 0)

# the factors that involve parameter v
factors_of <- function(M, v) Filter(function(f) v %in% f$vars, M$factors)

# how parameter v can be updated by a Gibbs step: "normal" (Gaussian full conditional),
# "variance" (conjugate inverse gamma), or a reason why not
gibbs_kind <- function(M, v) {
  P <- M$pars[[v]]
  if (P$dist %in% c("half_t", "inv_gamma")) {
    for (f in factors_of(M, v)) {
      if (!f$observed && f$y == v) next
      if (f$dist != "normal" || !identical(f$sd, call("sqrt", as.name(v))) || v %in% all.vars(f$r))
        return(sprintf("%s is used in '%s', not only as sd = sqrt(%s) of a normal", v, deparse1(f$decl), v))
    }
    return("variance")
  }
  for (f in factors_of(M, v)) {
    if (f$dist == "normal" && v %in% f$sd_vars) return(sprintf("the sd in '%s' depends on %s", deparse1(f$decl), v))
    if (!is_zero(tryCatch(stats::D(f$d1[[v]], v), error = function(e) 1)))
      return(sprintf("'%s' is not linear in %s", deparse1(f$decl), v))
  }
  "normal"
}

# ---- steps -------------------------------------------------------------------------------------------

#' Steps for declared models
#'
#' Steps name the parameters they update; [gibbs()] runs them in the order given, once per
#' iteration. What a step does is derived from the [block_model()]:
#'
#' * `step_gibbs(...)`: draws each parameter from its full conditional, one after the other. A
#'   parameter with a normal prior gets its normal full conditional: the precision and linear term
#'   are added up over every part of the model that uses it (normal parts directly; a logistic part
#'   through its Pólya–Gamma variables, which are drawn automatically when needed). A variance with
#'   a `half_t()` or `inv_gamma()` prior gets its conjugate update. The model must be linear in the
#'   parameter (otherwise use `step_mala()`).
#'   With `cut = "logT"` the parts of the data named there do not inform this step (a cut; Plummer,
#'   2015): e.g. `step_gibbs(theta, cut = "logT")` draws the ability from the responses only, while the
#'   response-time parameters still use the ability. The chain then targets the cut distribution, not
#'   the posterior (so the Geweke check of [check_sampler()] does not apply), and `ecm()` and
#'   `vi()` refuse a cut.
#' * `step_mala(..., step = 0.1)`: one MALA step (Roberts & Tweedie, 1996) per person or item, for
#'   the named parameters jointly (they must have the same dimension); the logistic parts use their
#'   exact likelihood (no Pólya–Gamma variables). A parameter with a lower bound is moved on the log
#'   scale with the Jacobian; the step sizes adapt during burn-in.
#' * `step_shift(...)`: moves several parameters together along a direction, e.g.
#'   `step_shift(theta = 1, b = 1)` (theta + c and b + c) or `step_shift(rho = -1, zeta = theta)`,
#'   with c drawn from its exact conditional (a generalized Gibbs step; Liu & Sabatti, 2000). The
#'   logistic parts must not change along the direction, and the normal parts must be linear in c.
#'   It speeds up directions that only the prior identifies.
#'
#' @param ... Parameter names (unquoted); for `step_shift()` `name = direction`.
#' @param step Starting step size of MALA.
#' @param cut Names of observed data whose parts are left out of this Gibbs step (a cut).
#' @return A step (a list) for [gibbs()].
#' @name block_steps
NULL

#' @rdname block_steps
#' @export
step_gibbs <- function(..., cut = character()) structure(list(type = "gibbs", vars = vapply(as.list(substitute(list(...)))[-1], deparse1, ""),
                                                             cut = as.character(cut)), class = "block_step")
#' @rdname block_steps
#' @export
step_mala <- function(..., step = 0.1) structure(list(type = "mala", vars = vapply(as.list(substitute(list(...)))[-1], deparse1, ""), step = step), class = "block_step")
#' @rdname block_steps
#' @export
step_shift <- function(...) {
  a <- as.list(substitute(list(...)))[-1]
  if (is.null(names(a)) || any(!nzchar(names(a)))) stop("step_shift(): give name = direction for each parameter, e.g. step_shift(theta = 1, b = 1)", call. = FALSE)
  structure(list(type = "shift", vars = names(a), along = a), class = "block_step")
}

# check the steps against the model; fill in defaults
prepare_steps <- function(M, steps) {
  pn <- names(M$pars)
  if (is.null(steps)) {
    kinds <- vapply(pn, function(v) gibbs_kind(M, v), "")
    g <- pn[kinds %in% c("normal", "variance")]; o <- setdiff(pn, g)
    steps <- c(if (length(g)) list(do.call(step_gibbs, lapply(g, as.name))),
               lapply(split(o, vapply(o, function(v) M$pars[[v]]$dim, "")), function(vs) do.call(step_mala, lapply(vs, as.name))))
  }
  if (inherits(steps, "block_step")) steps <- list(steps)
  if (!all(vapply(steps, inherits, TRUE, "block_step"))) stop("steps must be a list of step_gibbs(), step_mala() and step_shift()", call. = FALSE)
  for (k in seq_along(steps)) {
    st <- steps[[k]]; bad <- setdiff(st$vars, pn)
    if (length(bad)) stop(sprintf("step %d (%s): not a parameter of the model: %s", k, st$type, paste(bad, collapse = ", ")), call. = FALSE)
    if (!length(st$vars)) stop(sprintf("step %d (%s) names no parameter", k, st$type), call. = FALSE)
    if (st$type == "gibbs" && length(bad <- setdiff(st$cut, M$used)))
      stop(sprintf("step %d: cut = %s is not observed data of the model", k, paste(bad, collapse = ", ")), call. = FALSE)
    if (st$type == "gibbs") for (v in st$vars) {
      kd <- gibbs_kind(M, v)
      if (!kd %in% c("normal", "variance")) stop(sprintf("step %d: no Gibbs step for %s: %s; use step_mala(%s)", k, v, kd, v), call. = FALSE)
    }
    if (st$type == "mala") {
      d <- unique(vapply(st$vars, function(v) M$pars[[v]]$dim, ""))
      if (length(d) > 1) stop(sprintf("step %d: step_mala() moves parameters of one dimension together (%s); use separate steps", k, paste(st$vars, collapse = ", ")), call. = FALSE)
      for (v in st$vars) {
        P <- M$pars[[v]]
        if (P$dist != "normal") stop(sprintf("step %d: step_mala() is for parameters with a normal prior; %s is a variance (use step_gibbs(%s))", k, v, v), call. = FALSE)
        if (is.finite(P$upper)) stop(sprintf("step %d: step_mala() does not handle the upper bound of %s", k, v), call. = FALSE)
      }
      steps[[k]]$dim <- d
    }
    if (st$type == "shift") {
      for (v in st$vars) {
        dv <- all.vars(st$along[[v]])
        if (any(dv %in% st$vars)) stop(sprintf("step %d: the direction of %s uses a shifted parameter", k, v), call. = FALSE)
        if (any(!dv %in% c(pn, names(M$data)))) stop(sprintf("step %d: unknown name in the direction of %s", k, v), call. = FALSE)
        if (M$pars[[v]]$dist != "normal") stop(sprintf("step %d: step_shift() moves parameters with a normal prior; %s is a variance", k, v), call. = FALSE)
      }
      for (f in M$factors) {
        sv <- intersect(f$vars, st$vars); if (!length(sv)) next
        if (f$dist == "normal" && length(intersect(f$sd_vars, st$vars))) stop(sprintf("step %d: the sd in '%s' depends on a shifted parameter", k, deparse1(f$decl)), call. = FALSE)
        for (u in sv) for (w in sv) if (!is_zero(tryCatch(stats::D(f$d1[[u]], w), error = function(e) 1)))
          stop(sprintf("step %d: '%s' is not linear along the shift (%s, %s); use step_mala()", k, deparse1(f$decl), u, w), call. = FALSE)
      }
    }
  }
  upd <- unique(unlist(lapply(Filter(function(s) s$type != "shift", steps), `[[`, "vars")))
  miss <- setdiff(pn, upd)
  if (length(miss)) stop(sprintf("no step updates %s (step_shift() moves only one direction); add step_gibbs() or step_mala()", paste(miss, collapse = ", ")), call. = FALSE)
  steps
}

# the logistic parts whose Polya-Gamma variables a Gibbs step for v uses
pg_factors <- function(M, v) Filter(function(f) f$dist == "bernoulli_logit" && v %in% f$vars, M$factors)

draw_omega_f <- function(f, s, M) {
  om <- draw_pg(fval(f, s, M, "eta")); om[!f$mask] <- 0; om
}

gibbs_normal <- function(v, s, M, cut = character()) {
  P <- M$pars[[v]]; prec <- 0; num <- 0
  s0 <- s; s0[[v]] <- 0 * s[[v]]                                # the rest of the expression: its value at v = 0
  for (f in factors_of(M, v)) {
    if (f$observed && f$y %in% cut) next                         # a cut: this data part does not inform v
    fast <- ml_part_f(f, v, s, M, if (f$dist == "normal") "normal" else "gibbs", om = if (f$dist == "bernoulli_logit") s$.omega[[f$id]])
    if (!is.null(fast)) { prec <- prec + fast[[1]]; num <- num + fast[[2]]; next }
    m <- fmask(f, M); c1 <- fval(f, s, M, "d1", v)
    if (f$dist == "normal") {
      w <- m / fval(f, s, M, "sd2"); r0 <- fval(f, s0, M, "r")
      prec <- prec + reduce_to(w * c1^2, P$dim, f); num <- num - reduce_to(w * c1 * r0, P$dim, f)
    } else {
      om <- s$.omega[[f$id]]; eta0 <- fval(f, s0, M, "eta")
      prec <- prec + reduce_to(om * c1^2, P$dim, f); num <- num + reduce_to(c1 * (f$kappa - om * eta0), P$dim, f)
    }
  }
  draw(normal_part(prec, num), lower = P$lower, upper = P$upper)
}

gibbs_variance <- function(v, s, M, cut = character()) {
  P <- M$pars[[v]]; n <- dim_len(M, P$dim); al <- numeric(n); be <- numeric(n)
  for (f in factors_of(M, v)) {
    if (!f$observed && f$y == v) next
    if (f$observed && f$y %in% cut) next
    m <- fmask(f, M); r <- fval(f, s, M, "r")
    al <- al + reduce_to(m * full_of(1, f, M), P$dim, f) / 2; be <- be + reduce_to(m * r^2, P$dim, f) / 2
  }
  if (P$dist == "inv_gamma") return(draw_invgamma(P$hyper[1] + al, P$hyper[2] + be))
  vapply(seq_len(n), function(k) draw_var(al[k], be[k], s[[v]][k], prior = c(df = unname(P$hyper[1]), scale = unname(P$hyper[2]))), 0)
}

# log posterior and gradient of the parameters vs (one dimension) on the MALA scale
mala_target <- function(vs, s, M) {
  Ps <- M$pars[vs]; d <- Ps[[1]]$dim
  lo <- vapply(Ps, `[[`, 0, "lower")
  facs <- Filter(function(f) length(intersect(f$vars, vs)) > 0, M$factors)
  function(x) {
    x <- matrix(x, ncol = length(vs)); s2 <- s
    for (k in seq_along(vs)) s2[[vs[k]]] <- if (is.finite(lo[k])) lo[k] + exp(x[, k]) else x[, k]
    lp <- 0; G <- matrix(0, nrow(x), length(vs))
    for (f in facs) {
      m <- fmask(f, M)
      if (f$dist == "normal") {
        r <- fval(f, s2, M, "r"); v2 <- fval(f, s2, M, "sd2")
        lp <- lp + reduce_to(m * (-r^2 / (2 * v2) - log(v2) / 2), d, f)
        for (k in seq_along(vs)) if (vs[k] %in% f$vars) G[, k] <- G[, k] + reduce_to(m * gnorm(f, s2, M, vs[k], r, v2), d, f)
      } else {
        eta <- fval(f, s2, M, "eta"); y <- fval(f, s2, M, "y")
        lp <- lp + reduce_to(m * (y * eta - softplus(eta)), d, f)
        for (k in seq_along(vs)) if (vs[k] %in% f$vars)
          G[, k] <- G[, k] + reduce_to(m * (y - stats::plogis(eta)) * fval(f, s2, M, "d1", vs[k]), d, f)
      }
    }
    for (k in seq_along(vs)) if (is.finite(lo[k])) { lp <- lp + x[, k]; G[, k] <- G[, k] * exp(x[, k]) + 1 }
    list(lp = lp, grad = G)
  }
}

do_shift <- function(st, s, M) {
  dir <- list()
  for (v in st$vars) {                                        # direction in the dimension of v
    P <- M$pars[[v]]; e <- list2env(c(s[names(M$pars)], M$data), parent = baseenv())
    dv <- eval(st$along[[v]], e); n <- dim_len(M, P$dim)
    if (!length(dv) %in% c(1, n)) stop(sprintf("step_shift(): the direction of %s needs 1 or %d values", v, n), call. = FALSE)
    dir[[v]] <- rep_len(as.numeric(dv), n)
  }
  prec <- 0; num <- 0; lo <- -Inf; hi <- Inf
  for (f in M$factors) {
    sv <- intersect(f$vars, st$vars); if (!length(sv)) next
    m <- fmask(f, M)
    bc <- function(x, v) { d <- M$pars[[v]]$dim; if (f$shape == "NK") switch(d, person = matrix(x, M$N, M$K), item = matrix(x, M$N, M$K, byrow = TRUE), x) else x }
    dc <- 0; for (v in sv) dc <- dc + fval(f, s, M, "d1", v) * bc(dir[[v]], v)
    if (f$dist == "bernoulli_logit") {
      if (any(abs(m * dc) > 1e-8)) stop(sprintf("step_shift(): the shift changes the logistic part '%s'; choose a direction that leaves it unchanged", deparse1(f$decl)), call. = FALSE)
      next
    }
    w <- m / fval(f, s, M, "sd2"); r <- fval(f, s, M, "r")
    prec <- prec + sum(w * dc^2); num <- num - sum(w * dc * r)
  }
  for (v in st$vars) {
    P <- M$pars[[v]]; d <- dir[[v]]; x <- s[[v]]
    for (k in which(d > 0)) { lo <- max(lo, (P$lower - x[k]) / d[k]); hi <- min(hi, (P$upper - x[k]) / d[k]) }
    for (k in which(d < 0)) { lo <- max(lo, (P$upper - x[k]) / d[k]); hi <- min(hi, (P$lower - x[k]) / d[k]) }
  }
  if (!(prec > 0)) stop("step_shift(): the model gives no information along the shift", call. = FALSE)
  c0 <- draw_tnorm(num / prec, 1 / sqrt(prec), lo, hi)
  for (v in st$vars) s[[v]] <- s[[v]] + c0 * dir[[v]]
  s
}

# one iteration of the stacked steps; Polya-Gamma variables are redrawn when stale
model_step <- function(steps, M) function(s, d) {
  if (!is.null(d$data) && !identical(d$data, M$data)) M <- model_data(M, d$data)
  stale <- rep(TRUE, length(M$factors))                         # omega is drawn once per iteration ...
  for (k in seq_along(steps)) {
    st <- steps[[k]]
    if (st$type == "gibbs") for (v in st$vars) {
      for (f in pg_factors(M, v)) if (stale[f$id] && !(f$y %in% st$cut)) { s$.omega[[f$id]] <- draw_omega_f(M$factors[[f$id]], s, M); stale[f$id] <- FALSE; s$.n_omega <- s$.n_omega + 1 }
      s[[v]] <- if (M$pars[[v]]$dist == "normal") gibbs_normal(v, s, M, st$cut) else gibbs_variance(v, s, M, st$cut)
    }
    if (st$type == "mala") {
      lo <- vapply(M$pars[st$vars], `[[`, 0, "lower")
      x <- vapply(seq_along(st$vars), function(q) { z <- s[[st$vars[q]]]; if (is.finite(lo[q])) log(z - lo[q]) else z }, numeric(dim_len(M, st$dim)))
      x <- matrix(x, ncol = length(st$vars))
      x <- draw_mala(x, mala_target(st$vars, s, M), s$.eps[[k]])
      acc <- attr(x, "accepted")
      if (isTRUE(s$.burn)) s$.eps[[k]] <- mcmc_adapt(s$.eps[[k]], acc, s$.iter)
      for (q in seq_along(st$vars)) s[[st$vars[q]]] <- if (is.finite(lo[q])) lo[q] + exp(x[, q]) else x[, q]
      for (f in M$factors) if (f$dist == "bernoulli_logit" && any(st$vars %in% f$vars)) stale[f$id] <- TRUE   # ... and after MALA moved it
    }
    if (st$type == "shift") s <- do_shift(st, s, M)                # the logistic parts are unchanged: omega stays fresh
  }
  s
}

# the model with new data of the same shape (simulated data for check_sampler())
model_data <- function(M, data) {
  M$data <- data; M$data0 <- lapply(data, function(x) { x[is.na(x)] <- 0; x })
  for (f in M$factors) if (f$observed) {
    M$factors[[f$id]]$mask <- !is.na(data[[f$y]])
    if (f$dist == "bernoulli_logit") M$factors[[f$id]]$kappa <- ifelse(M$factors[[f$id]]$mask, data[[f$y]] - 0.5, 0)
  }
  M
}

init_model_state <- function(M, steps, init = NULL, from_data = TRUE) {
  s <- list()
  for (v in names(M$pars)) {
    P <- M$pars[[v]]; n <- dim_len(M, P$dim)
    s[[v]] <- if (!is.null(init[[v]])) rep_len(as.numeric(init[[v]]), n) else if (P$dist == "normal") {
      m0 <- if (is.numeric(P$mean) && length(P$mean) == 1) P$mean else 0
      if (is.finite(P$lower) && m0 <= P$lower) m0 <- P$lower + 1
      if (is.finite(P$upper) && m0 >= P$upper) m0 <- P$upper - 1
      x <- m0 + stats::rnorm(n, 0, 0.1)
      if (is.finite(P$lower)) x <- pmax(x, P$lower + 0.05)
      if (is.finite(P$upper)) x <- pmin(x, P$upper - 0.05)
      x
    } else rep(1, n)
  }
  if (from_data) s <- data_start(M, s, names(init))
  s$.omega <- vector("list", length(M$factors)); s$.n_omega <- 0
  s$.eps <- lapply(steps, function(st) if (st$type == "mala") rep(st$step, dim_len(M, st$dim)) else NULL)
  s
}

# Starting values from the data for intercept-like parameters (those whose derivative does not change with
# the person parameters, e.g. lambda in lambda - zeta - rho * theta, b in a * (theta - b)): with the person
# parameters at their prior means, move each one so that the mean residual of a normal part is 0, or the
# linear predictor of a logistic part is the logit of the observed proportion. One parameter per part.
# Far starting values (e.g. lambda at its prior mean 0 for log times near 4) make EM crawl along the
# direction that the person parameters can absorb.
data_start <- function(M, s, given = character()) {
  person <- names(M$pars)[vapply(M$pars, function(P) P$dim == "person", TRUE)]
  s0 <- s; for (v in person) s0[[v]] <- 0 * s[[v]]
  s1 <- s; for (v in person) s1[[v]] <- s[[v]] * 0 + 1
  for (f in M$factors) {
    if (!f$observed || !(f$dist %in% c("normal", "bernoulli_logit"))) next
    key <- if (f$dist == "normal") "r" else "eta"
    cand <- setdiff(intersect(f$vars, names(M$pars)), c(person, given))
    for (v in cand) {
      P <- M$pars[[v]]
      if (P$dist != "normal" || !(P$dim %in% c("item", "scalar"))) next
      d0 <- tryCatch(as.matrix(fval(f, s0, M, "d1", v)), error = function(e) NULL)
      d1 <- tryCatch(as.matrix(fval(f, s1, M, "d1", v)), error = function(e) NULL)
      if (is.null(d0) || is.null(d1) || any(!is.finite(d0)) || max(abs(d0 - d1)) > 1e-10 || any(abs(d0) < 1e-8)) next
      y <- as.matrix(fval(f, s0, M, "y")); x <- as.matrix(fval(f, s0, M, key))
      m <- !is.na(y) & as.matrix(fmask(f, M)) > 0                     # observed cells only
      x[!m] <- NA; y[!m] <- NA
      dc <- colMeans(d0)
      delta <- if (f$dist == "normal") -colMeans(x, na.rm = TRUE) / dc else {
        p <- pmin(pmax(colMeans(y, na.rm = TRUE), 0.02), 0.98)
        (stats::qlogis(p) - colMeans(x, na.rm = TRUE)) / dc
      }
      if (P$dim == "scalar") delta <- mean(delta)
      x <- s[[v]] + delta
      if (is.finite(P$lower)) x <- pmax(x, P$lower + 0.05)
      if (is.finite(P$upper)) x <- pmin(x, P$upper - 0.05)
      if (all(is.finite(x))) s[[v]] <- x
      break                                                            # one intercept per part
    }
  }
  s
}

#' Sample a declared model with stacked steps
#'
#' Runs the steps once per iteration in the order given (see [block_steps]) with [run_sampler()]
#' (chains, burn-in, storage, summary). Without `steps`, every parameter that has a Gibbs step gets
#' one (in the order of the declarations) and the others get MALA.
#'
#' @param model A [block_model()].
#' @param steps A list of steps, e.g. `list(step_mala(theta), step_gibbs(a, b))`.
#' @param n_iter,n_burnin,n_chain,n_thin,seed As in [run_sampler()].
#' @param init Optional starting values: a named list (one value or one per person/item). The
#'   default starts near the prior means; with two latent variables give data-based values (e.g.
#'   standardized sum scores and mean log times), since a chain can settle in another mode (the
#'   ability taking over the role of speed); `print(fit)` warns when the chains disagree.
#' @param monitor Parameters to store; by default all but the person parameters, whose posterior
#'   means and SDs are kept instead (`fit$mean`).
#' @return A `block_fit` (see [run_sampler()]) with the model and the steps; `algorithm(fit)` lists
#'   the steps with their full conditionals.
#' @export
gibbs.block_model <- function(model, steps = NULL, n_iter = 2000, n_burnin = floor(n_iter / 2), n_chain = 2, n_thin = 1, seed = 1,
                          init = NULL, monitor = NULL, ...) {
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  if (!inherits(model, "block_model")) stop("model must come from block_model()", call. = FALSE)
  steps <- prepare_steps(model, steps)
  pn <- names(model$pars)
  person <- pn[vapply(pn, function(v) model$pars[[v]]$dim == "person", TRUE)]
  monitor <- monitor %||% setdiff(pn, person)
  keep <- setdiff(person, monitor)
  fit <- run_sampler(model, function(d) init_model_state(model, steps, init), model_step(steps, model), monitor = monitor,
                     n_iter = n_iter, n_burnin = n_burnin, n_chain = n_chain, n_thin = n_thin, seed = seed,
                     keep_mean = if (length(keep)) keep)
  fit$model <- model; fit$steps <- steps
  fit$n_omega <- vapply(fit$last, function(s) s$.n_omega, 0)
  class(fit) <- c("block_model_fit", class(fit))
  fit
}


# ---- prior simulation and the Geweke check for declared models --------------------------------------

# draw the parameters from their priors and the data given them (declarations in dependency order)
simulate_model <- function(M, what = c("parameters", "data"), s = NULL) {
  what <- match.arg(what)
  if (what == "parameters") {
    s <- list(); todo <- names(M$pars)
    while (length(todo)) {
      ready <- todo[vapply(todo, function(v) all(intersect(c(all.vars(M$pars[[v]]$mean %||% 0), all.vars(M$pars[[v]]$sd %||% 0)), names(M$pars)) %in% names(s)), TRUE)]
      if (!length(ready)) stop("the priors refer to each other in a circle", call. = FALSE)
      for (v in ready) {
        P <- M$pars[[v]]; n <- dim_len(M, P$dim)
        s[[v]] <- if (P$dist == "normal") {
          e <- list2env(c(s, M$data), parent = baseenv())
          mu <- rep_len(eval(P$mean, e), n); sdv <- rep_len(eval(P$sd, e), n)
          if (is.finite(P$lower) || is.finite(P$upper)) draw_tnorm(mu, sdv, P$lower, P$upper) else stats::rnorm(n, mu, sdv)
        } else if (P$dist == "half_t") (P$hyper[2] * abs(stats::rt(n, P$hyper[1])))^2 else 1 / stats::rgamma(n, P$hyper[1], P$hyper[2])
      }
      todo <- setdiff(todo, ready)
    }
    return(s)
  }
  data <- M$data
  for (f in Filter(function(f) f$observed, M$factors)) {
    e <- fenv(f, s, M)
    data[[f$y]] <- if (f$dist == "bernoulli_logit") { p <- stats::plogis(full_of(eval(f$eta, e), f, M)); matrix(stats::rbinom(length(p), 1, p), M$N, M$K) }
                   else { mu <- full_of(eval(f$mean, e), f, M); sdv <- full_of(eval(f$sd, e), f, M); matrix(stats::rnorm(length(mu), mu, sdv), M$N, M$K) }
    data[[f$y]][is.na(M$data[[f$y]])] <- NA
  }
  data
}

check_model_steps <- function(M, steps, n_iter, probs, n_prior, monitor) {
  steps <- prepare_steps(M, steps)
  pn <- names(M$pars)
  monitor <- monitor %||% pn[vapply(pn, function(v) M$pars[[v]]$dim != "person", TRUE)]
  prior <- function() { s <- simulate_model(M); s0 <- init_model_state(M, steps); s0[pn] <- s[pn]; s0$.iter <- 0; s0 }
  one <- model_step(steps, M)
  step <- function(s, d) { s$.iter <- s$.iter + 1; s$.burn <- s$.iter <= n_iter / 10; one(s, d) }   # MALA adapts in the dropped 10%
  check_step(step, prior, function(s) list(data = simulate_model(M, "data", s)), monitor = monitor,
               n_iter = n_iter, n_prior = n_prior, probs = probs)
}
