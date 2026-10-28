# The built-in models by name, as in mirt and TAM:
#   fit_ecm <- ecm(dat, model = "cross")
#   fit_gibbs <- gibbs(dat, model = "cross", init = fit_ecm)
# `model` picks the structure and `itemtype` the item model; the engine settings that are rarely
# changed live in ecm_control() and gibbs_control(). Both build the model object of the old
# constructors (new_model()) and run the same engine, so old and new code give the same fit.

# model names -> internal model types (NA: not available yet)
model_types <- c(irt = "mlirt", null = "rtirt_null", latreg = "rtirt_latreg",
                 latent = "rtirt_latent", cross = "rtirt_cross", mixture = NA)

# one line per available model, for the messages and ?ecm
model_help <- c(irt = "2PL, no response times (with covariates: ability regressed on them)",
                null = "responses and times; ability and speed correlated (van der Linden, 2007)",
                latreg = "as \"null\", with ability and speed regressed on the covariates",
                latent = "speed regressed on ability and the covariates",
                cross = "item-specific effects of ability on the log times (cross-relations)")

available_models <- function() names(model_types)[!is.na(model_types)]

# the user-facing name of a model type: "rtirt_cross" -> "cross"
model_name_of <- function(type) {
  name <- names(model_types)[match(type, model_types)]
  if (is.na(name)) type else name
}

# 'model = "cross"', for headers and messages
model_call_text <- function(type) sprintf("model = \"%s\"", model_name_of(type))

model_list_text <- function() {
  available <- available_models()
  padded <- formatC(sprintf("\"%s\"", available), width = -max(nchar(available) + 2))
  paste(sprintf("  %s %s", padded, model_help[available]), collapse = "\n")
}

# the name given as model = : a public name (any case), or an error with the list of models and
# a did-you-mean
model_name <- function(model) {
  available <- available_models()
  if (is.null(model)) {
    stop(sprintf("name the model: model = one of\n%s", model_list_text()), call. = FALSE)
  }
  if (!is.character(model) || length(model) != 1 || is.na(model)) {
    stop(sprintf("model must be one name: one of\n%s", model_list_text()), call. = FALSE)
  }
  key <- tolower(trimws(model))
  if (key == "mixture") {
    stop("model = \"mixture\" (Wang & Xu, 2015) is not available yet", call. = FALSE)
  }
  if (key %in% available) {
    return(key)
  }
  stop(sprintf("unknown model \"%s\"%s\nmodel = one of\n%s", model, model_guess(key, available), model_list_text()),
       call. = FALSE)
}

# "; did you mean "cross"?" for a misspelled model, the public name for a class name ("mlirt"),
# and a pointer to itemtype = for "2pl"
model_guess <- function(key, available) {
  class_name <- match(key, model_types)
  if (!is.na(class_name)) {
    return(sprintf(": that is the class of model = \"%s\"; did you mean \"%s\"?", names(model_types)[class_name],
                   names(model_types)[class_name]))
  }
  if (key %in% c("2pl", "1pl", "rasch", "grm")) {
    return(sprintf(": that is an item type (itemtype = \"%s\"); model names the structure, e.g. model = \"irt\" for the 2PL alone",
                   if (key == "rasch") "1pl" else key))
  }
  guess <- closest_name(key, available)
  if (is.null(guess)) "" else sprintf("; did you mean \"%s\"?", guess)
}

check_itemtype <- function(itemtype, n_item) {
  allowed <- c("2pl", "1pl")
  if (is.numeric(itemtype)) {
    stop(sprintf("itemtype = %s: the third argument of ecm() and gibbs() is itemtype (\"2pl\" or \"1pl\"); name the arguments after model, e.g. n_iter = %s",
                 itemtype[1], itemtype[1]), call. = FALSE)
  }
  right_length <- length(itemtype) %in% c(1, n_item)
  if (!is.character(itemtype) || !right_length) {
    stop(sprintf("itemtype must be one value, or one per item (%d)", n_item), call. = FALSE)
  }
  itemtype <- tolower(itemtype)
  if ("grm" %in% itemtype) {
    stop("itemtype = \"grm\" is not available yet", call. = FALSE)
  }
  unknown <- setdiff(itemtype, allowed)
  if (length(unknown)) {
    hint <- if ("rasch" %in% unknown) "; the Rasch model is itemtype = \"1pl\"" else ""
    stop(sprintf("unknown itemtype: %s (use \"2pl\" or \"1pl\")%s", paste(unknown, collapse = ", "), hint), call. = FALSE)
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
  if (!is.null(cov) && !is.character(cov)) {
    in_data <- colnames(data$X)
    have <- if (length(in_data)) sprintf(" (here: %s)", paste(in_data, collapse = ", ")) else " (here: none)"
    stop(sprintf("cov = takes the names of covariates already in the data%s, not their values; add a covariate with input_data(..., cov = ), e.g. input_data(resp = items, time = times, cov = \"age_c\", data = raw)",
                 have), call. = FALSE)
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

# quantile = for a model without a quantile version: an error, as ecm() gives for quantile =
check_quantile <- function(model, type, quantile) {
  if (is.null(quantile) || .model_info[[type]]$qr) {
    return(invisible())
  }
  with_qr <- names(model_types)[model_types %in% names(.model_info)[vapply(.model_info, `[[`, TRUE, "qr")]]
  stop(sprintf("quantile = is available for model = %s only; model = \"%s\" is a mean model (drop quantile)",
               paste0("\"", with_qr, "\"", collapse = " or "), model), call. = FALSE)
}

build_model <- function(data, model, itemtype, speed_var, cov, quantile = NULL, intercept = FALSE) {
  model <- model_name(model)
  type <- model_types[[model]]
  check_quantile(model, type, quantile)
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
#' @param ... Not used: a name that is not a setting is an error that names the closest
#'   setting, or the function it belongs to (`gibbs_control(n_iter = 200)`: an argument of
#'   [gibbs()]).
#' @return A list of settings of class `ecm_control` or `gibbs_control`; [ecm()] refuses a
#'   `gibbs_control()` and [gibbs()] an `ecm_control()`. A plain list with the names of the
#'   settings is also accepted.
#' @examples
#' fit_settings <- gibbs_control(n_burnin = 500, collapse = TRUE)
#' fit_settings
#' @export
ecm_control <- function(max_iter = 1000, tol = 1e-10, nodes = 41, accelerate = TRUE, verbose = FALSE, ...) {
  check_control_dots("ecm", list(...))
  settings <- list(max_iter = max_iter, tol = tol, nodes = nodes, accelerate = accelerate, verbose = verbose)
  structure(settings, class = "ecm_control")
}

#' @rdname ecm_control
#' @export
gibbs_control <- function(n_burnin = NULL, n_thin = 1, collapse = FALSE, fit_indices = TRUE, verbose = TRUE, ...) {
  check_control_dots("gibbs", list(...))
  settings <- list(n_burnin = n_burnin, n_thin = n_thin, collapse = collapse, fit_indices = fit_indices, verbose = verbose)
  structure(settings, class = "gibbs_control")
}

#' @rdname ecm_control
#' @param x Settings from `ecm_control()` or `gibbs_control()`.
#' @export
print.ecm_control <- function(x, ...) {
  print_control(x)
}

#' @rdname ecm_control
#' @export
print.gibbs_control <- function(x, ...) {
  print_control(x)
}

print_control <- function(x) {
  values <- vapply(unclass(x), function(v) if (is.null(v)) "NULL (default)" else paste(format(v), collapse = ", "), "")
  cat(sprintf("%s()\n", class(x)[1]))
  cat(sprintf("  %s = %s\n", names(values), values), sep = "")
  invisible(x)
}

# names in ecm_control(...) / gibbs_control(...) that are not settings
check_control_dots <- function(fun, dots) {
  if (!length(dots)) {
    return(invisible())
  }
  control_fun <- paste0(fun, "_control")
  other <- if (fun == "ecm") "gibbs" else "ecm"
  own_control <- setdiff(names(formals(engine_control(fun))), "...")
  own_args <- setdiff(names(formals(engine_method(fun))), c("data", "...", "control"))
  other_control <- setdiff(names(formals(engine_control(other))), c("...", own_control))
  other_args <- setdiff(names(formals(engine_method(other))), c("data", "...", "control", own_args, own_control))
  arg_names <- names(dots) %||% rep("", length(dots))
  lines <- vapply(arg_names, function(name) {
    if (!nzchar(name)) {
      return(sprintf("unnamed value in %s(): name the settings, e.g. %s(%s = ...)", control_fun, control_fun, own_control[1]))
    }
    if (name %in% own_args) {
      return(sprintf("%s is an argument of %s(), not a setting of %s(): %s(data, model = ..., %s = ...)", name, fun,
                     control_fun, fun, name))
    }
    if (name %in% other_control) {
      return(sprintf("%s is a setting of %s_control() for %s(), not of %s()", name, other, other, control_fun))
    }
    if (name %in% other_args) {
      return(sprintf("%s is an argument of %s(), not a setting of %s()", name, other, control_fun))
    }
    guess <- closest_name(name, own_control)
    if (is.null(guess)) {
      sprintf("unknown setting %s of %s(); the settings are %s", name, control_fun, paste(own_control, collapse = ", "))
    } else {
      sprintf("unknown setting %s of %s(); did you mean %s?", name, control_fun, guess)
    }
  }, "")
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

# control = of ecm() and gibbs(): settings of its own engine (a plain list goes through the
# control function, which checks the names)
as_control <- function(control, fun) {
  control_fun <- paste0(fun, "_control")
  other <- if (fun == "ecm") "gibbs" else "ecm"
  if (inherits(control, control_fun)) {
    return(control)
  }
  if (inherits(control, paste0(other, "_control"))) {
    stop(sprintf("control = %s_control() holds the settings of %s(); %s() takes control = %s(...)", other, other, fun,
                 control_fun), call. = FALSE)
  }
  if (is.list(control) && !is.object(control) && (!length(control) || !is.null(names(control)))) {
    return(do.call(engine_control(fun), control))
  }
  stop(sprintf("control must come from %s(), e.g. %s(dat, model = \"cross\", control = %s(...))", control_fun, fun,
               control_fun), call. = FALSE)
}

# ---- arguments that ecm() and gibbs() do not have ------------------------------------------------
# The `...` of ecm(data, ...) and gibbs(data, ...) takes nothing. Before anything is fitted, an
# engine setting written directly (ecm(dat, model = "cross", nodes = 61)) is an error that shows
# the control = form, an argument of the other engine says so, and other names get the closest
# argument (did you mean ...?).

engine_method <- function(fun) if (fun == "ecm") ecm.input_data else gibbs.input_data
engine_control <- function(fun) if (fun == "ecm") ecm_control else gibbs_control

# the closest of `candidates` to a misspelled name (n_iter also matches "iter"), or NULL
closest_name <- function(name, candidates) {
  squash <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  variants <- c(candidates, sub("^n_", "", candidates))
  owners <- c(candidates, candidates)
  distance <- utils::adist(squash(name), squash(variants))[1, ]
  best <- which.min(distance)
  if (distance[best] <= max(1, floor(nchar(name) / 4))) owners[best] else NULL
}

check_dots <- function(fun, call, dots) {
  if (!length(dots)) {
    return(invisible())
  }
  other <- if (fun == "ecm") "gibbs" else "ecm"
  control_fun <- paste0(fun, "_control")
  own_args <- setdiff(names(formals(engine_method(fun))), c("data", "..."))
  own_control <- setdiff(names(formals(engine_control(fun))), "...")
  other_args <- setdiff(names(formals(engine_method(other))), c("data", "...", own_args, own_control))
  other_control <- setdiff(names(formals(engine_control(other))), c("...", own_args, own_control))
  arg_names <- names(dots) %||% rep("", length(dots))
  values <- vapply(dots, function(v) paste(deparse(v, width.cutoff = 500L), collapse = " "), "")
  lines <- character()
  to_control <- character()                                  # name = value pairs for control = ...()
  for (i in seq_along(dots)) {
    name <- arg_names[i]
    if (!nzchar(name)) {
      lines <- c(lines, sprintf("unnamed argument %s: name the arguments after data, e.g. model = \"cross\"", values[i]))
    } else if (name %in% own_control) {
      to_control <- c(to_control, sprintf("%s = %s", name, values[i]))
    } else if (name %in% other_args) {
      lines <- c(lines, sprintf("%s is an argument of %s(), not of %s()", name, other, fun))
    } else if (name %in% other_control) {
      lines <- c(lines, sprintf("%s is a setting of %s_control() for %s(), not of %s()", name, other, other, fun))
    } else {
      guess <- closest_name(name, c(own_args, own_control))
      if (is.null(guess)) {
        lines <- c(lines, sprintf("unknown argument %s; the arguments of %s() are %s (settings in control = %s())", name, fun,
                                  paste(own_args, collapse = ", "), control_fun))
      } else if (guess %in% own_control) {
        lines <- c(lines, sprintf("unknown argument %s; did you mean %s (a setting of %s())?", name, guess, control_fun))
        to_control <- c(to_control, sprintf("%s = %s", guess, values[i]))
      } else {
        lines <- c(lines, sprintf("unknown argument %s; did you mean %s?", name, guess))
      }
    }
  }
  if (length(to_control)) {
    direct <- arg_names[arg_names %in% own_control]
    hint <- control_hint(fun, call, to_control, direct)
    lines <- if (length(direct)) c(hint, lines) else c(lines, hint)
  }
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

# "nodes is a setting of ecm_control(): use ecm(dat, model = "cross", control = ecm_control(nodes = 61))",
# keeping the settings of a control = ecm_control(...) that was given
control_hint <- function(fun, call, to_control, direct) {
  control_fun <- paste0(fun, "_control")
  given <- call$control
  if (is.call(given) && identical(given[[1]], as.name(control_fun)) && length(given) > 1) {
    kept <- as.list(given)[-1]
    kept <- kept[!names(kept) %in% sub(" = .*", "", to_control)]
    kept_text <- vapply(kept, function(v) paste(deparse(v, width.cutoff = 500L), collapse = " "), "")
    to_control <- c(sprintf("%s = %s", names(kept), kept_text), to_control)
  }
  deparse_arg <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  data_text <- if (is.null(call$data)) "data" else deparse_arg(call$data)
  model_text <- if (is.null(call$model)) "" else sprintf(", model = %s", deparse_arg(call$model))
  usage <- sprintf("%s(%s%s, control = %s(%s))", fun, data_text, model_text, control_fun, paste(to_control, collapse = ", "))
  if (!length(direct)) {
    return(sprintf("use %s", usage))
  }
  what <- if (length(direct) == 1) "is a setting of %s(), not an argument" else "are settings of %s(), not arguments"
  sprintf("%s %s of %s(): use %s", paste(direct, collapse = ", "), sprintf(what, control_fun), fun, usage)
}

#' @rdname ecm
#' @export
ecm.input_data <- function(data, model = NULL, itemtype = "2pl", speed_var = c("free", "fixed"), cov = NULL,
                           prior = c("ml", "map"), priors = NULL, se = TRUE, control = ecm_control(), ...) {
  if ("quantile" %in% names(list(...))) {
    stop("ecm() fits the normal models; quantile models use gibbs(data, model = ..., quantile = )", call. = FALSE)
  }
  check_dots("ecm", match.call(), match.call(expand.dots = FALSE)$...)
  control <- as_control(control, "ecm")
  speed_var <- match.arg(speed_var)
  object <- build_model(data, model, itemtype, speed_var, cov)
  ecm.rtirt(object, max_iter = control$max_iter, tol = control$tol, nodes = control$nodes, se = se,
            prior = prior, priors = priors, accelerate = control$accelerate, verbose = control$verbose, ...)
}

# an ecm() fit as the start of gibbs(): it must be the same data and model; the error says what
# differs
attach_ecm <- function(object, fit) {
  problems <- init_differences(object, fit)
  if (length(problems)) {
    type <- class(object)[1]
    stop(sprintf("%s; init must be an ecm() fit of the same data and model: init = ecm(data, %s%s), or init = \"random\"",
                 paste(problems, collapse = "; "), model_call_text(type), init_settings_text(object)), call. = FALSE)
  }
  object$ecm <- fit$ecm
  object
}

init_differences <- function(object, fit) {
  type <- class(object)[1]
  if (!identical(class(fit)[1], type)) {
    return(sprintf("init is an ecm() fit of %s, but gibbs() was asked for %s", model_call_text(class(fit)[1]),
                   model_call_text(type)))
  }
  data <- object$data
  fit_data <- fit$data
  problems <- character()
  if (!identical(dim(fit_data$Y), dim(data$Y))) {
    problems <- c(problems, sprintf("init was fitted to %d persons x %d items, gibbs() has %d x %d", nrow(fit_data$Y),
                                    ncol(fit_data$Y), nrow(data$Y), ncol(data$Y)))
  } else {
    if (!identical(fit_data$Y, data$Y)) problems <- c(problems, "init was fitted to other responses")
    if (!identical(fit_data$log_t, data$log_t)) problems <- c(problems, "init was fitted to other response times")
  }
  fit_cov <- colnames(fit_data$X) %||% character()
  cov <- colnames(data$X) %||% character()
  if (!identical(fit_cov, cov)) {
    cov_text <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
    problems <- c(problems, sprintf("init was fitted with the covariates %s, gibbs() uses %s", cov_text(fit_cov), cov_text(cov)))
  } else if (length(problems) == 0 && !identical(fit_data$X, data$X)) {
    problems <- c(problems, "init was fitted to other covariate values")
  }
  for (key in c("speed_var", "itemtype")) {
    if (!identical(fit$settings[[key]], object$settings[[key]])) {
      problems <- c(problems, sprintf("init was fitted with %s = \"%s\", gibbs() was asked for %s = \"%s\"", key,
                                      fit$settings[[key]], key, object$settings[[key]]))
    }
  }
  problems
}

# the settings of the model that differ from the defaults, for the init = ecm(...) call
init_settings_text <- function(object) {
  settings <- object$settings
  parts <- c(if (!identical(settings$speed_var, "free")) sprintf(", speed_var = \"%s\"", settings$speed_var),
             if (!identical(settings$itemtype, "2pl")) sprintf(", itemtype = \"%s\"", settings$itemtype))
  paste(parts, collapse = "")
}

# init = of gibbs(): "random", "ecm" (any case) or an ecm() fit; quantile models start at random
check_init <- function(init, quantile) {
  quantile_text <- function(what) {
    sprintf("%s starts the chains of a mean model; quantile models (quantile = %s) start from random values: drop init",
            what, paste(format(quantile), collapse = ", "))
  }
  if (inherits(init, "rtirt")) {
    if (!has_em(init)) {
      stop("init must be an ecm() fit, \"random\" or \"ecm\"", call. = FALSE)
    }
    if (!is.null(quantile)) {
      stop(quantile_text("init = <ecm() fit>"), call. = FALSE)
    }
    return(init)
  }
  if (is.character(init) && length(init) == 1 && !is.na(init) && tolower(init) %in% c("random", "ecm")) {
    init <- tolower(init)
    if (init == "ecm" && !is.null(quantile)) {
      stop(quantile_text("init = \"ecm\""), call. = FALSE)
    }
    return(init)
  }
  shown <- paste(deparse(init, width.cutoff = 60L)[1], collapse = "")
  stop(sprintf("init must be \"random\", \"ecm\" or an ecm() fit of the same data and model (init = fit_ecm), not %s",
               shown), call. = FALSE)
}

#' @rdname gibbs
#' @export
gibbs.input_data <- function(data, model = NULL, itemtype = "2pl", speed_var = c("free", "fixed"), quantile = NULL,
                             cov = NULL, intercept = FALSE, priors = NULL, n_iter = 5000, n_chain = 4, seed = NULL,
                             init = "random", control = gibbs_control(), ...) {
  check_dots("gibbs", match.call(), match.call(expand.dots = FALSE)$...)
  control <- as_control(control, "gibbs")
  speed_var <- match.arg(speed_var)
  object <- build_model(data, model, itemtype, speed_var, cov, quantile, intercept)
  init <- check_init(init, quantile)
  if (inherits(init, "rtirt")) {
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
