# equations(): the model written out as equations (Unicode or LaTeX) -----------------------
#
# Same layout as equations() in birt: sections of lines (lhs, relation, rhs, note), rendered
# by a format object F. Persons j, items i (y_ij: item i, person j); the person variables are ability and speed.

#' Model equations
#'
#' Writes the model as equations, with the parameter names of [summary()] and the priors of
#' [gibbs()]. The default is Unicode text for the console; `format = "latex"` gives an
#' `aligned` block for a paper or an appendix (e.g. `writeLines(as.character(equations(fit, "latex")), "model.tex")`).
#' Works before estimation too: give the data and the model as in [ecm()] and [gibbs()].
#' @param object A fit of [ecm()] or [gibbs()], data from [input_data()] (with `model`), or a model
#'   object of the deprecated constructors ([models]) or an `rtirt_qset`.
#' @param format `"unicode"` or `"latex"`.
#' @param priors Include the priors (omitted for an [ecm()] fit, unless `prior = "map"`).
#' @param x An `rtirt_equations` object.
#' @param ... Unused.
#' @return An `rtirt_equations` object (printed); `as.character()` gives the lines.
#' @examples
#' cond <- set_cond(n_subj = 100, n_item = 5)
#' dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' equations(dat, model = "cross")
#' equations(dat, model = "cross", quantile = 0.25, format = "latex")
#' @export
equations <- function(object, ...) {
  UseMethod("equations")
}

# LaTeX special characters in names and text
tex_esc <- function(x) {
  gsub("([&%#$_{}])", "\\\\\\1", x)
}

eq_fmt <- function(format) {
  u <- format == "unicode"
  greek <- c(alpha = "\u03b1", beta = "\u03b2", gamma = "\u03b3", delta = "\u03b4", epsilon = "\u03b5", zeta = "\u03b6",
             eta = "\u03b7", theta = "\u03b8", kappa = "\u03ba", lambda = "\u03bb", nu = "\u03bd", xi = "\u03be",
             pi = "\u03c0", rho = "\u03c1", sigma = "\u03c3", tau = "\u03c4", psi = "\u03c8", mu = "\u03bc", phi = "\u03c6",
             omega = "\u03c9", iota = "\u03b9", chi = "\u03c7", upsilon = "\u03c5")
  plus_sup <- if (u) "\u207a" else "^{+}"

  list(u = u, greek = names(greek),
    g = function(n) if (u) greek[[n]] else paste0("\\", n),
    v = function(n) eq_variable(n, u),
    s = function(x, s) if (u) paste0(x, unicode_subscript(s)) else sprintf("%s_{%s}", x, s),
    sim = if (u) "\u223c" else "\\sim",
    eq = "=",
    dot = if (u) "\u00b7" else "\\,",
    minus = if (u) "\u2212" else "-",
    sq = function(x) if (x == "1") "1" else if (u) paste0(x, "\u00b2") else paste0(x, "^2"),
    plus_sup = plus_sup,
    prime = if (u) "\u2032" else "'",
    inset = if (u) " \u2208 " else " \\in ",
    N = function(m, v, plus = FALSE) {
      template <- if (u) "N%s(%s, %s)" else "\\mathcal{N}%s(%s, %s)"
      sprintf(template, if (plus) plus_sup else "", m, v)
    },
    fun = function(f, x) if (u) sprintf("%s(%s)", f, x) else sprintf("\\%s(%s)", f, x),
    op = function(f) if (u) f else if (f %in% c("log", "exp")) paste0("\\", f) else sprintf("\\operatorname{%s}", f),
    text = function(t) if (u) t else sprintf("\\text{%s}", tex_esc(t)),
    set = function(items) eq_set(items, u))
}

# Unicode subscript of s, or [s] when a character has no subscript form
unicode_subscript <- function(s) {
  subscripts <- c(i = "\u1d62", j = "\u2c7c", k = "\u2096", p = "\u209a", `0` = "\u2080", `1` = "\u2081", `2` = "\u2082",
                  `3` = "\u2083", `4` = "\u2084", `5` = "\u2085", `6` = "\u2086", `7` = "\u2087", `8` = "\u2088", `9` = "\u2089",
                  `-` = "\u208b", `+` = "\u208a", `,` = ",")
  chars <- strsplit(as.character(s), "")[[1]]
  if (all(chars %in% names(subscripts))) paste(subscripts[chars], collapse = "") else paste0("[", s, "]")
}

# a variable name: as is in Unicode and for one-letter names in LaTeX, otherwise \mathrm{}
eq_variable <- function(n, u) {
  if (u) {
    return(n)
  }
  if (nchar(n) == 1 && !grepl("[&%#$_{}]", n)) {
    return(n)
  }
  sprintf("\\mathrm{%s}", tex_esc(n))
}

# a set of item names, shortened to the first two, an ellipsis and the last with more than 4
eq_set <- function(items, u) {
  shown <- items
  if (length(items) > 4) {
    shown <- c(items[1:2], if (u) "\u2026" else "\\ldots", items[length(items)])
  }
  if (u) {
    return(sprintf("{%s}", paste(shown, collapse = ", ")))
  }
  latex_item <- function(x) if (x == "\\ldots") x else sprintf("\\mathrm{%s}", tex_esc(x))
  sprintf("\\{%s\\}", paste(vapply(shown, latex_item, ""), collapse = ", "))
}

eq_line <- function(lhs, rel, rhs, note = "") {
  list(lhs = lhs, rel = rel, rhs = rhs, note = note)
}

eq_object <- function(sections, format) {
  non_empty <- Filter(function(s) length(s$lines) > 0, sections)
  structure(list(sections = non_empty, format = format), class = "rtirt_equations")
}

#' @rdname equations
#' @export
as.character.rtirt_equations <- function(x, ...) {
  if (x$format == "unicode") {
    return(eq_unicode_lines(x$sections))
  }
  eq_latex_lines(x$sections)
}

# Unicode: a title per section, the lines aligned at the relation
eq_unicode_lines <- function(sections) {
  out <- character()
  for (section in sections) {
    out <- c(out, if (length(out)) "", paste0(section$title, ":"))
    lines <- section$lines
    lhs_width <- max(nchar(vapply(lines, `[[`, "", "lhs"), type = "width"))
    for (line in lines) {
      pad <- strrep(" ", lhs_width - nchar(line$lhs, type = "width"))
      note <- if (nzchar(line$note)) paste0("    ", line$note) else ""
      text <- sprintf("  %s%s %s %s%s", pad, line$lhs, line$rel, line$rhs, note)
      out <- c(out, sub("\\s+$", "", text))
    }
  }
  out
}

# LaTeX: one aligned environment, section titles as text lines
eq_latex_lines <- function(sections) {
  latex_line <- function(line) {
    rel <- if (nzchar(line$rel)) line$rel else "\\quad"
    note <- if (nzchar(line$note)) sprintf(" \\quad %s", line$note) else ""
    sprintf("  %s &%s %s%s \\\\", line$lhs, rel, line$rhs, note)
  }
  latex_section <- function(section) {
    c(sprintf("  & \\text{%s:} \\\\", section$title), vapply(section$lines, latex_line, ""))
  }
  body <- unlist(lapply(sections, latex_section))
  body[length(body)] <- sub(" \\\\\\\\$", "", body[length(body)])
  c("\\begin{aligned}", body, "\\end{aligned}")
}

#' @rdname equations
#' @export
print.rtirt_equations <- function(x, ...) {
  cat(as.character(x), sep = "\n")
  invisible(x)
}

#' @rdname equations
#' @export
equations.rtirt <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format)
  F <- eq_fmt(format)
  ctx <- eq_context(object, F)
  sections <- list(list(title = "Measurement", lines = eq_measurement(ctx, F)),
                   list(title = "Person variables", lines = eq_persons(ctx, F)))
  ecm_only <- is.null(object$post) && !is.null(object$ecm)
  map <- isTRUE(object$ecm$map)
  if (priors && (!ecm_only || map)) {
    title <- if (ecm_only) "Priors (item parameters; prior = \"map\")" else "Priors"
    sections[[3]] <- list(title = title, lines = eq_priors(object, ctx, F, ecm_only))
  }
  if (ecm_only) {
    method <- if (map) "posterior mode in the item parameters, maximum likelihood in the others" else "maximum marginal likelihood; persons integrated out"
    estimation <- eq_line(F$text("ECM"), "", F$text(method))
    sections[[length(sections) + 1]] <- list(title = "Estimation", lines = list(estimation))
  }
  eq_object(sections, format)
}

# what the equations of a fit depend on: engine, settings and the symbols of the person variables
eq_context <- function(object, F) {
  info <- .model_info[[class(object)[1]]]
  settings <- object$settings
  cov_names <- colnames(object$data$X) %||% character()
  list(info = info, eng = info$engine, settings = settings, qr = isTRUE(object$qr), q = object$q,
       fixed = identical(settings$speed_var, "fixed"), onepl = identical(settings$itemtype, "1pl"),
       A = F$s(F$v("ability"), "j"), Z = F$s(F$v("speed"), "j"), items = object$data$items,
       has_cov = length(cov_names) > 0)
}

eq_ald <- function(ctx, F, scale) {
  template <- if (F$u) "ALD[q = %s](0, %s)" else "\\mathrm{ALD}_{q=%s}(0, %s)"
  sprintf(template, format(ctx$q), scale)
}

# x_j' b, or NULL without covariates
eq_xb <- function(ctx, F, b) {
  if (ctx$has_cov) paste0(F$s(F$v("x"), "j"), F$prime, b) else NULL
}

eq_measurement <- function(ctx, F) {
  A <- ctx$A
  discrimination <- if (ctx$onepl) "" else F$s("a", "i")
  note <- sprintf("i%s%s", F$inset, F$set(ctx$items))
  lines <- list(eq_line(sprintf("%s P(%s = 1)", F$op("logit"), F$s("y", "ij")), F$eq,
                        sprintf("%s(%s %s %s)", discrimination, A, F$minus, F$s("b", "i")), note))
  if (!ctx$info$rt) {
    return(lines)
  }
  cross <- if (ctx$eng == "cross") c(F$minus, paste0(F$s(F$g("rho"), "i"), F$dot, A))
  linear <- paste(c(F$s(F$g("lambda"), "i"), F$minus, ctx$Z, cross), collapse = " ")
  epsilon <- F$s(F$g("epsilon"), "ij")
  sigma <- F$s(F$g("sigma"), "i")
  ald <- ctx$qr && ctx$eng == "cross"
  residual <- if (ald) eq_ald(ctx, F, sigma) else F$N("0", F$sq(sigma))
  residual_note <- if (ald) F$text(sprintf("quantile regression: the %s-quantile of log t is linear", format(ctx$q))) else ""
  c(lines, list(eq_line(paste(F$op("log"), F$s("t", "ij")), F$eq, paste(linear, "+", epsilon), note),
                eq_line(epsilon, F$sim, residual, residual_note)))
}

eq_persons <- function(ctx, F) {
  lines <- switch(ctx$eng,
    mlirt = eq_persons_mlirt(ctx, F),
    rtirt = eq_persons_rtirt(ctx, F),
    latent = eq_persons_latent(ctx, F),
    cross = eq_persons_cross(ctx, F))
  if (isTRUE(ctx$settings$intercept)) {
    lines[[length(lines) + 1]] <- eq_line(F$text("x"), "", F$text("includes an intercept column"), "")
  }
  lines
}

eq_persons_mlirt <- function(ctx, F) {
  e1 <- F$s("e", "1j")
  list(eq_line(ctx$A, F$eq, paste(c(eq_xb(ctx, F, F$g("beta")), e1), collapse = " + "), ""),
       eq_line(e1, F$sim, F$N("0", "1"), ""))
}

eq_persons_rtirt <- function(ctx, F) {
  e1 <- F$s("e", "1j")
  e2 <- F$s("e", "2j")
  speed <- paste(c(eq_xb(ctx, F, F$s(F$g("beta"), "2")), paste0("c", F$dot, e1), e2), collapse = " + ")
  e2_var <- if (ctx$fixed) paste0("1 ", F$minus, " c", if (F$u) "\u00b2" else "^2") else "v"
  correlation <- if (ctx$fixed) "c" else if (F$u) "c / \u221a(c\u00b2 + v)" else "c / \\sqrt{c^2 + v}"
  correlation_note <- if (ctx$fixed) {
    F$text("Var(speed | x) = 1")
  } else if (F$u) {
    "var_speed = c\u00b2 + v"
  } else {
    "\\text{var\\_speed} = c^2 + v"
  }
  list(eq_line(ctx$A, F$eq, paste(c(eq_xb(ctx, F, F$s(F$g("beta"), "1")), e1), collapse = " + "), ""),
       eq_line(ctx$Z, F$eq, speed, ""),
       eq_line(e1, F$sim, F$N("0", "1"), ""),
       eq_line(e2, F$sim, F$N("0", e2_var), ""),
       eq_line(sprintf("%s(%s, %s)", F$op("Corr"), ctx$A, ctx$Z), F$eq, correlation, correlation_note))
}

eq_persons_latent <- function(ctx, F) {
  u <- F$s("u", "j")
  speed <- paste(c(eq_xb(ctx, F, F$g("beta")), paste0(F$s("b", "ability"), F$dot, ctx$A), u), collapse = " + ")
  residual <- if (ctx$qr) eq_ald(ctx, F, "s") else F$N("0", "s")
  residual_note <- if (ctx$qr) F$text(sprintf("quantile regression of speed (q = %s)", format(ctx$q))) else F$text("s = var_speed")
  list(eq_line(ctx$A, F$sim, F$N("0", "1"), ""),
       eq_line(ctx$Z, F$eq, speed, ""),
       eq_line(u, F$sim, residual, residual_note))
}

eq_persons_cross <- function(ctx, F) {
  speed_note <- if (ctx$fixed) "" else F$text("s = var_speed")
  list(eq_line(ctx$A, F$sim, F$N("0", "1"), ""),
       eq_line(ctx$Z, F$sim, F$N("0", if (ctx$fixed) "1" else "s"), speed_note))
}

# the priors of the fit (item parameters only for ecm(prior = "map"))
eq_priors <- function(object, ctx, F, ecm_only) {
  P <- fit_priors(object)
  pf <- eq_prior_fmt(F)
  lines <- list(
    if (!ctx$onepl) eq_line(F$s("a", "i"), F$sim, F$N(pf$num(P$a[1]), pf$var(P$a[2]), TRUE), ""),
    eq_line(F$s("b", "i"), F$sim, paste0(F$N(pf$num(P$b[["mean"]]), pf$var(P$b[["sd"]])), pf$trunc(P$b[["lower"]], P$b[["upper"]])), ""))
  if (ctx$info$rt) {
    lines <- c(lines, eq_priors_rt(P, object$data$log_t, F, pf))
  }
  if (!ecm_only) {
    if (ctx$has_cov) {
      lines[[length(lines) + 1]] <- eq_line(F$g("beta"), F$sim, F$N("0", pf$var(P$beta[["sd"]], "I")), F$text("each coefficient"))
    }
    lines <- c(lines, eq_priors_structural(ctx, P, F, pf))
  }
  Filter(Negate(is.null), lines)
}

# number, variance, truncation and variance-prior formats of the prior lines
eq_prior_fmt <- function(F) {
  num <- function(x) {
    value <- format(signif(x, 3), trim = TRUE, drop0trailing = TRUE)
    if (F$u) gsub("-", "\u2212", value) else value
  }
  variance <- function(sd, one = "1") {
    if (sd == 1) {
      return(one)
    }
    if (one == "1") F$sq(num(sd)) else paste0(F$sq(num(sd)), one)
  }
  name <- function(n) if (F$u) n else sprintf("\\mathrm{%s}", n)
  truncation <- function(lo, hi) {
    if (!is.finite(lo) && !is.finite(hi)) {
      return("")
    }
    template <- if (F$u) " T(%s, %s)" else "\\,T(%s, %s)"
    sprintf(template, if (is.finite(lo)) num(lo) else "", if (is.finite(hi)) num(hi) else "")
  }
  variance_prior <- function(x, family) {
    distribution <- if (family == "half-t") name("half-t") else name("IG")
    sprintf("%s(%s, %s)", distribution, num(x[[1]]), num(x[[2]]))
  }
  # half-t on the SD, IG on the variance
  variance_lhs <- function(v2, v1, family) if (family == "half-t") v1 else v2
  sq_root <- function(x) if (F$u) paste0("\u221a", x) else sprintf("\\sqrt{%s}", x)
  list(num = num, var = variance, trunc = truncation, vprior = variance_prior, vlhs = variance_lhs, sqrt = sq_root)
}

eq_priors_rt <- function(P, L, F, pf) {
  lambda <- lambda_prior(P, L)
  lambda_i <- F$s(F$g("lambda"), "i")
  lambda_line <- if (is.null(P$lambda)) {
    eq_line(lambda_i, F$sim, F$N(format(round(lambda[1], 2)), F$sq(format(round(lambda[2], 2))), TRUE),
            F$text("centred at the mean and SD of log t"))
  } else {
    eq_line(lambda_i, F$sim, paste0(F$N(pf$num(lambda[1]), pf$var(lambda[2])), pf$trunc(lambda[3], Inf)), "")
  }
  sigma <- F$s(F$g("sigma"), "i")
  sigma_lhs <- pf$vlhs(F$sq(sigma), sigma, P$sigma2t_family)
  list(lambda_line, eq_line(sigma_lhs, F$sim, pf$vprior(P$sigma2t, P$sigma2t_family), ""))
}

eq_priors_structural <- function(ctx, P, F, pf) {
  speed_var_line <- function(name) {
    lhs <- pf$vlhs(name, pf$sqrt(name), P$var_speed_family)
    eq_line(lhs, F$sim, pf$vprior(P$var_speed, P$var_speed_family), "")
  }
  if (ctx$eng == "rtirt") {
    if (ctx$fixed) {
      uniform <- if (F$u) "U(\u22121, 1)" else "\\mathcal{U}(-1, 1)"
      return(list(eq_line("c", F$sim, uniform, F$text("c = Corr(ability, speed); Metropolis"))))
    }
    return(list(eq_line("c", F$sim, F$N("0", pf$var(P$cov[["sd"]])), ""), speed_var_line("v")))
  }
  if (ctx$eng == "latent") {
    return(list(eq_line(F$s("b", "ability"), F$sim, F$N("0", pf$var(P$beta[["sd"]])), ""), speed_var_line("s")))
  }
  if (ctx$eng == "cross") {
    rho_line <- eq_line(F$s(F$g("rho"), "i"), F$sim, F$N("0", pf$var(P$rho[["sd"]])), "")
    return(c(list(rho_line), if (!ctx$fixed) list(speed_var_line("s"))))
  }
  list()
}

#' @rdname equations
#' @export
equations.rtirt_qset <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format)
  F <- eq_fmt(format)
  eq <- equations(object$models[[1]], format, priors)
  levels <- eq_line(F$v("q"), trimws(F$inset), F$set(format(object$quantile)), F$text("one model per quantile"))
  eq$sections <- c(eq$sections, list(list(title = "Quantiles", lines = list(levels))))
  eq
}
