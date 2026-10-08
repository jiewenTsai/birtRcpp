# equations(): the model written out as equations (Unicode or LaTeX) -----------------------
#
# Same layout as equations() in birt: sections of lines (lhs, relation, rhs, note), rendered
# by a format object F. Persons j, items i (y_ij: item i, person j); the person variables are ability and speed.

#' Model equations
#'
#' Writes the model as equations, with the parameter names of [summary()] and the priors of
#' [gibbs()]. The default is Unicode text for the console; `format = "latex"` gives an
#' `aligned` block for a paper or an appendix (e.g. `writeLines(as.character(equations(fit, "latex")), "model.tex")`).
#' Works for unfitted models too.
#' @param object A model object (see [models]), sampled or not, or an `rtirt_qset`.
#' @param format `"unicode"` or `"latex"`.
#' @param priors Include the priors (omitted for an [ecm()] fit, unless `prior = "map"`).
#' @param x An `rtirt_equations` object.
#' @param ... Unused.
#' @return An `rtirt_equations` object (printed); `as.character()` gives the lines.
#' @examples
#' cond <- set_cond(n_subj = 100, n_item = 5)
#' m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
#' equations(m)
#' equations(m, "latex")
#' @export
equations <- function(object, ...) UseMethod("equations")

# LaTeX special characters in names and text
tex_esc <- function(x) gsub("([&%#$_{}])", "\\\\\\1", x)

eq_fmt <- function(format) {
  u <- format == "unicode"
  gr <- c(alpha = "\u03b1", beta = "\u03b2", gamma = "\u03b3", delta = "\u03b4", epsilon = "\u03b5", zeta = "\u03b6",
          eta = "\u03b7", theta = "\u03b8", kappa = "\u03ba", lambda = "\u03bb", nu = "\u03bd", xi = "\u03be",
          pi = "\u03c0", rho = "\u03c1", sigma = "\u03c3", tau = "\u03c4", psi = "\u03c8", mu = "\u03bc", phi = "\u03c6",
          omega = "\u03c9", iota = "\u03b9", chi = "\u03c7", upsilon = "\u03c5")
  usub <- c(i = "\u1d62", j = "\u2c7c", k = "\u2096", p = "\u209a", `0` = "\u2080", `1` = "\u2081", `2` = "\u2082",
            `3` = "\u2083", `4` = "\u2084", `5` = "\u2085", `6` = "\u2086", `7` = "\u2087", `8` = "\u2088", `9` = "\u2089",
            `-` = "\u208b", `+` = "\u208a", `,` = ",")
  sub1 <- function(s) {
    ch <- strsplit(as.character(s), "")[[1]]
    if (all(ch %in% names(usub))) paste(usub[ch], collapse = "") else paste0("[", s, "]")
  }
  list(u = u, greek = names(gr),
    g = function(n) if (u) gr[[n]] else paste0("\\", n),
    v = function(n) if (u) n else if (nchar(n) == 1 && !grepl("[&%#$_{}]", n)) n else sprintf("\\mathrm{%s}", tex_esc(n)),
    s = function(x, s) if (u) paste0(x, sub1(s)) else sprintf("%s_{%s}", x, s),
    sim = if (u) "\u223c" else "\\sim", eq = "=", dot = if (u) "\u00b7" else "\\,", minus = if (u) "\u2212" else "-",
    sq = function(x) if (x == "1") "1" else if (u) paste0(x, "\u00b2") else paste0(x, "^2"),
    plus_sup = if (u) "\u207a" else "^{+}", prime = if (u) "\u2032" else "'", inset = if (u) " \u2208 " else " \\in ",
    N = function(m, v, plus = FALSE) sprintf(if (u) "N%s(%s, %s)" else "\\mathcal{N}%s(%s, %s)", if (!plus) "" else if (u) "\u207a" else "^{+}", m, v),
    fun = function(f, x) if (u) sprintf("%s(%s)", f, x) else sprintf("\\%s(%s)", f, x),
    op = function(f) if (u) f else if (f %in% c("log", "exp")) paste0("\\", f) else sprintf("\\operatorname{%s}", f),
    text = function(t) if (u) t else sprintf("\\text{%s}", tex_esc(t)),
    set = function(items) {
      it <- if (length(items) > 4) c(items[1:2], if (u) "\u2026" else "\\ldots", items[length(items)]) else items
      if (u) sprintf("{%s}", paste(it, collapse = ", ")) else sprintf("\\{%s\\}", paste(vapply(it, function(x) if (x == "\\ldots") x else sprintf("\\mathrm{%s}", tex_esc(x)), ""), collapse = ", "))
    })
}

eq_line <- function(lhs, rel, rhs, note = "") list(lhs = lhs, rel = rel, rhs = rhs, note = note)

eq_object <- function(sections, format) structure(list(sections = Filter(function(s) length(s$lines) > 0, sections), format = format),
                                                  class = "rtirt_equations")

#' @rdname equations
#' @export
as.character.rtirt_equations <- function(x, ...) {
  if (x$format == "unicode") {
    out <- character()
    for (s in x$sections) {
      out <- c(out, if (length(out)) "", paste0(s$title, ":"))
      L <- s$lines; wl <- max(nchar(vapply(L, `[[`, "", "lhs"), type = "width"))
      for (l in L) {
        pad <- strrep(" ", wl - nchar(l$lhs, type = "width"))
        out <- c(out, sub("\\s+$", "", sprintf("  %s%s %s %s%s", pad, l$lhs, l$rel, l$rhs, if (nzchar(l$note)) paste0("    ", l$note) else "")))
      }
    }
    return(out)
  }
  body <- unlist(lapply(x$sections, function(s)
    c(sprintf("  & \\text{%s:} \\\\", s$title),
      vapply(s$lines, function(l) sprintf("  %s &%s %s%s \\\\", l$lhs, if (nzchar(l$rel)) l$rel else "\\quad", l$rhs,
                                           if (nzchar(l$note)) sprintf(" \\quad %s", l$note) else ""), ""))))
  body[length(body)] <- sub(" \\\\\\\\$", "", body[length(body)])
  c("\\begin{aligned}", body, "\\end{aligned}")
}

#' @rdname equations
#' @export
print.rtirt_equations <- function(x, ...) { cat(as.character(x), sep = "\n"); invisible(x) }

#' @rdname equations
#' @export
equations.rtirt <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format); F <- eq_fmt(format)
  info <- .model_info[[class(object)[1]]]; eng <- info$engine; S <- object$settings
  qr <- isTRUE(object$qr); q <- object$q; fixed <- identical(S$speed_var, "fixed"); onepl <- identical(S$itemtype, "1pl")
  A <- F$s(F$v("ability"), "j"); Z <- F$s(F$v("speed"), "j"); items <- object$data$items
  xn <- colnames(object$data$X) %||% character()
  xb <- function(b) if (length(xn)) paste0(F$s(F$v("x"), "j"), F$prime, b) else NULL
  aj <- if (onepl) "" else F$s("a", "i")
  ald <- function(s) sprintf(if (F$u) "ALD[q = %s](0, %s)" else "\\mathrm{ALD}_{q=%s}(0, %s)", format(q), s)
  note <- sprintf("i%s%s", F$inset, F$set(items))
  meas <- list(eq_line(sprintf("%s P(%s = 1)", F$op("logit"), F$s("y", "ij")), F$eq, sprintf("%s(%s %s %s)", aj, A, F$minus, F$s("b", "i")), note))
  if (info$rt) {
    lin <- paste(c(F$s(F$g("lambda"), "i"), F$minus, Z, if (eng == "cross") c(F$minus, paste0(F$s(F$g("rho"), "i"), F$dot, A))), collapse = " ")
    meas[[2]] <- eq_line(paste(F$op("log"), F$s("t", "ij")), F$eq, paste(lin, "+", F$s(F$g("epsilon"), "ij")), note)
    meas[[3]] <- eq_line(F$s(F$g("epsilon"), "ij"), F$sim,
                         if (qr && eng == "cross") ald(F$s(F$g("sigma"), "i")) else F$N("0", F$sq(F$s(F$g("sigma"), "i"))),
                         if (qr && eng == "cross") F$text(sprintf("quantile regression: the %s-quantile of log t is linear", format(q))) else "")
  }
  e1 <- F$s("e", "1j"); e2 <- F$s("e", "2j")
  st <- switch(eng,
    mlirt = list(eq_line(A, F$eq, paste(c(xb(F$g("beta")), e1), collapse = " + "), ""), eq_line(e1, F$sim, F$N("0", "1"), "")),
    rtirt = c(list(eq_line(A, F$eq, paste(c(xb(F$s(F$g("beta"), "1")), e1), collapse = " + "), ""),
                   eq_line(Z, F$eq, paste(c(xb(F$s(F$g("beta"), "2")), paste0("c", F$dot, e1), e2), collapse = " + "), ""),
                   eq_line(e1, F$sim, F$N("0", "1"), ""),
                   eq_line(e2, F$sim, F$N("0", if (fixed) paste0("1 ", F$minus, " c", if (F$u) "\u00b2" else "^2") else "v"), ""),
                   eq_line(sprintf("%s(%s, %s)", F$op("Corr"), A, Z), F$eq, if (fixed) "c" else if (F$u) "c / \u221a(c\u00b2 + v)" else "c / \\sqrt{c^2 + v}",
                           if (fixed) F$text("Var(speed | x) = 1") else (if (F$u) "var_speed = c\u00b2 + v" else "\\text{var\\_speed} = c^2 + v")))),
    latent = list(eq_line(A, F$sim, F$N("0", "1"), ""),
                  eq_line(Z, F$eq, paste(c(xb(F$g("beta")), paste0(F$s("b", "ability"), F$dot, A), F$s("u", "j")), collapse = " + "), ""),
                  eq_line(F$s("u", "j"), F$sim, if (qr) ald("s") else F$N("0", "s"),
                          if (qr) F$text(sprintf("quantile regression of speed (q = %s)", format(q))) else F$text("s = var_speed"))),
    cross = list(eq_line(A, F$sim, F$N("0", "1"), ""), eq_line(Z, F$sim, F$N("0", if (fixed) "1" else "s"), if (fixed) "" else F$text("s = var_speed"))))
  if (isTRUE(S$intercept)) st[[length(st) + 1]] <- eq_line(F$text("x"), "", F$text("includes an intercept column"), "")
  secs <- list(list(title = "Measurement", lines = meas), list(title = "Person variables", lines = st))
  em <- is.null(object$post) && !is.null(object$ecm)
  if (priors && (!em || isTRUE(object$ecm$map))) {
    L <- object$data$log_t; P <- fit_priors(object)
    nn <- function(m, v, plus = FALSE) F$N(m, v, plus)
    num <- function(x) { v <- format(signif(x, 3), trim = TRUE, drop0trailing = TRUE); if (F$u) gsub("-", "\u2212", v) else v }
    vr <- function(sd, one = "1") if (sd == 1) one else if (one == "1") F$sq(num(sd)) else paste0(F$sq(num(sd)), one)
    fn <- function(name) if (F$u) name else sprintf("\\mathrm{%s}", name)
    tr <- function(lo, hi) if (!is.finite(lo) && !is.finite(hi)) "" else
      sprintf(if (F$u) " T(%s, %s)" else "\\,T(%s, %s)", if (is.finite(lo)) num(lo) else "", if (is.finite(hi)) num(hi) else "")
    vprior <- function(x, fam) if (fam == "half-t") sprintf("%s(%s, %s)", fn("half-t"), num(x[[1]]), num(x[[2]]))
                               else sprintf("%s(%s, %s)", fn("IG"), num(x[[1]]), num(x[[2]]))
    vlhs <- function(v2, v1, fam) if (fam == "half-t") v1 else v2          # half-t on the SD, IG on the variance
    sq_root <- function(x) if (F$u) paste0("\u221a", x) else sprintf("\\sqrt{%s}", x)
    lp <- lambda_prior(P, L)
    pr <- list(if (!onepl) eq_line(F$s("a", "i"), F$sim, nn(num(P$a[1]), vr(P$a[2]), TRUE), ""),
               eq_line(F$s("b", "i"), F$sim, paste0(nn(num(P$b[["mean"]]), vr(P$b[["sd"]])), tr(P$b[["lower"]], P$b[["upper"]])), ""))
    if (info$rt) pr <- c(pr, list(
      if (is.null(P$lambda)) eq_line(F$s(F$g("lambda"), "i"), F$sim, nn(format(round(lp[1], 2)), F$sq(format(round(lp[2], 2))), TRUE),
                                     F$text("centred at the mean and SD of log t"))
      else eq_line(F$s(F$g("lambda"), "i"), F$sim, paste0(nn(num(lp[1]), vr(lp[2])), tr(lp[3], Inf)), ""),
      eq_line(vlhs(F$sq(F$s(F$g("sigma"), "i")), F$s(F$g("sigma"), "i"), P$sigma2t_family), F$sim, vprior(P$sigma2t, P$sigma2t_family), "")))
    if (!em) {
      if (length(xn)) pr[[length(pr) + 1]] <- eq_line(F$g("beta"), F$sim, nn("0", vr(P$beta[["sd"]], "I")), F$text("each coefficient"))
      pr <- c(pr, switch(eng,
        rtirt = if (fixed) list(eq_line("c", F$sim, if (F$u) "U(\u22121, 1)" else "\\mathcal{U}(-1, 1)", F$text("c = Corr(ability, speed); Metropolis"))) else
                  list(eq_line("c", F$sim, nn("0", vr(P$cov[["sd"]])), ""), eq_line(vlhs("v", sq_root("v"), P$var_speed_family), F$sim, vprior(P$var_speed, P$var_speed_family), "")),
        latent = list(eq_line(F$s("b", "ability"), F$sim, nn("0", vr(P$beta[["sd"]])), ""), eq_line(vlhs("s", sq_root("s"), P$var_speed_family), F$sim, vprior(P$var_speed, P$var_speed_family), "")),
        cross = c(list(eq_line(F$s(F$g("rho"), "i"), F$sim, nn("0", vr(P$rho[["sd"]])), "")), if (!fixed) list(eq_line(vlhs("s", sq_root("s"), P$var_speed_family), F$sim, vprior(P$var_speed, P$var_speed_family), ""))),
        list()))
    }
    secs[[3]] <- list(title = if (em) "Priors (item parameters; prior = \"map\")" else "Priors", lines = Filter(Negate(is.null), pr))
  }
  if (em) secs[[length(secs) + 1]] <- list(title = "Estimation", lines = list(eq_line(F$text("ECM"), "",
    F$text(if (isTRUE(object$ecm$map)) "posterior mode in the item parameters, maximum likelihood in the others" else "maximum marginal likelihood; persons integrated out"))))
  eq_object(secs, format)
}

#' @rdname equations
#' @export
equations.rtirt_qset <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format); F <- eq_fmt(format)
  eq <- equations(object$models[[1]], format, priors)
  eq$sections <- c(eq$sections, list(list(title = "Quantiles", lines = list(eq_line(F$v("q"), trimws(F$inset), F$set(format(object$quantile)), F$text("one model per quantile"))))))
  eq
}
