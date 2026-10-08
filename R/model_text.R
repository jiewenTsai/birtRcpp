# Equations of declared models (block_model()) in Unicode or LaTeX, with the layout of equations()
# for the built-in models: persons j (rows), items i (columns), y_ji. Also the CAVI algorithm written
# the mean-field way: the model p(), the factorized q(), and each factor of q.

# a symbol with its index: theta -> theta_j, a -> a_i, Y -> Y_ji; sigma2 -> sigma_i^2
bt_sym <- function(n, F, M, index = TRUE) {
  d <- M$pars[[n]]$dim %||% (if (n %in% names(M$ddim)) M$ddim[[n]])
  ix <- if (!index || is.null(d)) "" else switch(d, person = "j", item = "i", NK = "ji", "")
  sq <- grepl("^[a-z]+2$", n) && sub("2$", "", n) %in% F$greek
  base <- if (n %in% F$greek) F$g(n) else if (sq) F$g(sub("2$", "", n)) else F$v(n)
  out <- if (nzchar(ix)) F$s(base, ix) else base
  if (sq) F$sq(out) else out
}

# drop + 0, - 0, 0 *, 1 * and redundant parentheses from an expression (a prior mean of 0 gives r = y - 0;
# a parameter set to 0 gives the rest of an expression), and move signs out of products
bt_simplify <- function(e) {
  if (!is.call(e)) return(e)
  for (k in seq_along(e)[-1]) if (!is.null(e[[k]])) e[[k]] <- bt_simplify(e[[k]])
  op <- as.character(e[[1]]); z <- function(x, v) is.numeric(x) && length(x) == 1 && x == v
  neg <- function(x) is.call(x) && identical(x[[1]], as.name("-")) && length(x) == 2
  if (op == "(" && (is.name(e[[2]]) || is.numeric(e[[2]]))) return(e[[2]])
  if (op == "-" && length(e) == 2) { if (z(e[[2]], 0)) return(0); if (neg(e[[2]])) return(e[[2]][[2]])
    if (is.numeric(e[[2]])) return(-e[[2]]) }
  if (op %in% c("+", "-") && length(e) == 3 && z(e[[3]], 0)) return(e[[2]])
  if (op == "+" && length(e) == 3 && z(e[[2]], 0)) return(e[[3]])
  if (op == "-" && length(e) == 3 && z(e[[2]], 0)) return(bt_simplify(call("-", e[[3]])))
  if (op == "*" && (z(e[[2]], 0) || z(e[[3]], 0))) return(0)
  if (op == "/" && z(e[[2]], 0)) return(0)
  if (op == "*" && z(e[[2]], 1)) return(e[[3]])
  if (op == "*" && z(e[[3]], 1)) return(e[[2]])
  if (op == "*" && z(e[[3]], -1)) return(bt_simplify(call("-", e[[2]])))
  if (op == "*" && neg(e[[2]])) return(bt_simplify(call("-", call("*", e[[2]][[2]], e[[3]]))))
  if (op == "*" && neg(e[[3]])) return(bt_simplify(call("-", call("*", e[[2]], e[[3]][[2]]))))
  if (op == "-" && length(e) == 3 && neg(e[[3]])) return(call("+", e[[2]], e[[3]][[2]]))
  if (op == "-" && length(e) == 3 && bt_lead_neg(e[[3]])) return(call("+", e[[2]], bt_negate(e[[3]])))
  if (op == "+" && length(e) == 3 && neg(e[[3]])) return(call("-", e[[2]], e[[3]][[2]]))
  if (op == "(" && neg(e[[2]]) && (is.name(e[[2]][[2]]) || is.numeric(e[[2]][[2]]))) return(e[[2]])
  e
}
# a sum whose first term is negative (x - (-a - b) is written x + a + b), and its negation
bt_lead_neg <- function(e) {
  while (is.call(e) && ((as.character(e[[1]]) %in% c("+", "-") && length(e) == 3) || identical(e[[1]], as.name("(")))) e <- e[[2]]
  is.call(e) && identical(e[[1]], as.name("-")) && length(e) == 2
}
bt_negate <- function(e) {
  if (is.call(e) && identical(e[[1]], as.name("("))) return(bt_negate(e[[2]]))
  if (is.call(e) && identical(e[[1]], as.name("+")) && length(e) == 3) return(call("-", bt_negate(e[[2]]), e[[3]]))
  if (is.call(e) && identical(e[[1]], as.name("-")) && length(e) == 3) return(call("+", bt_negate(e[[2]]), e[[3]]))
  if (is.call(e) && identical(e[[1]], as.name("-")) && length(e) == 2) return(e[[2]])
  call("-", e)
}
# an expression with parameter v set to 0
bt_zero <- function(e, v) bt_simplify(do.call(substitute, list(e, stats::setNames(list(0), v))))

# an R expression of the model (indices added)
bt_expr <- function(e, F, M) {
  strip <- function(x) if (is.call(x) && identical(x[[1]], as.name("("))) x[[2]] else x
  rec <- function(e) {
    if (is.numeric(e)) { v <- format(e, trim = TRUE); return(if (F$u) gsub("-", F$minus, v) else v) }
    if (is.name(e)) return(bt_sym(as.character(e), F, M))
    if (!is.call(e)) return(paste(deparse(e), collapse = ""))
    op <- as.character(e[[1]]); a <- function(k) rec(e[[k]])
    sum_ <- function(x) is.call(x) && as.character(x[[1]]) %in% c("+", "-")            # needs parentheses after - and in products
    p <- function(k) if (sum_(e[[k]])) F$paren(a(k)) else a(k)
    switch(op,
      "(" = F$paren(a(2)),
      "+" = if (length(e) == 2) a(2) else sprintf("%s + %s", a(2), a(3)),
      "-" = if (length(e) == 2) paste0(F$minus, p(2)) else sprintf("%s %s %s", a(2), F$minus, if (sum_(e[[3]]) && length(e[[3]]) == 3) F$paren(a(3)) else a(3)),
      "*" = { l <- p(2); r <- p(3)
              if (is.call(e[[3]]) && (identical(e[[3]][[1]], as.name("(")) || sum_(e[[3]]))) paste0(l, r) else paste0(l, F$dot, r) },
      "/" = F$div(rec(strip(e[[2]])), rec(strip(e[[3]]))),
      "^" = if (F$u && identical(e[[3]], 2)) paste0(a(2), "\u00b2") else if (F$u) sprintf("%s^%s", a(2), a(3)) else sprintf("{%s}^{%s}", a(2), a(3)),
      sqrt = if (F$u) sprintf("\u221a(%s)", rec(strip(e[[2]]))) else sprintf("\\sqrt{%s}", rec(strip(e[[2]]))),
      exp = , log = F$fun(op, rec(strip(e[[2]]))),
      paste(deparse(e), collapse = ""))
  }
  rec(bt_simplify(e))
}
# an expression to be squared or multiplied: in parentheses unless it is one symbol or already in them
bt_wrap <- function(e, F, M) {
  e <- bt_simplify(e)
  atom <- is.name(e) || is.numeric(e) || (is.call(e) && as.character(e[[1]]) %in% c("(", "exp", "log", "sqrt"))
  if (atom) bt_expr(e, F, M) else F$paren(bt_expr(e, F, M))
}

# the variance of a normal factor, as an expression (sd = sqrt(x) -> x; a number -> its square)
bt_var <- function(sd) if (is.call(sd) && identical(sd[[1]], as.name("sqrt"))) sd[[2]] else if (is.numeric(sd)) sd^2 else call("^", call("(", sd), 2)

# the pieces of notation that equations() of the built-in models does not need
bt_fmt <- function(F) {
  u <- F$u
  sb <- function(ix) gsub("\\[|\\]", "", F$s("", ix))
  given <- function(g) if (length(g)) paste0(if (u) " | " else " \\mid ", paste(g, collapse = ", ")) else ""
  c(F, list(
    paren = function(x) sprintf(if (u) "(%s)" else "\\left(%s\\right)", x),
    div = function(a, b) if (u) sprintf("%s/%s", a, b) else sprintf("\\frac{%s}{%s}", a, b),
    E = function(x, minus = "") if (u) sprintf("E%s[%s]", if (minus == "x") "\u208b\u2093" else "", x)
        else sprintf("\\mathbb{E}%s\\left[%s\\right]", if (nzchar(minus)) sprintf("_{-%s}", minus) else "", x),
    sum = function(ix) if (!nzchar(ix)) "" else if (u) paste0("\u03a3", sb(ix), " ") else sprintf("\\sum_{%s} ", ix),
    prod = function(ix) if (!nzchar(ix)) "" else if (u) paste0("\u220f", sb(ix), " ") else sprintf("\\prod_{%s} ", ix),
    rest = if (u) "\u00b7" else "\\cdot", prop = if (u) "\u221d" else "\\propto",
    half = if (u) "\u00bd" else "\\tfrac{1}{2}", log = if (u) "log" else "\\log", exp = if (u) "exp" else "\\exp",
    const = if (u) "const" else "\\text{const}", star = if (u) "*" else "^{*}",
    dist = function(n, args) sprintf(if (u) "%s(%s)" else "\\mathrm{%s}(%s)", n, paste(args, collapse = ", ")),
    p = function(x, g = character()) sprintf("p(%s%s)", x, given(g)),
    q = function(x, g = character()) sprintf("q(%s%s)", x, given(g))))
}

# numbers and truncation bounds in equations() and algorithm(): LaTeX, or Unicode with the notation F (F$u)
fmt_num <- function(x, F = NULL) { v <- format(signif(x, 3), trim = TRUE, drop0trailing = TRUE); if (isTRUE(F$u)) gsub("-", F$minus, v) else v }
fmt_trunc <- function(lo, hi, F = NULL) {
  if (!is.finite(lo) && !is.finite(hi)) return("")
  sprintf(if (isTRUE(F$u)) " T(%s, %s)" else "\\,T(%s, %s)", if (is.finite(lo)) fmt_num(lo, F) else "", if (is.finite(hi)) fmt_num(hi, F) else "")
}
# a normal on (0, Inf) is written N+ (the half flag of F$N), without T(0, )
is_half <- function(lo, hi) identical(lo, 0) && !is.finite(hi)

# one line per factor of the model: data first, then the priors
bt_model_lines <- function(M, F) {
  X <- function(e) bt_expr(e, F, M)
  one <- function(f) {
    lhs <- bt_sym(f$y, F, M)
    switch(f$dist,
      bernoulli_logit = eq_line(sprintf("%s P(%s = 1)", F$op("logit"), lhs), F$eq, X(f$eta)),
      normal = eq_line(lhs, F$sim, paste0(F$N(X(f$mean), X(bt_var(f$sd)), is_half(f$lower, f$upper)), if (!is_half(f$lower, f$upper)) fmt_trunc(f$lower, f$upper, F))),
      half_t = eq_line(if (grepl("2$", f$y) && sub("2$", "", f$y) %in% F$greek) F$s(F$g(sub("2$", "", f$y)), switch(M$pars[[f$y]]$dim, item = "i", person = "j", ""))
                       else (if (F$u) sprintf("\u221a%s", lhs) else sprintf("\\sqrt{%s}", lhs)),
                       F$sim, F$dist("half-t", c(fmt_num(f$hyper[1], F), fmt_num(f$hyper[2], F))), F$text("half-t prior on the SD")),
      inv_gamma = eq_line(lhs, F$sim, F$dist("IG", c(fmt_num(f$hyper[1], F), fmt_num(f$hyper[2], F)))))
  }
  list(data = lapply(Filter(function(f) f$observed, M$factors), one), priors = lapply(Filter(function(f) !f$observed, M$factors), one))
}

#' @rdname block_model
#' @param format `"unicode"` (default) or `"latex"` (an `aligned` block, as [equations()] for the
#'   built-in models).
#' @export
equations.block_model <- function(object, format = c("unicode", "latex"), ...) {
  format <- match.arg(format); F <- bt_fmt(eq_fmt(format)); L <- bt_model_lines(object, F)
  eq_object(list(list(title = sprintf("Data (j = person, N = %d; i = item, K = %d)", object$N, object$K), lines = L$data),
                 list(title = "Priors", lines = L$priors)), format)
}

# pieces shared by algorithm() of gibbs(), ecm() and vi() fits
bt_idx <- function(d) switch(d, person = "j", item = "i", NK = "j,i", "")
bt_over <- function(f, d) {                                              # the index summed over: the shape of f minus that of the parameter
  fi <- switch(f$shape, person = "j", item = "i", NK = c("j", "i"), character()); vi <- switch(d, person = "j", item = "i", character())
  paste(setdiff(fi, vi), collapse = ",")
}
bt_grp <- function(items, dims, F) {                                     # one product per index: data cells, persons, items, the rest
  out <- character()
  for (d in c("NK", "person", "item", "scalar")) if (any(dims == d)) out <- c(out, paste0(F$prod(bt_idx(d)), paste(items[dims == d], collapse = " ")))
  paste(out, collapse = " ")
}
# omega, its expectation (bar = TRUE) and kappa = y - 1/2 of a logistic part (with the data name when there are several)
bt_logi_sym <- function(f, M, F, base) {
  ix <- if (sum(vapply(M$factors, function(g) g$dist == "bernoulli_logit", TRUE)) > 1) paste0(f$y, ",ji") else "ji"
  F$s(base, ix)
}
bt_om <- function(f, M, F, bar = FALSE) bt_logi_sym(f, M, F, if (!bar) F$g("omega") else if (F$u) "\u03c9\u0304" else "\\bar{\\omega}")
bt_kappa <- function(f, M, F) bt_logi_sym(f, M, F, F$g("kappa"))
bt_aux <- function(v, M, F) {                                            # the Huang-Wand auxiliary variable of a half-t prior
  hw <- names(M$pars)[vapply(M$pars, function(P) P$dist == "half_t", TRUE)]; d <- M$pars[[v]]$dim
  same <- sum(vapply(hw, function(w) M$pars[[w]]$dim == d, TRUE)) > 1
  ix <- switch(d, person = "j", item = "i", ""); base <- if (same) F$v(paste0("u.", v)) else "u"; if (nzchar(ix)) F$s(base, ix) else base
}
bt_p_section <- function(M, F) {                                         # the model p: a product of its parts, then each part
  S <- function(n) bt_sym(n, F, M); L <- bt_model_lines(M, F)
  pp <- vapply(M$factors, function(f) F$p(S(f$y), vapply(setdiff(f$vars, f$y), S, "")), "")
  c(list(eq_line("p", F$eq, bt_grp(pp, vapply(M$factors, function(f) if (f$observed) "NK" else M$pars[[f$y]]$dim, ""), F))), L$data, L$priors)
}

#' @rdname vi.block_model
#' @param format `"unicode"` (default) or `"latex"`.
#' @export
algorithm.block_model_vi <- function(object, format = c("unicode", "latex"), ...) {
  format <- match.arg(format); F <- bt_fmt(eq_fmt(format))
  E <- bm_setup(object$model, object$steps, object$ecm$latent, object$ecm$nodes, "map", NULL, TRUE, object$ecm$quadrature %||% "adaptive", from_data = FALSE); M <- E$M
  X <- function(e) bt_expr(e, F, M); W <- function(e) bt_wrap(e, F, M); S <- function(n) bt_sym(n, F, M)
  over <- bt_over; grp <- function(items, dims) bt_grp(items, dims, F)
  logi <- Filter(function(f) f$dist == "bernoulli_logit", M$factors)
  om <- function(f) bt_om(f, M, F)
  grid <- vapply(E$grid, S, "")
  hw <- E$fixed[vapply(E$fixed, function(v) M$pars[[v]]$dist == "half_t", TRUE)]
  aux <- function(v) bt_aux(v, M, F)
  sec_p <- bt_p_section(M, F)                                            # p(): the model
  # q(): the mean-field factorization and the update rule
  qf <- c(F$q(paste(vapply(E$latent, S, ""), collapse = ", ")), vapply(logi, function(f) F$q(om(f), grid), ""),
          vapply(E$fixed, function(v) F$q(S(v)), ""), vapply(hw, function(v) F$q(aux(v)), ""))
  qd <- c("person", rep("NK", length(logi)), vapply(E$fixed, function(v) M$pars[[v]]$dim, ""), vapply(hw, function(v) M$pars[[v]]$dim, ""))
  sec_q <- list(eq_line("q", F$eq, grp(qf, qd), F$text("mean field: one factor per block")),
                eq_line(sprintf("%s q%s(x)", F$log, F$star), F$eq, sprintf("%s + %s", F$E(sprintf("%s %s", F$log, F$p(F$rest)), "x"), F$const),
                        F$text("each block x in turn, expectation over the others (Bishop, 2006, Eq. 10.9)")))
  # the factors of q, in the order of one sweep
  qs <- list(eq_line(F$q(paste(grid, collapse = ", ")), F$prop, sprintf("%s %s", F$exp, F$E(sprintf("%s %s", F$log, F$p(F$rest)), paste(grid, collapse = ","))),
                     F$text(sprintf("free form on %s %s; expectation over the other blocks", paste(E$Q, collapse = " x "),
                                    if (E$quadrature == "ba81") "common node(s), the same for every person (Bock-Aitkin)" else "adaptive node(s) per person"))))
  if (length(E$collapse)) qs[[length(qs) + 1]] <- eq_line(F$text(paste(vapply(E$collapse, function(id) M$factors[[id]]$y, ""), collapse = ", ")), "",
                                                          F$text("summed over persons at each node (artificial data): the updates work on nodes x items"), "")
  for (u in E$ana) qs[[length(qs) + 1]] <- eq_line(F$q(S(u), grid), F$eq, F$N(F$div(F$s("h", "j"), F$s("p", "j")), F$div("1", F$s("p", "j"))),
                                                     F$text(if (E$sigma) sprintf("normal: exact given the grid variables; the updates use each person's mean and covariance of (%s) through %d points (exact: these parts are quadratic in them)",
                                                                                 paste(E$latent, collapse = ", "), 2L * length(E$latent)) else "normal: exact given the grid variables"))
  for (f in logi) {
    qs[[length(qs) + 1]] <- eq_line(F$q(om(f), grid), F$eq, F$dist("PG", c("1", F$s("c", "ji"))), F$text("Durante & Rigon (2019)"))
    qs[[length(qs) + 1]] <- eq_line(F$sq(F$s("c", "ji")), F$eq, F$E(F$sq(W(f$eta))), "")
  }
  sp <- if (F$u) " " else "\\,"
  for (st in E$steps) for (v in st$vars) {
    P <- M$pars[[v]]; d <- P$dim; ii <- switch(d, person = "j", item = "i", ""); sub_ <- function(x) if (nzchar(ii)) F$s(x, ii) else x
    if (P$dist == "normal") {
      terms <- vapply(factors_of(M, v), function(f) {
        if (!f$observed && f$y == v) return(if (is.numeric(f$sd)) (if (f$sd == 1) "1" else F$div("1", fmt_num(f$sd^2, F))) else F$E(F$div("1", X(bt_var(f$sd)))))
        c1 <- f$d1[[v]]; if (is.call(c1) && identical(c1[[1]], as.name("-")) && length(c1) == 2) c1 <- c1[[2]]
        cc <- if (is.numeric(c1) && abs(c1) == 1) "" else F$sq(W(c1))
        s <- F$sum(over(f, d))
        if (f$dist == "bernoulli_logit") return(paste0(s, F$E(paste(c(om(f), if (nzchar(cc)) cc), collapse = sp))))
        w <- if (is.numeric(f$sd)) (if (f$sd == 1) "" else F$div("1", fmt_num(f$sd^2, F))) else F$E(F$div("1", X(bt_var(f$sd))))
        body <- paste(c(if (nzchar(w)) w, if (nzchar(cc)) F$E(cc)), collapse = sp)
        paste0(s, if (nzchar(body)) body else "1")
      }, "")
      qs[[length(qs) + 1]] <- eq_line(F$q(S(v)), F$eq, paste0(F$N(F$div(sub_("m"), sub_("p")), F$div("1", sub_("p")), is_half(P$lower, P$upper)),
                                                              if (!is_half(P$lower, P$upper)) fmt_trunc(P$lower, P$upper, F)),
                                      F$text("the Gibbs normal part with expected coefficients; m: its linear terms"))
      qs[[length(qs) + 1]] <- eq_line(sub_("p"), F$eq, paste(terms, collapse = " + "), "")
    } else {
      res <- vapply(Filter(function(f) !(f$y == v && !f$observed), factors_of(M, v)), function(f) paste0(F$sum(over(f, d)), F$E(F$sq(W(f$r)))), "")
      a0 <- if (P$dist == "inv_gamma") fmt_num(P$hyper[1], F) else F$div(fmt_num(P$hyper[1], F), "2")
      b0 <- if (P$dist == "inv_gamma") fmt_num(P$hyper[2], F) else paste0(fmt_num(P$hyper[1], F), F$dot, F$E(F$div("1", aux(v))))
      qs[[length(qs) + 1]] <- eq_line(F$q(S(v)), F$eq, F$dist("IG", c(sprintf("%s + %s", a0, F$div(sub_("n"), "2")), sprintf("%s + %s%s", b0, F$half, F$paren(paste(res, collapse = " + "))))),
                                      F$text("n: the number of terms in the sum"))
      if (P$dist == "half_t") qs[[length(qs) + 1]] <- eq_line(F$q(aux(v)), F$eq, F$dist("IG", c(fmt_num((P$hyper[1] + 1) / 2, F),
                                                                 sprintf("%s%s%s + %s", fmt_num(P$hyper[1], F), F$dot, F$E(F$div("1", S(v))), fmt_num(1 / P$hyper[2]^2, F)))),
                                                               F$text(sprintf("auxiliary of the half-t(%s, %s) prior (Huang & Wand, 2013)", fmt_num(P$hyper[1], F), fmt_num(P$hyper[2], F))))
    }
  }
  sec_s <- list(eq_line(F$text("one sweep"), "", F$text("the factors of q above, in this order; each update raises the ELBO"), ""),
                eq_line(F$text("rounds"), "", F$text("the grid is recentred between rounds, as in ecm(); stop when the ELBO changes by less than tol"), ""))
  eq_object(list(list(title = "Model p", lines = sec_p), list(title = "Approximation q", lines = sec_q),
                 list(title = "Factors of q", lines = qs), list(title = "CAVI", lines = sec_s)), format)
}

# ---- algorithm() of gibbs() and ecm() fits: each step as its full conditional or conditional mode -----------

# The normal part of a parameter v with a normal prior: v | . ~ N(m/p, 1/p) (Gibbs) or v <- m/p (ECM), with
# p and m summed over the parts facs of the model that use v. For a part with c = d(expression)/dv and the
# rest of the expression at v = 0 (eta0, r0): logistic p += omega c^2, m += c (kappa - omega eta0); normal
# p += c^2 / variance, m -= c r0 / variance; the prior of v p += 1 / sd^2, m += mean / sd^2. For ECM
# (lat: the parts with latent variables) omega is its expectation and the terms are averaged over the E-step.
bt_normal_part <- function(M, v, F, facs, lat = function(f) FALSE, bar = FALSE) {
  X <- function(e) bt_expr(e, F, M); W <- function(e) bt_wrap(e, F, M)
  Ws <- function(e) if (is.call(e) && as.character(e[[1]]) %in% c("+", "-")) F$paren(X(e)) else X(e)    # parentheses only around sums
  d <- M$pars[[v]]$dim; sp <- if (F$u) " " else "\\,"
  dv <- function(num, den, prod = FALSE) if (identical(den, "1")) num else if (prod) F$div(num, den) else bt_frac(num, den, F)
  vt <- function(f) { x <- bt_var(f$sd); if (is.numeric(x)) fmt_num(x, F) else X(x) }
  join <- function(t) { t <- t[nzchar(t)]; if (!length(t)) return("0"); out <- paste(t, collapse = " + ")
    gsub(paste0("\\+ ", F$minus), paste0(F$minus, " "), out) }
  pt <- character(); mt <- character()
  for (f in facs) {
    if (!f$observed && f$y == v) {                                     # the prior of v
      pt <- c(pt, dv("1", vt(f)))
      mu <- bt_simplify(f$mean); if (!(is.numeric(mu) && mu == 0)) mt <- c(mt, dv(X(mu), vt(f)))
      next
    }
    c1 <- f$d1[[v]]; neg <- is.call(c1) && identical(c1[[1]], as.name("-")) && length(c1) == 2; if (neg) c1 <- c1[[2]]
    c1 <- bt_simplify(c1); one <- is.numeric(c1) && c1 == 1
    E <- function(x) if (lat(f)) F$E(x) else x
    s <- F$sum(bt_over(f, d))
    if (f$dist == "bernoulli_logit") {
      om <- bt_om(f, M, F, bar); ka <- bt_kappa(f, M, F); e0 <- bt_zero(f$eta, v)
      en <- is.call(e0) && identical(e0[[1]], as.name("-")) && length(e0) == 2      # kappa - omega (-x) is written kappa + omega x
      inner <- if (is.numeric(e0) && e0 == 0) ka else sprintf("%s %s %s%s%s", ka, if (en) "+" else F$minus, om, sp, Ws(if (en) e0[[2]] else e0))
      pt <- c(pt, paste0(s, E(paste(c(om, if (!one) F$sq(W(c1))), collapse = sp))))
      mt <- c(mt, paste0(if (neg) F$minus, s, E(if (one) inner else if (identical(inner, ka)) paste0(W(c1), sp, ka) else paste0(W(c1), F$paren(inner)))))
    } else {
      r0 <- bt_zero(f$r, v)
      pt <- c(pt, paste0(s, dv(if (one) "1" else E(F$sq(W(c1))), vt(f))))
      if (!(is.numeric(r0) && r0 == 0))
        mt <- c(mt, paste0(if (!neg) F$minus, s, dv(E(if (one) X(r0) else paste0(W(c1), Ws(r0))), vt(f), prod = !one || lat(f))))
    }
  }
  list(p = join(pt), m = join(mt))
}

# a fraction; in Unicode a compound numerator or denominator (a sum, or a number times a symbol) gets parentheses
bt_frac <- function(num, den, F) {
  if (!F$u) return(F$div(num, den))
  cp <- function(x) (grepl("[ +\u2212-]", x) || grepl("^[0-9.]+[^0-9.]", x)) && !grepl("^\\(.*\\)$", x) && !grepl("^E\\[.*\\]$", x)
  F$div(if (cp(num)) F$paren(num) else num, if (cp(den)) F$paren(den) else den)
}

# the sum of squared residuals of the parts in which the variance v is the variance (sd = sqrt(v))
bt_sq_resid <- function(M, v, F, facs, lat = function(f) FALSE) {
  d <- M$pars[[v]]$dim
  t <- vapply(Filter(function(f) !(f$y == v && !f$observed), facs), function(f) {
    r <- F$sq(bt_wrap(bt_simplify(f$r), F, M)); paste0(F$sum(bt_over(f, d)), if (lat(f)) F$E(r) else r) }, "")
  paste(t, collapse = " + ")
}

# the lines of the steps: Gibbs (E = NULL) or the conditional maximizations of ECM (E from bm_setup())
bt_step_lines <- function(M, steps, F, E = NULL) {
  ecm <- !is.null(E); S <- function(n) bt_sym(n, F, M); X <- function(e) bt_expr(e, F, M)
  bar_ <- if (F$u) " | " else " \\mid "; tx <- function(x) if (F$u) x else sprintf("\\text{%s}", x)
  sub_ <- function(x, d) { ix <- switch(d, person = "j", item = "i", ""); if (nzchar(ix)) F$s(x, ix) else x }
  lat <- function(f) ecm && any(f$vars %in% E$latent)
  inq <- function(f) !ecm || bm_in_q(E, f)
  out <- list(); stale <- rep(TRUE, length(M$factors)); why <- rep("drawn once per iteration", length(M$factors)); moved <- vector("list", length(M$factors))
  L <- function(...) out[[length(out) + 1]] <<- eq_line(...)
  for (st in steps) {
    if (st$type == "gibbs") for (v in st$vars) {
      P <- M$pars[[v]]; d <- P$dim; facs <- Filter(inq, factors_of(M, v))
      if (length(st$cut)) facs <- Filter(function(f) !(f$observed && f$y %in% st$cut), facs)
      if (!ecm) for (f in pg_factors(M, v)) if (stale[f$id] && !(f$y %in% st$cut)) {
        L(sprintf("%s%s%s", bt_om(f, M, F), bar_, F$rest), F$sim, F$dist("PG", c("1", X(f$eta))), F$text(sprintf("automatic: %s", why[f$id])))
        L(bt_kappa(f, M, F), F$eq, sprintf("%s %s %s", S(f$y), F$minus, F$half), "")
        stale[f$id] <- FALSE; moved[f$id] <- list(NULL)
      }
      if (P$dist == "normal") {
        np <- bt_normal_part(M, v, F, facs, lat, bar = ecm); half <- is_half(P$lower, P$upper)
        dist <- paste0(F$N(F$div(sub_("m", d), sub_("p", d)), F$div("1", sub_("p", d)), half), if (!half) fmt_trunc(P$lower, P$upper, F))
        cut <- if (length(st$cut)) sprintf(" (cut: %s left out)", paste(st$cut, collapse = ", ")) else ""
        if (ecm) L(S(v), if (F$u) "\u2190" else "\\leftarrow",
                   sprintf("%s%s", F$div(sub_("m", d), sub_("p", d)), if (is.finite(P$lower) || is.finite(P$upper)) tx(sprintf(", clipped to [%s, %s]", fmt_num(P$lower), fmt_num(P$upper))) else ""),
                   F$text("conditional mode: the mode of the Gibbs step, its parts averaged over the E-step"))
        else L(sprintf("%s%s%s", S(v), bar_, F$rest), F$sim, dist, F$text(sprintf("Gibbs: normal full conditional from %d part(s) of the model%s", length(facs), cut)))
        L(sub_("p", d), F$eq, np$p, ""); L(sub_("m", d), F$eq, np$m, "")
      } else {
        rs <- bt_sq_resid(M, v, F, facs, lat); n <- sub_("n", d); nu <- P$hyper[1]; A <- P$hyper[2]
        note_n <- sprintf("%s: the number of observed terms in the sum", if (F$u) n else sprintf("n%s", switch(d, person = "_j", item = "_i", "")))
        if (!ecm && P$dist == "inv_gamma")
          L(sprintf("%s%s%s", S(v), bar_, F$rest), F$sim, F$dist("IG", c(sprintf("%s + %s", fmt_num(P$hyper[1], F), F$div(n, "2")), sprintf("%s + %s%s", fmt_num(P$hyper[2], F), F$half, F$paren(rs)))),
            F$text(sprintf("Gibbs: conjugate; %s", note_n)))
        if (!ecm && P$dist == "half_t") {
          u <- bt_aux(v, M, F)
          L(sprintf("%s%s%s", u, bar_, S(v)), F$sim, F$dist("IG", c(fmt_num((nu + 1) / 2, F), sprintf("%s + %s", F$div(fmt_num(nu, F), S(v)), fmt_num(1 / A^2, F)))),
            F$text(sprintf("auxiliary of the half-t(%s, %s) prior of the SD (Huang & Wand, 2013), drawn from the current %s", fmt_num(nu, F), fmt_num(A, F), v)))
          L(sprintf("%s%s%s, %s", S(v), bar_, u, F$rest), F$sim, F$dist("IG", c(sprintf("%s + %s", F$div(fmt_num(nu, F), "2"), F$div(n, "2")), sprintf("%s + %s%s", F$div(fmt_num(nu, F), u), F$half, F$paren(rs)))),
            F$text(sprintf("Gibbs: conjugate; %s", note_n)))
        }
        if (ecm) {
          ar <- if (F$u) "\u2190" else "\\leftarrow"
          rsum <- F$paren(rs); half_rs <- paste0(F$half, rsum)
          if (E$prior == "ml") L(S(v), ar, bt_frac(rsum, n, F), F$text(sprintf("conditional mode of the variance (maximum likelihood); %s", note_n)))
          else if (P$dist == "inv_gamma") L(S(v), ar, bt_frac(sprintf("%s + %s", fmt_num(P$hyper[2], F), half_rs), sprintf("%s + %s + 1", fmt_num(P$hyper[1], F), F$div(n, "2")), F),
                                            F$text(sprintf("conditional mode of the variance (inverse-gamma); %s", note_n)))
          else {
            rt <- if (F$u) sprintf("\u221a%s", S(v)) else sprintf("\\sqrt{%s}", S(v))
            tnu <- if (F$u) paste0("t", eq_sub(fmt_num(nu, F))) else sprintf("t_{%s}", fmt_num(nu, F))
            L(S(v), ar, sprintf("%s %s%s %s %s %s %s %s(%s)", F$op("argmax"), F$minus, bt_frac(sprintf("%s + 1", n), "2", F), F$fun("log", S(v)), F$minus,
                                bt_frac(half_rs, S(v), F), "+", paste(F$op("log"), tnu), if (A == 1) rt else F$div(rt, fmt_num(A, F))),
              F$text(sprintf("conditional mode with the half-t(%s, %s) prior of the SD (1-D maximization); %s", fmt_num(nu, F), fmt_num(A, F), note_n)))
          }
        }
      }
    }
    if (st$type == "mala") {
      tr <- vapply(st$vars, function(v) { P <- M$pars[[v]]; if (is.finite(P$lower)) F$fun("log", if (P$lower == 0) S(v) else sprintf("%s %s %s", S(v), F$minus, fmt_num(P$lower, F))) else S(v) }, "")
      jac <- st$vars[vapply(st$vars, function(v) is.finite(M$pars[[v]]$lower), TRUE)]
      lhs <- if (length(tr) > 1) F$paren(paste(tr, collapse = ", ")) else tr
      if (ecm) L(lhs, if (F$u) "\u2190" else "\\leftarrow", sprintf("%s Q(%s)", tx("gradient ascent on"), paste(vapply(st$vars, S, ""), collapse = ", ")),
                 F$text(sprintf("generalized EM (the step_mala() counterpart): a few steps with step halving on the expected log %s; log scale for bounded parameters (no Jacobian)",
                                if (E$prior == "map") "posterior" else "likelihood")))
      else {
        L(lhs, if (F$u) "\u2190" else "\\leftarrow", sprintf("%s %s %s", tx("MALA on"), F$op("log"), F$p(paste(vapply(st$vars, S, ""), collapse = ", "), F$rest)),
          F$text(sprintf("one target per %s; exact likelihood (no omega)%s; step sizes adapt during burn-in", if (st$dim == "scalar") "model" else st$dim,
                         if (length(jac)) sprintf("; Jacobian of the log scale added for %s", paste(jac, collapse = ", ")) else "")))
        for (f in M$factors) if (f$dist == "bernoulli_logit" && any(st$vars %in% f$vars)) {
          moved[[f$id]] <- union(if (stale[f$id] && !is.null(moved[[f$id]])) moved[[f$id]], intersect(st$vars, f$vars))
          stale[f$id] <- TRUE; why[f$id] <- sprintf("redrawn because MALA moved %s", paste(moved[[f$id]], collapse = ", "))
        }
      }
    }
    if (st$type == "shift" && !ecm) {
      new <- vapply(st$vars, function(v) X(call("+", as.name(v), call("*", as.name(".c"), st$along[[v]]))), "")
      L(F$paren(paste(vapply(st$vars, S, ""), collapse = ", ")), if (F$u) "\u2190" else "\\leftarrow", sprintf("%s, c%s%s %s N", F$paren(paste(gsub(".c", "c", new, fixed = TRUE), collapse = ", ")), bar_, F$rest, F$sim),
        F$text("shift (generalized Gibbs; Liu & Sabatti, 2000): the logistic parts do not change"))
    }
  }
  out
}
eq_sub <- function(x) { m <- c(`0` = "\u2080", `1` = "\u2081", `2` = "\u2082", `3` = "\u2083", `4` = "\u2084", `5` = "\u2085", `6` = "\u2086",
                                `7` = "\u2087", `8` = "\u2088", `9` = "\u2089", `.` = ".")
  ch <- strsplit(x, "")[[1]]; if (all(ch %in% names(m))) paste(m[ch], collapse = "") else paste0("_", x) }

#' @rdname gibbs.block_model
#' @param object A `block_model_fit`, or a `block_model` (then give `steps`).
#' @param format `"unicode"` (default) or `"latex"` (an `aligned` block).
#' @param ... Unused.
#' @export
algorithm.block_model <- function(object, steps = NULL, format = c("unicode", "latex"), ...) {
  format <- match.arg(format); F <- bt_fmt(eq_fmt(format)); M <- object; steps <- prepare_steps(M, steps)
  eq_object(list(list(title = "Model p", lines = bt_p_section(M, F)),
                 list(title = "Gibbs: one iteration, the full conditionals in this order", lines = bt_step_lines(M, steps, F)),
                 list(title = "Notation", lines = list(eq_line(F$rest, "", F$text("all other parameters and the data (their current values)"), ""),
                                                       eq_line("m, p", "", F$text("the linear term and the precision of a normal full conditional, summed over the parts of the model"), ""))))
            , format)
}

#' @rdname gibbs.block_model
#' @export
algorithm.block_model_fit <- function(object, format = c("unicode", "latex"), ...) algorithm.block_model(object$model, steps = object$steps, format = format)

#' @rdname ecm.block_model
#' @param object A fit of `ecm()`.
#' @param format `"unicode"` (default) or `"latex"` (an `aligned` block).
#' @param ... Unused.
#' @export
algorithm.block_model_ecm <- function(object, format = c("unicode", "latex"), ...) {
  format <- match.arg(format); F <- bt_fmt(eq_fmt(format)); method <- object$ecm$method %||% "ecm"
  E <- bm_setup(object$model, object$steps, object$ecm$latent, object$ecm$nodes, object$ecm$prior, NULL, isTRUE(object$ecm$analytic %||% TRUE),
                object$ecm$quadrature %||% "adaptive", from_data = FALSE)
  M <- E$M; S <- function(n) bt_sym(n, F, M); X <- function(e) bt_expr(e, F, M)
  lat <- vapply(E$latent, S, ""); latj <- paste(lat, collapse = ", "); phi <- F$g("phi"); ar <- if (F$u) "\u2190" else "\\leftarrow"
  yj <- F$s("y", "j"); lp <- if (E$prior == "map") sprintf(" + %s %s", F$op("log"), F$p(phi)) else ""
  others <- paste(vapply(E$fixed, S, ""), collapse = ", ")
  obj <- list(eq_line(sprintf("%s(%s)", if (F$u) "\u2113" else "\\ell", phi), F$eq,
                      sprintf("%s%s %s %s d%s%s", F$sum("j"), F$op("log"), if (F$u) "\u222b" else "\\int", F$p(paste(c(yj, latj), collapse = ", "), phi),
                              paste(lat, collapse = if (F$u) " d" else "\\,d"), lp),
                      F$text(sprintf("objective: log marginal %s; %s = (%s)", if (E$prior == "map") "posterior" else "likelihood", "phi", paste(E$fixed, collapse = ", ")))))
  es <- list(eq_line(sprintf("Q(%s%s%s%s)", phi, if (F$u) " | " else " \\mid ", phi, F$prime), F$eq,
                     paste0(F$sum("j"), F$E(sprintf("%s %s", F$op("log"), F$p(paste(c(yj, latj), collapse = ", "), phi))), lp),
                     F$text(sprintf("expected complete-data log %s; E over p(%s | y_j, phi') at the current phi'", if (E$prior == "map") "posterior" else "likelihood",
                                    paste(E$latent, collapse = ", ")))))
  for (l in seq_along(E$grid)) {
    v <- E$grid[l]; zk <- F$s("z", "k")
    es[[length(es) + 1]] <- if (E$quadrature == "ba81")
      eq_line(S(v), F$eq, sprintf("%s + %s%s", "m", "s", paste0(F$dot, zk)), F$text(sprintf("Bock-Aitkin: the same %d Gauss-Hermite node(s) z_k for every person; m, s: the marginal posterior mean and SD (recentred between rounds)", E$Q[l])))
    else eq_line(S(v), F$eq, sprintf("%s + %s%s", F$s("m", "j"), F$s("s", "j"), paste0(F$dot, zk)),
                 F$text(sprintf("adaptive Gauss-Hermite: %d node(s) z_k per person (%d on the grid); m_j, s_j: posterior mean and SD, fixed within a round and recentred between rounds", E$Q[l], E$ng)))
  }
  if (length(E$collapse)) es[[length(es) + 1]] <- eq_line(F$text(paste(vapply(E$collapse, function(id) M$factors[[id]]$y, ""), collapse = ", ")), "",
                                                          F$text("summed over persons at each node (artificial data, Bock-Aitkin): the M-step works on nodes x items"), "")
  for (u in E$ana) es[[length(es) + 1]] <- eq_line(sprintf("%s%s%s, %s", S(u), if (F$u) " | " else " \\mid ", paste(vapply(E$grid, S, ""), collapse = ", "), yj), F$sim,
                                                   F$N(F$div(F$s("h", "j"), F$s("p", "j")), F$div("1", F$s("p", "j"))),
                                                   F$text(sprintf("integrated exactly at each node: normal and linear in %s (completing the square)%s", u,
                                                                  if (E$sigma) "; p and h are sums over items per person" else "")))
  if (E$sigma) es[[length(es) + 1]] <- eq_line(F$E(sprintf("g(%s)", latj)), "", F$text(sprintf("in the M-step for the parts of %s: each person's posterior mean and covariance of (%s) through %d points mean -/+ sqrt(%d) x column of chol(covariance), exact for these quadratic parts (Julier & Uhlmann, 2004)",
                                                            paste(E$ana, collapse = ", "), paste(E$latent, collapse = ", "), 2L * length(E$latent), length(E$latent))), "")
  if (method == "ecm") for (f in M$factors) if (f$dist == "bernoulli_logit") {
    om <- bt_om(f, M, F, TRUE); et <- F$s(F$g("eta"), "ji")
    es[[length(es) + 1]] <- eq_line(om, F$eq, bt_frac(sprintf("%s(%s)", F$op("tanh"), F$div(et, "2")), sprintf("2%s", et), F),
                                    F$text(sprintf("E[omega | eta] at each node, eta = %s: the Polya-Gamma minorizer (Polson, Scott & Windle, 2013)", gsub("\\s+", " ", deparse1(f$eta)))))
    es[[length(es) + 1]] <- eq_line(bt_kappa(f, M, F), F$eq, sprintf("%s %s %s", S(f$y), F$minus, F$half), "")
  }
  secs <- list(list(title = "Model p", lines = bt_p_section(M, F)), list(title = "Objective", lines = obj), list(title = "E-step", lines = es))
  if (method == "ecm") {
    secs[[4]] <- list(title = "M-step: conditional maximizations, in this order", lines = bt_step_lines(M, E$steps, F, E))
    secs[[5]] <- list(title = "Then", lines = list(
      eq_line(F$text("E"), "", F$text(sprintf("expectation over the E-step (the nodes and their posterior weights) of the parts that use %s", paste(E$latent, collapse = ", "))), ""),
      eq_line(F$text("acceleration"), "", F$text("SQUAREM (run_em()); an iteration never lowers the objective"), ""),
      eq_line(F$text("standard errors"), "", F$text("from the Hessian of the objective by Louis' formula"), "")))
  } else secs[[4]] <- list(title = "Maximization", lines = list(
    eq_line(phi, ar, sprintf("%s %s(%s)", F$op("argmax"), if (F$u) "\u2113" else "\\ell", phi),
            F$text(sprintf("nlminb() within the bounds, all parameters together; gradient by Fisher's identity (the expected complete-data score), %s (Glas & van der Linden, 2010, Sec. 3.1)",
                           if (method == "newton") "Hessian by Louis' formula (Newton)" else "quasi-Newton Hessian"))),
    eq_line(F$text("standard errors"), "", F$text("from the Hessian of the objective by Louis' formula"), "")))
  eq_object(secs, format)
}
