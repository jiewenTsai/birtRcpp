# Multilinear form of the expressions of a declared model: a sum of monomials coef * x1 * x2 * ...
# Linear predictors such as a * (theta - b) or lambda - zeta - rho * theta are multilinear; their
# value, and the derivative with respect to any parameter, are computed in C++ (.ml_eval) from the
# person and item vectors directly. Expressions that are not multilinear (exp(), log(), sqrt() of a
# parameter, division by a parameter) keep the R evaluation.

ml_mono <- function(coef, syms = character()) list(coef = coef, syms = syms)

ml_compile <- function(e) {
  if (is.numeric(e) && length(e) == 1) return(list(ml_mono(as.numeric(e))))
  if (is.name(e)) return(list(ml_mono(1, as.character(e))))
  if (!is.call(e)) return(NULL)
  op <- as.character(e[[1]])
  scale <- function(a, k) if (is.null(a)) NULL else lapply(a, function(m) { m$coef <- m$coef * k; m })
  prod2 <- function(a, b) if (is.null(a) || is.null(b)) NULL else
    unlist(lapply(a, function(x) lapply(b, function(y) ml_mono(x$coef * y$coef, c(x$syms, y$syms)))), recursive = FALSE)
  if (op == "(") return(ml_compile(e[[2]]))
  if (op %in% c("+", "-") && length(e) == 2) return(scale(ml_compile(e[[2]]), if (op == "-") -1 else 1))
  if (op %in% c("+", "-")) {
    a <- ml_compile(e[[2]]); b <- ml_compile(e[[3]])
    if (is.null(a) || is.null(b)) return(NULL)
    return(ml_simplify(c(a, scale(b, if (op == "-") -1 else 1))))
  }
  if (op == "*") return(ml_simplify(prod2(ml_compile(e[[2]]), ml_compile(e[[3]]))))
  if (op == "/" && is.numeric(e[[3]])) return(scale(ml_compile(e[[2]]), 1 / e[[3]]))
  if (op == "^" && is.numeric(e[[3]]) && e[[3]] %in% 0:4) {
    a <- ml_compile(e[[2]]); if (is.null(a)) return(NULL)
    out <- list(ml_mono(1)); for (k in seq_len(e[[3]])) out <- prod2(out, a)
    return(ml_simplify(out))
  }
  NULL
}

ml_simplify <- function(ir) {
  if (is.null(ir)) return(NULL)
  key <- vapply(ir, function(m) paste(sort(m$syms), collapse = "*"), "")
  out <- lapply(split(seq_along(ir), factor(key, levels = unique(key))), function(i)
    ml_mono(sum(vapply(ir[i], `[[`, 0, "coef")), sort(ir[[i[1]]]$syms)))
  out <- Filter(function(m) m$coef != 0, unname(out))
  if (!length(out)) list(ml_mono(0)) else out
}

# derivative with respect to v (multilinear again)
ml_deriv <- function(ir, v) {
  out <- list()
  for (m in ir) { k <- sum(m$syms == v)
    if (k) out[[length(out) + 1]] <- ml_mono(m$coef * k, m$syms[-match(v, m$syms)]) }
  ml_simplify(if (length(out)) out else list(ml_mono(0)))
}

# a compiled expression: the monomials plus what .ml_eval() needs (symbols, their dimension codes, ids)
ml_pack <- function(ir, M) {
  if (is.null(ir)) return(NULL)
  syms <- unique(unlist(lapply(ir, `[[`, "syms")))
  code <- c(scalar = 0L, person = 1L, item = 2L, NK = 3L)
  list(ir = ir, coef = vapply(ir, `[[`, 0, "coef"), syms = syms, ids = lapply(ir, function(m) match(m$syms, syms) - 1L),
       dims = unname(vapply(syms, function(n) code[[if (n %in% names(M$pars)) M$pars[[n]]$dim else M$ddim[[n]]]], 0L)),
       is_par = syms %in% names(M$pars))
}

# value of a compiled expression for a factor of the given shape
ml_eval <- function(pk, s, M, shape) {
  vals <- lapply(seq_along(pk$syms), function(k) { x <- if (pk$is_par[k]) s[[pk$syms[k]]] else M$data0[[pk$syms[k]]]
    if (is.double(x)) x else as.numeric(x) })                       # the C++ side reads the doubles in place
  nr <- if (shape %in% c("NK", "person")) M$N else 1L; nc <- if (shape %in% c("NK", "item")) M$K else 1L
  out <- .ml_eval(pk$coef, pk$ids, vals, pk$dims, nr, nc)
  switch(shape, NK = out, person = out[, 1], item = out[1, ], scalar = out[1, 1])
}

# the compiled parts of a factor: r = y - mean (normal) or eta (logistic), its derivatives, and the variance
ml_factor <- function(f, M) {
  ir <- list()
  main <- if (f$dist == "normal") f$r else if (f$dist == "bernoulli_logit") f$eta else NULL
  if (is.null(main)) return(NULL)
  key <- if (f$dist == "normal") "r" else "eta"
  ir[[key]] <- ml_compile(main)
  if (f$dist == "normal") {
    ir$mean <- ml_compile(f$mean)
    sd <- f$sd
    if (is.call(sd) && identical(sd[[1]], as.name("sqrt"))) ir$sd2 <- ml_compile(sd[[2]])
    else { c1 <- ml_compile(sd); if (!is.null(c1)) ir$sd2 <- ml_simplify(unlist(lapply(c1, function(x) lapply(c1, function(y) ml_mono(x$coef * y$coef, c(x$syms, y$syms)))), recursive = FALSE)) }
  }
  out <- lapply(ir, ml_pack, M = M)
  if (!is.null(ir[[key]])) out$d1 <- lapply(stats::setNames(f$vars, f$vars), function(v) ml_pack(ml_deriv(ir[[key]], v), M))
  if (!is.null(ir$sd2)) out$dsd2 <- lapply(stats::setNames(f$sd_vars, f$sd_vars), function(v) ml_pack(ml_deriv(ir$sd2, v), M))
  out
}

# derivative of the log density of a normal factor with respect to v, elementwise (r and the variance s2 given)
gnorm <- function(f, s, M, v, r, s2) {
  g <- -(r / s2) * fval(f, s, M, "d1", v)
  if (v %in% f$sd_vars) g <- g + if (!is.null(f$ir$dsd2[[v]])) (r^2 / (2 * s2^2) - 1 / (2 * s2)) * ml_eval(f$ir$dsd2[[v]], s, M, f$shape)
                                 else { sdv <- sqrt(s2); (r^2 / sdv^3 - 1 / sdv) * fval(f, s, M, "dsd", v) }
  g
}

# value of a part of factor f: "eta", "r", "mean", "sd", "d1" (with v), "dsd" (with v) or "y"
fval <- function(f, s, M, key, v = NULL) {
  ir <- f$ir
  if (!is.null(ir)) {
    if (key == "d1" && !is.null(ir$d1[[v]])) return(ml_eval(ir$d1[[v]], s, M, f$shape))
    if (key %in% c("eta", "r", "mean") && !is.null(ir[[key]])) return(ml_eval(ir[[key]], s, M, f$shape))
    if (key == "sd" && !is.null(ir$sd2)) return(sqrt(ml_eval(ir$sd2, s, M, f$shape)))
    if (key == "sd2" && !is.null(ir$sd2)) return(ml_eval(ir$sd2, s, M, f$shape))
  }
  if (key == "y" && f$observed) return(M$data0[[f$y]])
  ev <- fenv(f, s, M)
  expr <- switch(key, eta = f$eta, r = f$r, mean = f$mean, sd = f$sd, d1 = f$d1[[v]], dsd = f$dsd[[v]], y = as.name(f$y),
                 sd2 = call("^", f$sd, 2))
  full_of(eval(expr, ev), f, M)
}

# products and parts of compiled expressions (for expectations under a factorized q)
ml_mul <- function(a, b) ml_simplify(unlist(lapply(a, function(x) lapply(b, function(y) ml_mono(x$coef * y$coef, c(x$syms, y$syms)))), recursive = FALSE))
ml_drop <- function(ir, v) { out <- Filter(function(m) !(v %in% m$syms), ir); if (length(out)) out else list(ml_mono(0)) }

# a compiled expression whose symbols are powers: x (mean) and x^2 (second moment)
ml_pack_mom <- function(ir, M) {
  ir2 <- lapply(ir, function(m) { tb <- table(m$syms); ml_mono(m$coef, if (length(tb)) ifelse(tb == 1, names(tb), paste0(names(tb), "^", tb)) else character()) })
  pk <- ml_pack(ir2, M_mom(M, ir2))
  pk$base <- sub("\\^.*$", "", pk$syms); pk$power <- rep(1L, length(pk$syms)); hp <- grepl("\\^", pk$syms); pk$power[hp] <- as.integer(sub("^.*\\^", "", pk$syms[hp]))
  pk$is_par <- pk$base %in% names(M$pars)
  if (any(pk$power > 2)) stop("vi(): an expression has a parameter to a power above 2 after squaring; it is not supported", call. = FALSE)
  pk
}
M_mom <- function(M, ir2) {                                            # dimension lookup that knows the x^2 names
  syms <- unique(unlist(lapply(ir2, `[[`, "syms"))); base <- sub("\\^.*$", "", syms)
  for (k in seq_along(syms)) if (syms[k] != base[k]) { if (base[k] %in% names(M$pars)) M$pars[[syms[k]]] <- M$pars[[base[k]]] else M$ddim[[syms[k]]] <- M$ddim[[base[k]]] }
  M
}
# E over q of a moment-packed expression: mom$m (means) and mom$m2 (second moments) per parameter
ml_eval_mom <- function(pk, mom, M, shape) {
  vals <- lapply(seq_along(pk$syms), function(k) {
    b <- pk$base[k]; x <- if (pk$is_par[k]) (if (pk$power[k] == 1) mom$m[[b]] else mom$m2[[b]]) else M$data0[[b]]^pk$power[k]
    if (is.double(x)) x else as.numeric(x) })
  nr <- if (shape %in% c("NK", "person")) M$N else 1L; nc <- if (shape %in% c("NK", "item")) M$K else 1L
  out <- .ml_eval(pk$coef, pk$ids, vals, pk$dims, nr, nc)
  switch(shape, NK = out, person = out[, 1], item = out[1, ], scalar = out[1, 1])
}

# a compiled expression with its values, in the form .ml_part() reads
ml_bind <- function(pk, s, M) {
  if (is.null(pk)) return(list())
  vals <- lapply(seq_along(pk$syms), function(k) { x <- if (pk$is_par[k]) s[[pk$syms[k]]] else M$data0[[pk$syms[k]]]
    if (is.double(x)) x else as.numeric(x) })
  list(pk$coef, pk$ids, vals, pk$dims)
}

# the normal part of parameter v from factor f in C++ (NULL when the factor is not compiled)
# mode: "normal", "gibbs" (omega given), "em" (E[omega]); roww: weights per row (stacked nodes)
ml_part_f <- function(f, v, s, M, mode, om = NULL, roww = NULL, red_dim = M$pars[[v]]$dim, keep_om = FALSE) {
  ir <- f$ir; key <- if (f$dist == "normal") "r" else "eta"
  if (is.null(ir) || is.null(ir[[key]]) || is.null(ir$d1[[v]]) || (f$dist == "normal" && is.null(ir$sd2))) return(NULL)
  s0 <- s; s0[[v]] <- 0 * s[[v]]
  nr <- if (f$shape %in% c("NK", "person")) M$N else 1L; nc <- if (f$shape %in% c("NK", "item")) M$K else 1L
  red <- if (red_dim == "scalar") 0L else if (f$shape == "NK") (if (red_dim == "person") 1L else 2L) else (if (f$shape == "person") 1L else 2L)
  mask <- if (is.null(f$mask)) NULL else { x <- f$mask; storage.mode(x) <- "double"; x }
  third <- switch(mode, normal = ml_bind(ir$sd2, s, M), em = ml_bind(ir$eta, s, M), list())
  md <- switch(mode, normal = 0L, gibbs = 1L, em = 2L)
  kap <- if (f$dist == "bernoulli_logit") f$kappa else NULL
  if (!is.null(om)) storage.mode(om) <- "double"
  .ml_part(ml_bind(ir$d1[[v]], s, M), ml_bind(ir[[key]], s0, M), third, om, mask, kap, roww, md, red, nr, nc, keep_om)
}
