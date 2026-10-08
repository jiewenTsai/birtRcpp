# Summary objects: tables computed once, printed in the layout of tppcm::irt_pars() and
# blavaan --------------------------------------------------------------------------------
#
# summary(fit) returns an "rtirt_summary": header lines and titled tables (data frames of
# class "rtirt_table"): fit, items and se (one row per item), parameters (one row per
# parameter, blocks a / b / lambda / sigma2t / rho / structural, with the prior of each for
# Gibbs fits), covariance and reliability. s$items etc. extract them, as.data.frame(s) gives
# the parameter table; print() shows each table under its name.

rtable <- function(df, title, cls = NULL) {
  if (is.null(df)) return(NULL)
  rownames(df) <- NULL
  structure(df, class = c(cls, "rtirt_table", "data.frame"), title = title)
}

#' Summary tables
#'
#' `summary()` of a sampled model (or one fitted by [ecm()]) returns its tables as one
#' object: `header`, `fit` (fit indices), `items` and `se` (one row per item: a, b, lambda,
#' sigma2t, rho and their posterior SDs or standard errors), `parameters` (one row per
#' parameter with its interval and, for Gibbs fits, its prior, as in blavaan), `covariance`
#' (person covariance) and `reliability`. Each element is a data frame; `as.data.frame()`
#' returns `parameters`. `print()` shows the tables under their names, so `s$items` gets the
#' one printed as `$items`.
#' @param x An `rtirt_summary`, or one of its tables.
#' @param tables Names of the tables to print (default all).
#' @param digits Decimals.
#' @param title Print the title of a table.
#' @param row.names,optional Unused.
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `as.data.frame()` the parameter table.
#' @name rtirt_summary
NULL

#' @rdname rtirt_summary
#' @export
print.rtirt_summary <- function(x, tables = NULL, digits = 3, ...) {
  cat(x$header, sep = "\n")
  nm <- setdiff(names(x), "header")
  for (n in if (is.null(tables)) nm else intersect(tables, nm)) {
    t <- x[[n]]
    if (is.null(t) || (is.data.frame(t) && !nrow(t))) next
    cat(sprintf("\n$%s%s\n", n, if (!is.null(attr(t, "title"))) paste0("  ", attr(t, "title")) else ""))
    print(t, digits = digits, title = FALSE)
  }
  invisible(x)
}

#' @rdname rtirt_summary
#' @export
as.data.frame.rtirt_summary <- function(x, row.names = NULL, optional = FALSE, ...) {
  p <- x$parameters; attr(p, "title") <- NULL; class(p) <- "data.frame"; p
}

#' @rdname rtirt_summary
#' @export
print.rtirt_table <- function(x, digits = 3, title = TRUE, ...) {
  if (title && !is.null(attr(x, "title"))) cat(attr(x, "title"), "\n")
  d <- x; attr(d, "title") <- NULL; class(d) <- "data.frame"
  print_table(d, digits)
  invisible(x)
}

# blocks of the parameter table
par_blocks <- c(a = "Discrimination (a)", b = "Difficulty (b)", lambda = "Time intensity (lambda)",
                sigma2t = "Residual variance of log RT (sigma2t)", rho = "Cross-relation of ability with log RT (rho)",
                structural = "Structural")

# long parameter table: block, item, est, se, lower, upper, then rhat, ess, prior (Gibbs) or
# z, pvalue (ECM; prior when prior = "map")
rtirt_partable <- function(e, bayes, ald = FALSE, title) {
  blk <- sub("\\[.*$", "", e$parameter); blk[!blk %in% names(par_blocks)[1:5]] <- "structural"
  item <- ifelse(blk == "structural", e$parameter, sub("^[^[]*\\[(.*)\\]$", "\\1", e$parameter))
  out <- data.frame(parameter = e$parameter, block = blk, item = item, est = e$est,
                    se = if (bayes) e$sd else e$se, lower = e$q025, upper = e$q975, stringsAsFactors = FALSE)
  if (bayes) { out$rhat <- e$rhat; out$ess <- e$ess }
  else { out$z <- e$z; out$pvalue <- 2 * stats::pnorm(-abs(e$z)) }
  if ("prior" %in% names(e)) out$prior <- e$prior
  out$fixed <- is.na(out$se) | out$se == 0
  p <- rtable(out, title, "rtirt_partable")
  attr(p, "bayes") <- bayes; attr(p, "ald") <- ald
  p
}

#' @export
print.rtirt_partable <- function(x, digits = 3, title = TRUE, ...) {
  if (title && !is.null(attr(x, "title"))) cat(attr(x, "title"), "\n")
  bayes <- isTRUE(attr(x, "bayes")); p <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  num <- if (bayes) c("est", "se", "lower", "upper", "rhat", "ess") else c("est", "se", "z", "pvalue", "lower", "upper")
  lab <- c(est = "Estimate", se = if (bayes) "Post.SD" else "Std.Err", lower = if (bayes) "pi.lower" else "ci.lower",
           upper = if (bayes) "pi.upper" else "ci.upper", rhat = "Rhat", ess = "ESS", z = "z-value", pvalue = "P(>|z|)")
  cells <- vapply(num, function(n) {
    v <- p[[n]]
    if (n == "ess") ifelse(is.na(v), "", formatC(round(v), format = "d")) else ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  }, character(nrow(p)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  cells[p$fixed, setdiff(num, "est")] <- ""
  pri <- if ("prior" %in% names(p)) ifelse(p$fixed, "", p$prior) else NULL
  w0 <- max(nchar(paste0("    ", p$item)), 20)
  w <- pmax(nchar(lab[num]), apply(cells, 2, function(v) max(nchar(v), 0)))
  line <- function(first, vals, prior = NULL) {
    s <- paste0(sprintf("%-*s", w0, first), paste(sprintf("%*s", w, vals), collapse = "  "))
    if (!is.null(prior) && nzchar(prior)) s <- paste0(s, "  ", prior)
    sub("\\s+$", "", s)
  }
  for (b in intersect(names(par_blocks), unique(p$block))) {
    i <- which(p$block == b)
    ttl <- if (b == "sigma2t" && isTRUE(attr(x, "ald"))) "ALD scale of log RT (sigma2t)" else par_blocks[[b]]
    cat(sprintf("\n%s:\n", ttl))
    cat(line("", lab[num], if (!is.null(pri)) "Prior"), "\n", sep = "")
    for (k in i) cat(line(paste0("    ", p$item[k]), cells[k, ], if (!is.null(pri)) pri[k]), "\n", sep = "")
  }
  invisible(x)
}

# wide item tables (est and se), from the parameter table
rtirt_item_wide <- function(p, items) {
  build <- function(col) {
    out <- data.frame(item = items, stringsAsFactors = FALSE)
    for (b in c("a", "b", "lambda", "sigma2t", "rho")) {
      r <- p[p$block == b, ]
      if (nrow(r)) out[[b]] <- r[[col]][match(items, r$item)]
    }
    out
  }
  list(est = build("est"), se = build("se"))
}

rtirt_header <- function(x) {
  D <- x$data; C <- x$cond
  c(model_label(x),
    sprintf("  %d persons x %d items%s", nrow(D$Y), ncol(D$Y),
            if (.model_info[[class(x)[1]]]$X) sprintf("; covariates: %s", paste(colnames(D$X), collapse = ", ")) else ""),
    if (has_em(x)) {
      f <- x$ecm
      sprintf("  %s (ECM, %d Gauss-Hermite nodes): %d iterations%s, %.1f s",
              if (f$map) "marginal posterior mode (item priors of gibbs())" else "maximum marginal likelihood", f$Q, f$iterations,
              if (f$converged) "" else " (not converged)", f$secs)
    } else sprintf("  %s: %d chains x %d iterations (burn-in %d, thin %d); %.1f s",
                     if (isTRUE(x$post$settings$collapse)) "Gibbs (partially collapsed)" else "Gibbs", C$n_chain, C$n_iter, C$n_burnin, C$n_thin, x$post$secs))
}
