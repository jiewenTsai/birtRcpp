# Summary objects: tables computed once, printed in the layout of tppcm::irt_pars() and
# blavaan --------------------------------------------------------------------------------
#
# summary(fit) returns an "rtirt_summary": header lines and titled tables (tibbles of
# class "rtirt_table"): fit, items and se (one row per item), parameters (one row per
# parameter, blocks a / b / lambda / sigma2t / rho / structural, with the prior of each for
# Gibbs fits), covariance and reliability. s$items etc. extract them, as.data.frame(s) gives
# the parameter table; print() shows each table under its name.

rtable <- function(df, title, cls = NULL) {
  if (is.null(df)) {
    return(NULL)
  }
  table <- tibble::as_tibble(as.data.frame(df))
  structure(table, class = c(cls, "rtirt_table", class(table)), title = title)
}

# a table as a plain data frame, without its title
plain_table <- function(x) {
  attr(x, "title") <- NULL
  class(x) <- "data.frame"
  x
}

#' Summary tables
#'
#' `summary()` of a sampled model (or one fitted by [ecm()]) returns its tables as one
#' object: `header`, `fit` (fit indices), `items` and `se` (one row per item: a, b, lambda,
#' sigma2t, rho and their posterior SDs or standard errors), `parameters` (one row per
#' parameter with its interval and, for Gibbs fits, its prior, as in blavaan), `covariance`
#' (person covariance) and `reliability`. Each element is a tibble; `as.data.frame()`
#' returns `parameters`. `print()` shows the tables under their names, so `s$items` gets the
#' one printed as `$items`.
#' @param x An `rtirt_summary`, or one of its tables.
#' @param tables Names of the tables to print (default all); an unknown name is an error.
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
  table_names <- setdiff(names(x), "header")
  if (!is.null(tables)) {
    check_table_names(tables, table_names)
    table_names <- intersect(tables, table_names)
  }
  cat(x$header, sep = "\n")
  for (name in table_names) {
    table <- x[[name]]
    if (is.null(table) || (is.data.frame(table) && !nrow(table))) {
      next
    }
    title <- attr(table, "title")
    title_text <- if (is.null(title)) "" else paste0("  ", title)
    cat(sprintf("\n$%s%s\n", name, title_text))
    print(table, digits = digits, title = FALSE)
  }
  invisible(x)
}

# print(s, tables = "item"): an error with the closest table name
check_table_names <- function(tables, table_names) {
  unknown <- setdiff(tables, table_names)
  if (!length(unknown)) {
    return(invisible())
  }
  guesses <- vapply(unknown, function(name) {
    guess <- closest_name(name, table_names)
    if (is.null(guess)) "" else sprintf(" (did you mean \"%s\"?)", guess)
  }, "")
  stop(sprintf("unknown table %s; the tables are %s", paste0("\"", unknown, "\"", guesses, collapse = ", "),
               paste0("\"", table_names, "\"", collapse = ", ")), call. = FALSE)
}

#' @rdname rtirt_summary
#' @export
as.data.frame.rtirt_summary <- function(x, row.names = NULL, optional = FALSE, ...) {
  plain_table(x$parameters)
}

#' @rdname rtirt_summary
#' @export
print.rtirt_table <- function(x, digits = 3, title = TRUE, ...) {
  print_title(x, title)
  print_table(plain_table(x), digits)
  invisible(x)
}

print_title <- function(x, title) {
  if (title && !is.null(attr(x, "title"))) {
    cat(attr(x, "title"), "\n")
  }
}

# blocks of the parameter table
par_blocks <- c(a = "Discrimination (a)", b = "Difficulty (b)", lambda = "Time intensity (lambda)",
                sigma2t = "Residual variance of log RT (sigma2t)", rho = "Cross-relation of ability with log RT (rho)",
                structural = "Structural")
item_blocks <- names(par_blocks)[1:5]

# long parameter table: block, item, est, se, lower, upper, then rhat, ess, prior (Gibbs) or
# z, pvalue (ECM; prior when prior = "map")
rtirt_partable <- function(e, bayes, ald = FALSE, title) {
  block <- sub("\\[.*$", "", e$parameter)
  block[!block %in% item_blocks] <- "structural"
  item <- ifelse(block == "structural", e$parameter, sub("^[^[]*\\[(.*)\\]$", "\\1", e$parameter))
  se <- if (bayes) e$sd else e$se
  table <- data.frame(parameter = e$parameter, block = block, item = item, est = e$est,
                      se = se, lower = e$q025, upper = e$q975, stringsAsFactors = FALSE)
  if (bayes) {
    table$rhat <- e$rhat
    table$ess <- e$ess
  } else {
    table$z <- e$z
    table$pvalue <- 2 * stats::pnorm(-abs(e$z))
  }
  if ("prior" %in% names(e)) {
    table$prior <- e$prior
  }
  table$fixed <- is.na(table$se) | table$se == 0
  partable <- rtable(table, title, "rtirt_partable")
  attr(partable, "bayes") <- bayes
  attr(partable, "ald") <- ald
  partable
}

# column labels of the printed parameter table
partable_labels <- function(bayes) {
  c(est = "Estimate",
    se = if (bayes) "Post.SD" else "Std.Err",
    lower = if (bayes) "pi.lower" else "ci.lower",
    upper = if (bayes) "pi.upper" else "ci.upper",
    rhat = "Rhat", ess = "ESS", z = "z-value", pvalue = "P(>|z|)")
}

# formatted cells: a character matrix, one column per numeric column of the table
partable_cells <- function(table, columns, digits) {
  format_column <- function(name) {
    v <- table[[name]]
    formatted <- if (name == "ess") formatC(round(v), format = "d") else formatC(v, digits = digits, format = "f")
    ifelse(is.na(v), "", formatted)
  }
  cells <- vapply(columns, format_column, character(nrow(table)))
  if (is.null(dim(cells))) {
    cells <- matrix(cells, nrow = 1)
  }
  cells[table$fixed, setdiff(columns, "est")] <- ""
  cells
}

#' @export
print.rtirt_partable <- function(x, digits = 3, title = TRUE, ...) {
  print_title(x, title)
  bayes <- isTRUE(attr(x, "bayes"))
  table <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  columns <- if (bayes) c("est", "se", "lower", "upper", "rhat", "ess") else c("est", "se", "z", "pvalue", "lower", "upper")
  labels <- partable_labels(bayes)[columns]
  cells <- partable_cells(table, columns, digits)
  priors <- if ("prior" %in% names(table)) ifelse(table$fixed, "", table$prior) else NULL

  first_width <- max(nchar(paste0("    ", table$item)), 20)
  widths <- pmax(nchar(labels), apply(cells, 2, function(v) max(nchar(v), 0)))
  line <- function(first, values, prior = NULL) {
    text <- paste0(sprintf("%-*s", first_width, first), paste(sprintf("%*s", widths, values), collapse = "  "))
    if (!is.null(prior) && nzchar(prior)) {
      text <- paste0(text, "  ", prior)
    }
    sub("\\s+$", "", text)
  }

  for (block in intersect(names(par_blocks), unique(table$block))) {
    rows <- which(table$block == block)
    block_title <- if (block == "sigma2t" && isTRUE(attr(x, "ald"))) "ALD scale of log RT (sigma2t)" else par_blocks[[block]]
    cat(sprintf("\n%s:\n", block_title))
    prior_label <- if (!is.null(priors)) "Prior"
    cat(line("", labels, prior_label), "\n", sep = "")
    for (k in rows) {
      prior_k <- if (!is.null(priors)) priors[k]
      cat(line(paste0("    ", table$item[k]), cells[k, ], prior_k), "\n", sep = "")
    }
  }
  invisible(x)
}

# wide item tables (est and se), from the parameter table
rtirt_item_wide <- function(p, items) {
  build <- function(column) {
    wide <- data.frame(item = items, stringsAsFactors = FALSE)
    for (block in item_blocks) {
      rows <- p[p$block == block, ]
      if (nrow(rows)) {
        wide[[block]] <- rows[[column]][match(items, rows$item)]
      }
    }
    wide
  }
  list(est = build("est"), se = build("se"))
}

rtirt_header <- function(x) {
  data <- x$data
  cond <- x$cond
  covariates <- if (.model_info[[class(x)[1]]]$X) sprintf("; covariates: %s", paste(colnames(data$X), collapse = ", ")) else ""
  size <- sprintf("  %d persons x %d items%s", nrow(data$Y), ncol(data$Y), covariates)
  c(model_label(x), size, fit_line(x, cond))
}

# how the model was fitted: ECM or Gibbs, with its settings and run time
fit_line <- function(x, cond) {
  if (has_em(x)) {
    fit <- x$ecm
    method <- if (fit$map) "marginal posterior mode (item priors of gibbs())" else "maximum marginal likelihood"
    converged <- if (fit$converged) "" else " (not converged)"
    return(sprintf("  %s (ECM, %d Gauss-Hermite nodes): %d iterations%s, %.1f s",
                   method, fit$nodes, fit$iterations, converged, fit$secs))
  }
  sampler <- if (isTRUE(x$post$settings$collapse)) "Gibbs (partially collapsed)" else "Gibbs"
  sprintf("  %s: %d chains x %d iterations (burn-in %d, thin %d); %.1f s",
          sampler, cond$n_chain, cond$n_iter, cond$n_burnin, cond$n_thin, x$post$secs)
}
