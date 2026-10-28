# Printed numbers ----------------------------------------------------------------------------
#
# Every print() and summary() of the package shows numbers the same way: 3 decimals, fixed
# notation (formatC(x, format = "f", digits = 3)); counts (persons, items, npar, Df, ESS, ...)
# without decimals; p-values with 3 decimals and "<0.001" below 0.001; NA blank in tables. Only
# the display is rounded: the objects keep full precision (x$est, coef(fit), as.data.frame(x)).
# Run times keep one decimal and percentages of print() one decimal; the tolerances of
# check_ecm() stay in scientific notation (1e-08 shown as 0.000 would say nothing).

# columns of whole numbers (shown without decimals when every value is whole)
count_columns <- c("n", "N", "npar", "Df", "df", "ess", "ess_tail", "nodes", "iterations", "cycles", "rep", "id",
                   "n_chain", "chains", "persons", "items")

# columns of p-values ("<0.001" below 0.001); p_flag (a posterior probability) and p_waic (a
# number of parameters) are not p-values
pvalue_columns <- c("p.value", "p.holm", "pvalue", "Pr(>Chisq)", "p_chi_pf", "p_Lambda_s", "p_L_s")

# x with `digits` decimals, no scientific notation, NA as `na`
fmt_fixed <- function(x, digits = 3, na = "") {
  out <- formatC(as.numeric(x), format = "f", digits = digits)
  out[is.na(x)] <- na
  out
}

# whole numbers without decimals
fmt_count <- function(x, na = "") {
  out <- formatC(round(as.numeric(x)), format = "d", big.mark = "")
  out[is.na(x)] <- na
  out
}

# p-values: 3 decimals, "<0.001" below 0.001
fmt_p <- function(p, digits = 3, na = "") {
  limit <- 10^-digits
  out <- ifelse(p < limit, paste0("<", formatC(limit, format = "f", digits = digits)), formatC(p, format = "f", digits = digits))
  out[is.na(p)] <- na
  out
}

# "p = 0.123" or "p < 0.001", for the lines of print()
p_text <- function(p, digits = 3) {
  text <- fmt_p(p, digits, na = "NA")
  if (startsWith(text, "<")) paste("p", sub("<", "< ", text, fixed = TRUE)) else paste("p =", text)
}

is_whole <- function(x) {
  values <- x[is.finite(x)]
  all(abs(values - round(values)) < 1e-8)
}

# the cells of one column of a printed table
fmt_column <- function(x, name = "", digits = 3) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (!is.numeric(x)) {
    return(ifelse(is.na(x), "", as.character(x)))
  }
  if (name %in% pvalue_columns) {
    return(fmt_p(x, digits))
  }
  if (is.integer(x) || (name %in% count_columns && is_whole(x))) {
    return(fmt_count(x))
  }
  fmt_fixed(x, digits)
}

# plain table (as in birt): numbers right-aligned with 3 decimals, text left-aligned; the lines
format_table <- function(df, digits = 3) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  header <- names(df)
  cells <- vapply(seq_along(df), function(j) fmt_column(df[[j]], header[j], digits), character(nrow(df)))
  cells <- matrix(cells, nrow = nrow(df), ncol = ncol(df))
  widths <- pmax(nchar(header), apply(rbind(cells, ""), 2, function(x) max(nchar(x))))
  left <- !vapply(df, is.numeric, TRUE)
  format_row <- function(x) {
    padded <- ifelse(left, sprintf("%-*s", widths, x), sprintf("%*s", widths, x))
    sub("\\s+$", "", paste(padded, collapse = "  "))
  }
  c(format_row(header), if (nrow(df)) apply(cells, 1, format_row))
}

print_table <- function(df, digits = 3) {
  cat(format_table(df, digits), sep = "\n")
}

# ---- tables returned to the user -------------------------------------------------------------

# estimates(), scores(), fit_indices(), compare_para(), person_fit(), ... return tibbles (or data
# frames) of class "rtirt_tbl": printed with 3 decimals, used as any tibble (x$est, dplyr,
# as.data.frame()) with full precision
as_rtirt_tbl <- function(x) {
  class(x) <- unique(c("rtirt_tbl", class(x)))
  x
}

#' Printing the tables of birtRcpp
#'
#' The tables of [estimates()], [scores()], [fit_indices()], [compare_para()], [person_fit()],
#' [preknowledge_test()], [dic()] and [run_simulation()] are tibbles (or data frames) of class
#' `rtirt_tbl`. They print numbers with 3 decimals (counts without decimals, p-values below 0.001
#' as `<0.001`, `NA` blank), and all rows of a table with up to 40 rows (else the first 20). The
#' values are not rounded: `x$est`, dplyr verbs and `as.data.frame(x)` have the full precision.
#' Named vectors and matrices of numbers ([coef()][coef.rtirt], [reliability()],
#' [marginal_loglik()]) print with 3 decimals as well.
#' @param x A table of birtRcpp.
#' @param digits Decimals.
#' @param n Rows to show (default: all up to 40 rows, else 20); `Inf` for all.
#' @param ... Unused.
#' @return `format()` the lines; `print()` returns `x` invisibly.
#' @examples
#' cond <- set_cond(n_subj = 200, n_item = 5)
#' fit <- ecm(sim_data(cond, sim_para(cond, "null"), "null"), model = "null", se = FALSE)
#' est <- estimates(fit)
#' est            # 3 decimals
#' est$est[1]     # full precision
#' as.data.frame(est)[1:2, ]
#' @name rtirt_tbl
NULL

#' @rdname rtirt_tbl
#' @export
format.rtirt_tbl <- function(x, digits = 3, n = NULL, ...) {
  n_row <- nrow(x)
  n <- n %||% if (n_row > 40) 20 else n_row
  shown <- min(n, n_row)
  kind <- if (inherits(x, "tbl_df")) "A tibble" else "A data frame"
  header <- sprintf("# %s: %d x %d", kind, n_row, ncol(x))
  table <- as.data.frame(x, stringsAsFactors = FALSE)[seq_len(shown), , drop = FALSE]
  more <- if (n_row > shown) sprintf("# ... %d more rows: print(x, n = Inf) shows all", n_row - shown)
  c(header, format_table(table, digits), more)
}

#' @rdname rtirt_tbl
#' @export
print.rtirt_tbl <- function(x, digits = 3, n = NULL, ...) {
  cat(format(x, digits = digits, n = n), sep = "\n")
  invisible(x)
}

# named vectors and matrices of numbers (coef(), reliability(), marginal_loglik())
as_rtirt_values <- function(x) {
  x <- unclass(x)
  class(x) <- c("rtirt_values", if (is.matrix(x)) c("matrix", "array") else "numeric")
  x
}

#' @rdname rtirt_tbl
#' @export
print.rtirt_values <- function(x, digits = 3, ...) {
  values <- unclass(x)
  if (!is.numeric(values)) {
    print(values, ...)
    return(invisible(x))
  }
  cells <- fmt_fixed(values, digits, na = "NA")
  attributes(cells) <- attributes(values)[intersect(names(attributes(values)), c("names", "dim", "dimnames"))]
  if (is.matrix(cells)) {                       # column labels right-aligned over the numbers
    labels <- colnames(cells) %||% sprintf("[,%d]", seq_len(ncol(cells)))
    widths <- pmax(nchar(labels), apply(rbind(cells, ""), 2, function(v) max(nchar(v))))
    colnames(cells) <- sprintf("%*s", widths, labels)
  }
  print(noquote(cells), right = TRUE)
  invisible(x)
}

#' @export
`[.rtirt_values` <- function(x, ...) {
  out <- NextMethod()
  if (is.numeric(out)) as_rtirt_values(out) else out
}
