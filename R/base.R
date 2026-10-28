# Containers mirroring Base.pl.jl: estimation/simulation conditions, data, parameters.

#' Conditions of a simulation (or estimation) run
#'
#' Mirrors `setCond()` of ExtendedRtIrtModeling.jl (Julia argument names are accepted by the
#' alias [setCond()]). For real data it is not needed: the dimensions come from the data and
#' the MCMC settings are arguments of [gibbs()]. In simulations it carries the dimensions for
#' [sim_para()] / [sim_data()] and default MCMC settings.
#'
#' @param n_subj,n_item,n_feat Numbers of persons, items and covariates (`NULL`: taken from
#'   the data when a model is built).
#' @param n_iter,n_chain Iterations per chain and number of (independent) chains.
#' @param n_burnin Burn-in per chain (default `n_iter / 2`; the Julia version ignored this argument).
#' @param n_thin Thinning interval for the stored draws.
#' @param n_rep Replications in [run_simulation()].
#' @param q_ra Kept for compatibility (unused, as in the Julia package).
#' @param q_rt Quantile level used by the Julia aliases `GibbsRtIrtLatentQr()`,
#'   `GibbsRtIrtCrossQr()` (the snake_case constructors take `quantile =`).
#' @return An object of class `sim_cond`.
#' @examples
#' cond <- set_cond(n_subj = 500, n_item = 10, n_feat = 2, n_iter = 2000, n_chain = 3)
#' cond
#' @export
set_cond <- function(n_subj = NULL, n_item = NULL, n_feat = NULL, n_iter = 5000, n_chain = 4, n_burnin = NULL,
                     n_thin = 1, n_rep = 10, q_ra = 0.5, q_rt = 0.5) {
  n_burnin <- n_burnin %||% round(n_iter / 2)
  stopifnot(n_burnin < n_iter, n_thin >= 1, q_rt > 0, q_rt < 1)
  cond <- list(n_subj = n_subj, n_item = n_item, n_feat = n_feat, n_iter = n_iter, n_chain = n_chain,
               n_burnin = n_burnin, n_thin = n_thin, n_rep = n_rep, q_ra = q_ra, q_rt = q_rt)
  structure(cond, class = "sim_cond")
}

#' @export
print.sim_cond <- function(x, ...) {
  from_data <- function(value) if (is.null(value)) "(from data)" else value
  cat("sim_cond\n")
  cat(sprintf("  persons %s, items %s, covariates %s\n", from_data(x$n_subj), from_data(x$n_item), from_data(x$n_feat)))
  cat(sprintf("  %d chains x %d iterations (burn-in %d, thin %d)\n", x$n_chain, x$n_iter, x$n_burnin, x$n_thin))
  invisible(x)
}

# a block of columns: names of columns in `data`, or a matrix / data frame; `hint` is added to
# the message about non-numeric columns
get_block <- function(x, data, what, hint = "") {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.character(x)) {
    if (is.null(data)) {
      stop(what, " is given as column names; pass the data frame as `data =`", call. = FALSE)
    }
    missing_cols <- setdiff(x, names(data))
    if (length(missing_cols)) {
      stop(what, ": columns not in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    x <- data[x]
  }
  if (is.data.frame(x)) {
    is_number <- vapply(x, function(v) is.numeric(v) || is.logical(v), TRUE)
    if (!all(is_number)) {
      stop(what, ": non-numeric columns: ", paste(names(x)[!is_number], collapse = ", "), hint, call. = FALSE)
    }
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

# default item names, zero-padded to sort correctly: item01, ..., item10 (item001 with 100 items or more)
item_names <- function(K) {
  sprintf("item%0*d", max(2, nchar(K)), seq_len(K))
}

#' Data container
#'
#' Collects the item responses, response times and covariates of one test (mirrors
#' `InputData()` of the Julia package). Each block is either a matrix / data frame or, with
#' `data`, a vector of column names:
#'
#' ```
#' dat <- input_data(resp = items, time = paste0(items, "_S"), cov = covs, id = "IDSTUD", data = raw)
#' dat <- input_data(resp = Y, time = RT, cov = X)
#' ```
#'
#' Item names (the column names of `resp`) label all output. Missing values are not
#' supported: drop the incomplete persons first.
#'
#' The data must be in wide format: one row per person, one column per item (for responses
#' and for times). Data with one row per person and item (long format, as exported by
#' PsychoPy or E-Prime) are refused, recognized by a single item or by ids that repeat on
#' most rows; reshape them first, e.g.
#' `tidyr::pivot_wider(long, id_cols = id, names_from = item, values_from = c(correct, rt))`
#' or `stats::reshape(long, idvar = "id", timevar = "item", direction = "wide")`. A few
#' duplicated ids give a warning. Covariates must be numeric: dummy-code categorical ones
#' first, e.g. `model.matrix(~ grade, data)[, -1]` (k - 1 columns for k categories), and
#' centre them.
#'
#' Responses coded as text or factors (`"correct"` / `"wrong"`) must be recoded to 0/1 first;
#' the message shows the values found and the recoding. `TRUE` / `FALSE` are accepted.
#'
#' **Units of time.** Seconds are expected, but the unit only shifts the time intensities:
#' times in milliseconds give every `lambda` larger by log(1000) = 6.91 (minutes: smaller by
#' log(60) = 4.09), and a, b, rho, the residual variances and the correlations do not change.
#' Report the unit with `lambda`. The default prior of `lambda` in [gibbs()] (N(0, 10^2)) is
#' meant for seconds, so a median time above 1000 gives a warning that suggests milliseconds.
#'
#' @param resp Binary (0/1) responses, persons x items.
#' @param time Response times in seconds (positive), same dimensions as `resp`; they are
#'   log-transformed. Give log times as `log_time` instead.
#' @param cov Person covariates, persons x covariates. The latent regressions have no
#'   intercept (unless `gibbs(intercept = TRUE)`), so centre (best: standardize) them.
#' @param data Optional data frame from which `resp`, `time`, `cov` and `id` are taken by name.
#' @param id Person identifiers (a vector, or a column name of `data`); default `1..N`.
#' @param log_time Log response times, instead of `time`.
#' @return An object of class `input_data` with `Y`, `kappa = Y - 1/2`, `T`, `log_t`, `X`,
#'   `items` and `id`.
#' @export
input_data <- function(resp, time = NULL, cov = NULL, data = NULL, id = NULL, log_time = NULL) {
  if (!is.null(data) && !is.data.frame(data)) {
    stop("`data` must be a data frame")
  }
  Y <- get_block(resp, data, "resp", resp_text_hint(resp, data))
  items <- colnames(Y) %||% item_names(ncol(Y))
  check_resp(Y, items)

  times <- read_times(time, log_time, data, Y, items)
  X <- read_cov(cov, data, nrow(Y))
  id_arg <- id
  id <- read_id(id, data, Y)
  check_wide(Y, id, list(resp = resp, time = time %||% log_time, id = id_arg), data)

  dimnames(Y) <- list(NULL, items)
  container <- list(Y = Y, kappa = Y - 0.5, T = times$T, log_t = times$log_t, X = X, items = items, id = id)
  structure(container, class = "input_data")
}

# responses given as text or factors: the values found and how to recode them to 0/1
resp_text_hint <- function(resp, data) {
  block <- if (is.data.frame(resp)) {
    resp
  } else if (is.character(resp) && is.data.frame(data) && all(resp %in% names(data))) {
    data[resp]
  } else {
    return("")
  }
  is_text <- !vapply(block, function(v) is.numeric(v) || is.logical(v), TRUE)
  values <- unique(unlist(lapply(block[is_text], function(v) as.character(unique(v)))))
  values <- values[!is.na(values)]
  if (!length(values)) {
    return("")
  }
  correct <- values[grepl("^(correct|right|true|yes|1|c)$", tolower(values))]
  example <- if (length(correct)) correct[1] else values[1]
  more <- if (length(values) > 4) ", ..." else ""
  sprintf(" (values: %s%s); recode them to 0/1, e.g. data[resp] <- lapply(data[resp], function(x) as.integer(x == \"%s\"))",
          paste0("\"", utils::head(values, 4), "\"", collapse = ", "), more, example)
}

# missing values of resp or time: how many, in how many persons, and how to drop them
stop_missing <- function(what, M) {
  stop(sprintf("%s has %d missing values in %d persons; birtRcpp needs complete data. Drop those persons first, e.g. data[stats::complete.cases(data[c(resp, time)]), ]",
               what, sum(is.na(M)), sum(rowSums(is.na(M)) > 0)), call. = FALSE)
}

# responses: complete, coded 0/1, every item with both values
check_resp <- function(Y, items) {
  if (anyNA(Y)) {
    stop_missing("resp", Y)
  }
  values <- sort(unique(c(Y)))
  if (!all(values %in% c(0, 1))) {
    if (length(values) == 2) {
      stop(sprintf("resp is coded %s/%s; recode it to 0/1 (e.g. resp - %s)", values[1], values[2], values[1]), call. = FALSE)
    }
    more <- if (length(values) > 6) ", ..." else ""
    stop("resp must be 0/1 (values found: ", paste(utils::head(values, 6), collapse = ", "), more, ")", call. = FALSE)
  }
  constant <- which(apply(Y, 2, function(v) length(unique(v)) == 1))
  if (length(constant)) {
    constant_items <- paste(sprintf("%s = %s", items[constant], Y[1, constant]), collapse = ", ")
    stop(sprintf("items without variation (everyone answers the same): %s; their parameters are not identified, drop them",
                 constant_items), call. = FALSE)
  }
}

# response times in seconds (`time`) or on the log scale (`log_time`): both versions, checked
read_times <- function(time, log_time, data, Y, items) {
  if (!is.null(time) && !is.null(log_time)) {
    stop("give either `time` (seconds) or `log_time`, not both")
  }
  T <- get_block(time, data, "time")
  log_t <- get_block(log_time, data, "log_time")
  check_same_dim(T, Y, "time")
  check_same_dim(log_t, Y, "log_time")

  if (!is.null(T)) {
    check_seconds(T)
    log_t <- log(T)
  } else if (!is.null(log_t)) {
    check_log_seconds(log_t)
    T <- exp(log_t)
  }

  if (!is.null(log_t)) {
    dimnames(T) <- list(NULL, items)
    dimnames(log_t) <- list(NULL, items)
    warn_flat_times(log_t, items)
  }
  list(T = T, log_t = log_t)
}

check_same_dim <- function(M, Y, what) {
  if (!is.null(M) && !identical(dim(M), dim(Y))) {
    stop(sprintf("%s is %d x %d but resp is %d x %d", what, nrow(M), ncol(M), nrow(Y), ncol(Y)), call. = FALSE)
  }
}

check_seconds <- function(T) {
  if (anyNA(T)) {
    stop_missing("time", T)
  }
  if (any(T <= 0)) {
    hint <- if (any(T < 0)) "; negative values suggest log times: use log_time =" else ""
    stop("time must be positive seconds", hint, call. = FALSE)
  }
  if (stats::median(T) > 1000) {
    warning(sprintf("time has median %.0f: if these are milliseconds, divide them by 1000 (seconds). The unit shifts only the time intensities (lambda by log(1000) = 6.91), but the default prior of lambda is meant for seconds",
                    stats::median(T)), call. = FALSE)
  }
  if (max(T) < 15 && stats::median(T) < 6) {
    warning(sprintf("time is in (%.2f, %.2f) with median %.2f: if these are already log times, pass them as log_time = (time is log-transformed)",
                    min(T), max(T), stats::median(T)), call. = FALSE)
  }
}

check_log_seconds <- function(log_t) {
  if (anyNA(log_t)) {
    stop_missing("log_time", log_t)
  }
  if (min(log_t) > 0 && stats::median(log_t) > 15) {
    warning(sprintf("log_time has median %.1f: these look like raw seconds; pass them as time =", stats::median(log_t)), call. = FALSE)
  }
}

warn_flat_times <- function(log_t, items) {
  flat <- which(apply(log_t, 2, stats::sd) < 1e-8)
  if (length(flat)) {
    warning(sprintf("response times without variation (everyone the same time): %s; their residual variance goes to 0, check the data",
                    paste(items[flat], collapse = ", ")), call. = FALSE)
  }
}

# covariates: complete, one row per person, named; a warning when they are not centred
read_cov <- function(cov, data, n_person) {
  hint <- paste0("; covariates must be numeric: dummy-code categorical ones first, e.g. model.matrix(~ grade, data)[, -1] ",
                 "(k - 1 columns for k categories) or data$female <- as.numeric(data$gender == \"F\"), then centre them")
  X <- get_block(cov, data, "cov", hint)
  if (is.null(X)) {
    return(NULL)
  }
  if (nrow(X) != n_person) {
    stop(sprintf("cov has %d rows; resp has %d persons", nrow(X), n_person), call. = FALSE)
  }
  if (anyNA(X)) {
    stop(sprintf("cov has %d missing values in %d persons; birtRcpp needs complete data. Drop those persons first, e.g. data[stats::complete.cases(data[c(resp, time, cov)]), ]",
                 sum(is.na(X)), sum(rowSums(is.na(X)) > 0)), call. = FALSE)
  }
  if (is.null(colnames(X))) {
    colnames(X) <- sprintf("x%d", seq_len(ncol(X)))
  }
  cov_mean <- colMeans(X)
  cov_sd <- apply(X, 2, stats::sd)
  off_centre <- abs(cov_mean) > 0.5 * pmax(cov_sd, 1e-8)
  if (any(off_centre)) {
    means <- paste(sprintf("%s mean %.2f", names(cov_mean)[off_centre], cov_mean[off_centre]), collapse = ", ")
    warning("covariates not centred (", means,
            "): the latent regressions have no intercept, so centre or standardize them ",
            "(or use gibbs(intercept = TRUE))", call. = FALSE)
  }
  X
}

# person identifiers: a vector, a column name of `data`, the row names of resp, or 1..N
read_id <- function(id, data, Y) {
  n_person <- nrow(Y)
  if (is.character(id) && length(id) == 1 && !is.null(data)) {
    if (!id %in% names(data)) {
      stop("id column '", id, "' not in `data`")
    }
    id <- data[[id]]
  }
  id <- id %||% rownames(Y) %||% seq_len(n_person)
  if (length(id) != n_person) {
    stop(sprintf("id has length %d; resp has %d persons", length(id), n_person))
  }
  id
}

# one row per person: long data (one row per person x item) show up as one item or as ids that
# repeat on most rows, and are refused with how to reshape them; a few repeated ids give a warning
# (scores() returns the ids for merge(), which needs them unique)
check_wide <- function(Y, id, given, data) {
  n_row <- nrow(Y)
  repeated <- duplicated(id)
  n_id <- sum(!repeated)
  if (ncol(Y) == 1 || n_id <= n_row / 2) {
    shape <- sprintf("%d item%s", ncol(Y), if (ncol(Y) > 1) "s" else "")
    if (any(repeated)) {
      shape <- sprintf("%s, %d distinct ids on %d rows", shape, n_id, n_row)
    }
    id_name <- long_name(given$id, "id")
    stop(sprintf(paste0("the data look like long format (%s): input_data() needs one row per person and one column per item. ",
                        "Reshape first, e.g. wide <- %s, then input_data(resp = ..., time = ..., id = \"%s\", data = wide); ",
                        "or stats::reshape(..., idvar = \"%s\", timevar = \"%s\", direction = \"wide\")"),
                 shape, long_reshape_call(given, data, id), id_name, id_name, long_item_column(given, data, id)),
         call. = FALSE)
  }
  if (any(repeated)) {
    shown <- paste(utils::head(unique(id[repeated]), 3), collapse = ", ")
    warning(sprintf("id has %d repeated value%s (%s): each person should be one row; fix or drop the duplicates (merge(scores(fit), data) by id needs unique ids)",
                    sum(repeated), if (sum(repeated) > 1) "s" else "", shown), call. = FALSE)
  }
}

long_name <- function(x, default) {
  if (is.character(x) && length(x) == 1) x else default
}

# the column of `data` that names the item: as many distinct values as rows per person
long_item_column <- function(given, data, id) {
  per_id <- max(table(id))
  if (is.null(data) || per_id < 2) {
    return("item")
  }
  used <- unlist(Filter(is.character, given))
  candidates <- setdiff(names(data), used)
  matches <- candidates[vapply(candidates, function(v) length(unique(data[[v]])) == per_id, TRUE)]
  if (length(matches)) matches[1] else "item"
}

long_reshape_call <- function(given, data, id) {
  values <- c(long_name(given$resp, "resp"), if (!is.null(given$time)) long_name(given$time, "rt"))
  values_text <- if (length(values) > 1) sprintf("c(%s)", paste(values, collapse = ", ")) else values
  sprintf("tidyr::pivot_wider(data, id_cols = %s, names_from = %s, values_from = %s)", long_name(given$id, "id"),
          long_item_column(given, data, id), values_text)
}

#' @export
print.input_data <- function(x, ...) {
  times <- if (is.null(x$log_t)) "" else ", response times"
  covariates <- if (is.null(x$X)) "" else sprintf(", %d covariates", ncol(x$X))
  cat(sprintf("input_data: %d persons x %d items%s%s\n", nrow(x$Y), ncol(x$Y), times, covariates))
  more_items <- if (length(x$items) > 6) "..." else ""
  cat("  items:", paste(utils::head(x$items, 6), collapse = ", "), more_items, "\n")
  if (!is.null(x$X)) {
    cat("  covariates:", paste(colnames(x$X), collapse = ", "), "\n")
  }
  invisible(x)
}

#' Parameter container
#'
#' Mirrors `InputPara`: a list of the model parameters (`theta`, `a`, `b`, `zeta`, `lambda`,
#' `sigma2t`, `nu`, `beta`, `rho`, `b_theta`, `sigma_p`). Used for true values in
#' simulations and for posterior means (`fit$post$mean`).
#' @param ... Named parameter values.
#' @export
input_para <- function(...) {
  structure(list(...), class = "input_para")
}

#' @export
print.input_para <- function(x, ...) {
  dims <- vapply(x, function(v) paste(dim(as.matrix(v)), collapse = "x"), "")
  cat("input_para:", paste(sprintf("%s[%s]", names(x), dims), collapse = ", "), "\n")
  invisible(x)
}

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
