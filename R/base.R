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
  if (is.null(n_burnin)) n_burnin <- round(n_iter / 2)
  stopifnot(n_burnin < n_iter, n_thin >= 1, q_rt > 0, q_rt < 1)
  structure(list(n_subj = n_subj, n_item = n_item, n_feat = n_feat, n_iter = n_iter, n_chain = n_chain,
                 n_burnin = n_burnin, n_thin = n_thin, n_rep = n_rep, q_ra = q_ra, q_rt = q_rt),
            class = "sim_cond")
}

#' @export
print.sim_cond <- function(x, ...) {
  d <- function(v) if (is.null(v)) "(from data)" else v
  cat("sim_cond\n")
  cat(sprintf("  persons %s, items %s, covariates %s\n", d(x$n_subj), d(x$n_item), d(x$n_feat)))
  cat(sprintf("  %d chains x %d iterations (burn-in %d, thin %d)\n", x$n_chain, x$n_iter, x$n_burnin, x$n_thin))
  invisible(x)
}

# a block of columns: names of columns in `data`, or a matrix / data frame
get_block <- function(x, data, what) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) {
    if (is.null(data)) stop(what, " is given as column names; pass the data frame as `data =`", call. = FALSE)
    miss <- setdiff(x, names(data))
    if (length(miss)) stop(what, ": columns not in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
    x <- data[x]
  }
  if (is.data.frame(x)) {
    bad <- names(x)[!vapply(x, function(v) is.numeric(v) || is.logical(v), TRUE)]
    if (length(bad)) stop(what, ": non-numeric columns: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  x <- as.matrix(x); storage.mode(x) <- "double"
  x
}

# default item names, zero-padded to sort correctly: item01, ..., item10 (item001 with 100 items or more)
item_names <- function(K) sprintf("item%0*d", max(2, nchar(K)), seq_len(K))

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
  if (!is.null(data) && !is.data.frame(data)) stop("`data` must be a data frame")
  Y <- get_block(resp, data, "resp")
  N <- nrow(Y); K <- ncol(Y)
  items <- colnames(Y) %||% item_names(K)
  if (anyNA(Y)) {
    rows <- rowSums(is.na(Y)) > 0
    stop(sprintf("resp has %d missing values in %d persons; birtRcpp needs complete data. Drop those persons first, e.g. data[stats::complete.cases(data[c(resp, time)]), ]",
                 sum(is.na(Y)), sum(rows)), call. = FALSE)
  }
  vals <- sort(unique(c(Y)))
  if (!all(vals %in% c(0, 1))) {
    if (length(vals) == 2) stop(sprintf("resp is coded %s/%s; recode it to 0/1 (e.g. resp - %s)", vals[1], vals[2], vals[1]), call. = FALSE)
    stop("resp must be 0/1 (values found: ", paste(utils::head(vals, 6), collapse = ", "), if (length(vals) > 6) ", ..." else "", ")", call. = FALSE)
  }
  const <- which(apply(Y, 2, function(v) length(unique(v)) == 1))
  if (length(const))
    stop(sprintf("items without variation (everyone answers the same): %s; their parameters are not identified, drop them",
                 paste(sprintf("%s = %s", items[const], Y[1, const]), collapse = ", ")), call. = FALSE)
  if (!is.null(time) && !is.null(log_time)) stop("give either `time` (seconds) or `log_time`, not both")
  T <- get_block(time, data, "time"); log_t <- get_block(log_time, data, "log_time")
  check_dim <- function(M, what) if (!is.null(M) && !identical(dim(M), dim(Y)))
    stop(sprintf("%s is %d x %d but resp is %d x %d", what, nrow(M), ncol(M), N, K), call. = FALSE)
  check_dim(T, "time"); check_dim(log_t, "log_time")
  if (!is.null(T)) {
    if (anyNA(T)) stop("time has missing values; drop those persons first", call. = FALSE)
    if (any(T <= 0)) stop("time must be positive seconds", if (any(T < 0)) "; negative values suggest log times: use log_time =" else "", call. = FALSE)
    if (max(T) < 15 && stats::median(T) < 6)
      warning(sprintf("time is in (%.2f, %.2f) with median %.2f: if these are already log times, pass them as log_time = (time is log-transformed)",
                      min(T), max(T), stats::median(T)), call. = FALSE)
    log_t <- log(T)
  } else if (!is.null(log_t)) {
    if (anyNA(log_t)) stop("log_time has missing values; drop those persons first", call. = FALSE)
    if (min(log_t) > 0 && stats::median(log_t) > 15)
      warning(sprintf("log_time has median %.1f: these look like raw seconds; pass them as time =", stats::median(log_t)), call. = FALSE)
    T <- exp(log_t)
  }
  if (!is.null(log_t)) {
    dimnames(T) <- dimnames(log_t) <- list(NULL, items)
    flat <- which(apply(log_t, 2, stats::sd) < 1e-8)
    if (length(flat)) warning(sprintf("response times without variation (everyone the same time): %s; their residual variance goes to 0, check the data",
                                      paste(items[flat], collapse = ", ")), call. = FALSE)
  }
  X <- get_block(cov, data, "cov")
  if (!is.null(X)) {
    if (nrow(X) != N) stop(sprintf("cov has %d rows; resp has %d persons", nrow(X), N), call. = FALSE)
    if (anyNA(X)) stop("cov has missing values; drop those persons first", call. = FALSE)
    if (is.null(colnames(X))) colnames(X) <- sprintf("x%d", seq_len(ncol(X)))
    m <- colMeans(X); s <- apply(X, 2, stats::sd)
    off <- abs(m) > 0.5 * pmax(s, 1e-8)
    if (any(off)) warning("covariates not centred (", paste(sprintf("%s mean %.2f", names(m)[off], m[off]), collapse = ", "),
                          "): the latent regressions have no intercept, so centre or standardize them ",
                          "(or use gibbs(intercept = TRUE))", call. = FALSE)
  }
  if (is.character(id) && length(id) == 1 && !is.null(data)) {
    if (!id %in% names(data)) stop("id column '", id, "' not in `data`")
    id <- data[[id]]
  }
  if (is.null(id)) id <- rownames(Y) %||% seq_len(N)
  if (length(id) != N) stop(sprintf("id has length %d; resp has %d persons", length(id), N))
  dimnames(Y) <- list(NULL, items)
  structure(list(Y = Y, kappa = Y - 0.5, T = T, log_t = log_t, X = X, items = items, id = id), class = "input_data")
}

#' @export
print.input_data <- function(x, ...) {
  cat(sprintf("input_data: %d persons x %d items%s%s\n", nrow(x$Y), ncol(x$Y),
              if (is.null(x$log_t)) "" else ", response times",
              if (is.null(x$X)) "" else sprintf(", %d covariates", ncol(x$X))))
  cat("  items:", paste(utils::head(x$items, 6), collapse = ", "), if (length(x$items) > 6) "..." else "", "\n")
  if (!is.null(x$X)) cat("  covariates:", paste(colnames(x$X), collapse = ", "), "\n")
  invisible(x)
}

#' Parameter container
#'
#' Mirrors `InputPara`: a list of the model parameters (`theta`, `a`, `b`, `zeta`, `lambda`,
#' `sigma2t`, `nu`, `beta`, `rho`, `b_theta`, `sigma_p`). Used for true values in
#' simulations and for posterior means (`fit$post$mean`).
#' @param ... Named parameter values.
#' @export
input_para <- function(...) structure(list(...), class = "input_para")

#' @export
print.input_para <- function(x, ...) {
  cat("input_para:", paste(sprintf("%s[%s]", names(x), vapply(x, function(v) paste(dim(as.matrix(v)), collapse = "x"), "")),
                           collapse = ", "), "\n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
