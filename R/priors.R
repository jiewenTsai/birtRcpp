# Priors of gibbs() and ecm(prior = "map") ----------------------------------------------------------

#' Priors of the Gibbs samplers
#'
#' Sets the prior distributions used by [gibbs()] (and, for the item parameters, by
#' `ecm(prior = "map")`). The defaults are common weakly informative priors that do not depend on the
#' data, and every full conditional of the samplers stays conjugate (pure Gibbs):
#'
#' | parameter | default | Gibbs step |
#' |---|---|---|
#' | a (discrimination) | N(1, 1) on (0, Inf) | truncated normal |
#' | b (difficulty) | N(0, 3^2), not truncated | normal |
#' | lambda (time intensity) | N(0, 10^2), not truncated | normal |
#' | sigma2t (residual variance of log RT) | its SD ~ half-t(3, 1) | two inverse-gamma steps |
#' | speed (residual) variance | its SD ~ half-t(3, 1) | two inverse-gamma steps |
#' | regression coefficients, rho, c | N(0, 1) | normal |
#'
#' The half-t prior on a standard deviation (Gelman, 2006) is the scale mixture of Huang & Wand
#' (2013): var | u ~ IG(df / 2, df / u), u ~ IG(1 / 2, 1 / scale^2). The samplers draw
#' u | var ~ IG((df + 1) / 2, df / var + 1 / scale^2) and then var | u, data ~ IG(df / 2 + n / 2,
#' df / u + SS / 2), both conjugate. The auxiliary u is drawn afresh from the current variance in
#' every iteration, so it is not stored.
#'
#' [rtirt_priors_julia()] gives the priors of ExtendedRtIrtModeling.jl (Tsai, 2025), with which the
#' dissertation results are reproduced: a ~ N+(1, 1), b ~ N(0, 1) on \[-4, 4\], lambda ~ N+ centred
#' at the mean and SD of the log response times, inverse-gamma(0.001, 0.001) variances.
#'
#' A variance with `c(df, scale)` has a half-t prior on its square root and with `c(shape, scale)`
#' an inverse-gamma prior on the variance (one conjugate step).
#' @references Gelman, A. (2006). Prior distributions for variance parameters in hierarchical
#'   models. *Bayesian Analysis, 1*, 515-534.
#'
#'   Huang, A., & Wand, M. P. (2013). Simple marginally noninformative prior distributions for
#'   covariance matrices. *Bayesian Analysis, 8*, 439-452.
#' @param a Discrimination, normal on (0, Inf): `c(mean, sd)`.
#' @param b Difficulty, normal: `c(mean, sd)` or `c(mean, sd, lower, upper)` (truncated).
#' @param lambda Time intensity, normal: `c(mean, sd)` or `c(mean, sd, lower)` (truncated below).
#'   Only [rtirt_priors_julia()] centres it on the data (mean and SD of all log response times,
#'   truncated at 0), to reproduce ExtendedRtIrtModeling.jl; every other prior is fixed in advance.
#' @param sigma2t Residual variance (or ALD scale) of the log response times: `c(df = , scale = )`
#'   (half-t on its square root) or `c(shape = , scale = )` (inverse gamma).
#' @param beta Regression coefficients of the latent regressions (including `b_ability` of
#'   `rtirt_latent()`), normal with mean 0: `c(sd)`.
#' @param rho Cross-relations of `rtirt_cross()`, normal with mean 0: `c(sd)`.
#' @param var_speed Speed variance (`rtirt_cross()`), residual speed variance (`rtirt_latreg()`,
#'   `rtirt_null()`) or its ALD scale (`rtirt_latent()`): `c(df = , scale = )` or
#'   `c(shape = , scale = )`, as `sigma2t`.
#' @param cov Covariance `c` of speed with ability in `rtirt_latreg()`/`rtirt_null()` with free
#'   speed variance, normal with mean 0: `c(sd)` (with `speed_var = "fixed"` the correlation is
#'   uniform on (-1, 1)).
#' @return An object of class `rtirt_priors`.
#' @examples
#' rtirt_priors()                                       # the defaults
#' rtirt_priors(a = c(mean = 1, sd = 0.5))              # a tighter discrimination prior
#' rtirt_priors_julia()                                 # the priors of the dissertation
#' @export
rtirt_priors <- function(a = c(mean = 1, sd = 1), b = c(mean = 0, sd = 3), lambda = c(mean = 0, sd = 10),
                         sigma2t = c(df = 3, scale = 1), beta = c(sd = 1), rho = c(sd = 1),
                         var_speed = c(df = 3, scale = 1), cov = c(sd = 1)) {
  variance_sets <- list(c("df", "scale"), c("shape", "scale"))
  a_prior <- named_prior(a, list(c("mean", "sd")), "a")
  b_values <- b_prior(b)
  lambda_values <- lambda_prior_values(lambda)
  sigma2t_prior <- named_prior(sigma2t, variance_sets, "sigma2t")
  var_speed_prior <- named_prior(var_speed, variance_sets, "var_speed")
  priors <- list(a = a_prior$x, b = b_values, lambda = lambda_values,
                 sigma2t = sigma2t_prior$x, sigma2t_family = variance_family(sigma2t_prior$set),
                 beta = sd_prior(beta, "beta"), rho = sd_prior(rho, "rho"),
                 var_speed = var_speed_prior$x, var_speed_family = variance_family(var_speed_prior$set),
                 cov = sd_prior(cov, "cov"))
  check_priors(priors)
  structure(priors, class = "rtirt_priors")
}

numeric_prior <- function(x, what) {
  if (!is.numeric(x) || anyNA(x)) {
    stop(what, " must be a numeric vector", call. = FALSE)
  }
  x
}

# a prior given by name (one of `sets` of element names) or by position (the first set)
named_prior <- function(x, sets, what) {
  x <- numeric_prior(x, what)
  given_names <- names(x)
  for (set in sets) {
    if (!is.null(given_names) && setequal(given_names, set)) {
      return(list(set = set, x = x[set]))
    }
  }
  if (is.null(given_names) && length(x) == length(sets[[1]])) {
    return(list(set = sets[[1]], x = stats::setNames(x, sets[[1]])))
  }
  choices <- vapply(sets, function(set) sprintf("c(%s)", paste(set, collapse = ", ")), "")
  stop(sprintf("%s needs the elements %s", what, paste(choices, collapse = " or ")), call. = FALSE)
}

# b: c(mean, sd) or c(mean, sd, lower, upper)
b_prior <- function(b) {
  b <- numeric_prior(b, "b")
  if (length(b) == 2) {
    b <- c(b, -Inf, Inf)
  }
  if (length(b) != 4) {
    stop("b must be c(mean, sd) or c(mean, sd, lower, upper)", call. = FALSE)
  }
  element_names <- c("mean", "sd", "lower", "upper")
  if (!is.null(names(b)) && all(nzchar(names(b))) && setequal(names(b), element_names)) {
    b <- b[element_names]
  }
  stats::setNames(as.numeric(b), element_names)
}

# lambda: c(mean, sd) or c(mean, sd, lower)
lambda_prior_values <- function(lambda) {
  if (is.null(lambda)) {
    stop("lambda must be c(mean, sd) or c(mean, sd, lower); a prior centred on the data is only used by rtirt_priors_julia()", call. = FALSE)
  }
  lambda <- numeric_prior(lambda, "lambda")
  if (length(lambda) == 2) {
    lambda <- c(lambda, -Inf)
  }
  if (length(lambda) != 3) {
    stop("lambda must be c(mean, sd) or c(mean, sd, lower)", call. = FALSE)
  }
  stats::setNames(as.numeric(lambda), c("mean", "sd", "lower"))
}

sd_prior <- function(x, what) {
  x <- numeric_prior(x, what)
  if (length(x) != 1) {
    stop(what, " must be c(sd)", call. = FALSE)
  }
  c(sd = unname(x))
}

variance_family <- function(set) {
  if (set[1] == "df") "half-t" else "inverse-gamma"
}

check_priors <- function(p) {
  sds <- c(p$a[2], p$b["sd"], p$lambda["sd"], p$beta, p$rho, p$cov)
  if (any(sds <= 0)) {
    stop("prior standard deviations must be positive", call. = FALSE)
  }
  if (any(c(p$sigma2t, p$var_speed) <= 0)) {
    stop("variance prior parameters (df, shape, scale) must be positive", call. = FALSE)
  }
  if (p$b["lower"] >= p$b["upper"] || p$b["mean"] < p$b["lower"] || p$b["mean"] > p$b["upper"]) {
    stop("b needs lower < mean < upper", call. = FALSE)
  }
}

#' @rdname rtirt_priors
#' @export
rtirt_priors_julia <- function() {
  priors <- rtirt_priors(a = c(mean = 1, sd = 1), b = c(mean = 0, sd = 1, lower = -4, upper = 4),
                         sigma2t = c(shape = 0.001, scale = 0.001), var_speed = c(shape = 0.001, scale = 0.001))
  priors["lambda"] <- list(NULL)                           # the only data-centred prior: see lambda_prior()
  priors
}

#' @export
print.rtirt_priors <- function(x, ...) {
  text <- prior_text(x, NULL)
  cat("Priors (gibbs(); item parameters also for ecm(prior = \"map\")):\n")
  width <- max(nchar(names(text)))
  for (name in names(text)) {
    cat(sprintf("  %-*s  %s\n", width, name, text[[name]]))
  }
  invisible(x)
}

# numeric vector in the order of the C++ struct Prior (lambda NaN: centred on the data in C++)
prior_vector <- function(p) {
  family_code <- function(family) if (family == "half-t") 1 else 0
  lambda <- if (is.null(p$lambda)) c(NA_real_, NA_real_, 0) else unname(p$lambda[c("mean", "sd", "lower")])
  c(unname(p$a[1]), unname(p$a[2]),
    p$b[["mean"]], p$b[["sd"]], p$b[["lower"]], p$b[["upper"]],
    lambda,
    family_code(p$sigma2t_family), unname(p$sigma2t[1]), unname(p$sigma2t[2]), p$beta[["sd"]], p$rho[["sd"]],
    family_code(p$var_speed_family), unname(p$var_speed[1]), unname(p$var_speed[2]), p$cov[["sd"]])
}

# lambda prior actually used: c(mean, sd, lower) as given, or centred at the mean and SD of log T
lambda_prior <- function(p, L) {
  if (!is.null(p$lambda)) {
    return(unname(p$lambda))
  }
  if (is.null(L)) {
    return(c(NA_real_, NA_real_, 0))
  }
  c(mean(L, na.rm = TRUE), stats::sd(as.vector(L), na.rm = TRUE), 0)
}

# text of each prior (for summaries and print); L = log times (NULL: describe the rule)
prior_text <- function(p, L) {
  number <- function(v) format(signif(v, 3), trim = TRUE, drop0trailing = TRUE)
  normal_args <- function(m, s) if (s == 1) sprintf("%s, 1", number(m)) else sprintf("%s, %s^2", number(m), number(s))
  truncation <- function(lo, hi) {
    if (!is.finite(lo) && !is.finite(hi)) {
      return("")
    }
    sprintf(" T(%s, %s)", if (is.finite(lo)) number(lo) else "", if (is.finite(hi)) number(hi) else "")
  }
  variance_text <- function(x, family) {
    if (family == "half-t") {
      return(sprintf("SD ~ half-t(%s, %s)", number(x[1]), number(x[2])))
    }
    sprintf("IG(%s, %s)", number(x[1]), number(x[2]))
  }
  lambda <- lambda_prior(p, L)
  lambda_text <- if (anyNA(lambda[1:2])) {
    "N+(mean of log T, SD of log T^2)"
  } else if (is.null(p$lambda)) {
    sprintf("N+(%.2f, %.2f^2)", lambda[1], lambda[2])
  } else {
    sprintf("N(%s)%s", normal_args(lambda[1], lambda[2]), truncation(lambda[3], Inf))
  }
  list(a = sprintf("N+(%s)", normal_args(p$a[1], p$a[2])),
       b = sprintf("N(%s)%s", normal_args(p$b[["mean"]], p$b[["sd"]]), truncation(p$b[["lower"]], p$b[["upper"]])),
       lambda = lambda_text,
       sigma2t = variance_text(p$sigma2t, p$sigma2t_family),
       beta = sprintf("N(%s)", normal_args(0, p$beta[["sd"]])),
       rho = sprintf("N(%s)", normal_args(0, p$rho[["sd"]])),
       var_speed = variance_text(p$var_speed, p$var_speed_family),
       cov = sprintf("N(%s)", normal_args(0, p$cov[["sd"]])))
}

# the priors stored with a fit (gibbs or em), or the defaults
fit_priors <- function(object) {
  object$post$settings$priors %||% object$ecm$priors %||% rtirt_priors()
}

# log prior density of a variance v and its derivative in v (ecm(prior = "map")): inverse gamma on v, or
# half-t on sqrt(v) written as a density of v (Jacobian 1 / (2 sqrt(v)))
var_lprior <- function(v, x, fam) {
  if (fam == "half-t") {
    df <- x[[1]]
    A <- x[[2]]
    lp <- -(df + 1) / 2 * log1p(v / (df * A^2)) - 0.5 * log(v)
    d <- -(df + 1) / 2 / (df * A^2 + v) - 0.5 / v
    return(list(lp = lp, d = d))
  }
  shape <- x[[1]]
  scale <- x[[2]]
  list(lp = -(shape + 1) * log(v) - scale / v, d = -(shape + 1) / v + scale / v^2)
}
