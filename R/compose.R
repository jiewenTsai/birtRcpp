# Helpers for composing samplers from the blocks: the steps that the examples ex08-ex12 used to write
# by hand. Notation: persons j = 1..N (rows), items i = 1..K (columns).

#' Helpers for composing samplers
#'
#' General steps that come back in many samplers written from the blocks ([blocks]); with them the
#' samplers of the examples `ex08`-`ex12` (testlet, multidimensional, conditional dependence, LNIRT,
#' variable speed) are a few lines per step. Persons j are the rows and items i the columns of the
#' matrices.
#'
#' * `calc_normal_part(x, z, w, rest, by)`: the part of a full conditional that comes from normal
#'   (pseudo-)observations \eqn{z_{ij} = x_{ij}\beta + \text{rest}_{ij} + e_{ij}}, \eqn{e_{ij} \sim
#'   N(0, 1/w_{ij})}, for a parameter \eqn{\beta} of each item (`by = "item"`, sums over persons) or of
#'   each person (`by = "person"`, sums over items): precision \eqn{\sum w x^2} and linear term
#'   \eqn{\sum w x (z - \text{rest})}. Add the prior and draw with [draw_gauss()] or [draw_tnorm()].
#'   Pólya–Gamma: `z = (Y - 1/2) / omega`, `w = omega`; probit: `z` from [draw_probit_z()], `w = 1`;
#'   log times: `z = log t`, `w = 1 / sigma2`. Missing `z` (NA) are skipped.
#' * `draw_mvn(prec, num, lower, upper)`: multivariate normal \eqn{N(Q^{-1}h, Q^{-1})} from the
#'   precision `prec` (Q) and linear term `num` (h); `prec` can be an array p x p x n and `num` a
#'   p x n matrix for n independent draws (one per item or person), optionally truncated to the box
#'   `[lower, upper]` (rejection; a slice that fails is NA).
#' * `draw_items_mh(Y, theta, a, b, prior_mean, prior_prec, link, form)`: item parameters (a, b) with
#'   the latent responses integrated out, Metropolis–Hastings with a Fisher-scoring proposal (the
#'   cure for the slow mixing of easy and difficult items under data augmentation); logit or probit,
#'   `form = "difficulty"` for a (theta - b) or `"intercept"` for a theta - b; prior
#'   N((a, b); `prior_mean` (a row per item or one row), `prior_prec`^-1). Returns `list(a, b)` with the
#'   attribute `accepted` (as [draw_items_joint()]); draw the latent responses anew afterwards.
#' * `draw_invwishart(df, scale)`: inverse Wishart.
#' * `draw_mean_mvn(X, Sigma, mean0, cov0)`: the common mean of the rows of `X` ~ N(mu, Sigma) with
#'   prior mu ~ N(mean0, cov0) (hyperparameters of item- or person-level distributions).
#' * `draw_cor(x1, x2, cor, steps, step, method)`: the correlation of two standard-normal variables
#'   (unit variances, the identification of several abilities), `steps` Metropolis steps on
#'   atanh(cor) with a uniform prior: random walk (`"rw"`) or MALA (`"mala"`, see [draw_mala()]).
#' * `calc_person_fit(L, fitted, sigma2)`: the RT person-fit statistic of Marianti et al. (2014)
#'   \eqn{l^t_j = \sum_i (\log t_{ij} - \text{fitted}_{ij})^2 / \sigma_i^2} per person, with the number
#'   of observed times `df` (missing times are skipped); chi-square with `df` degrees of freedom.
#' @param x,z,w,rest Numbers or N x K matrices (persons x items).
#' @param by `"item"` or `"person"`.
#' @param prec,num Precision matrix and linear term (or an array of precisions and a matrix of
#'   linear terms, one slice / column per draw).
#' @param lower,upper Bounds of the truncation box (recycled to the dimension).
#' @param max_tries Rejection attempts per draw.
#' @param Y Responses (0/1, NA for missing), N x K.
#' @param theta Abilities.
#' @param a,b Current item parameters.
#' @param prior_mean,prior_prec Prior mean (K x 2 or length 2) and prior precision (2 x 2) of (a, b).
#' @param link `"logit"` or `"probit"`.
#' @param form `"difficulty"` or `"intercept"`.
#' @param df,scale Degrees of freedom and scale matrix.
#' @param X Matrix whose rows share the mean.
#' @param Sigma Covariance of the rows.
#' @param mean0,cov0 Prior mean and covariance of the common mean.
#' @param x1,x2 The two variables (standardized by the model).
#' @param cor Current correlation.
#' @param steps,step Number and size of the Metropolis steps.
#' @param method `"rw"` (random walk) or `"mala"`.
#' @param L,fitted,sigma2 Log times, their fitted means (N x K) and the residual variances (per item).
#' @return `calc_normal_part()`: a normal part (`prec`, `num`; see [normal_parts]); `draw_mvn()`: a vector or a p x n
#'   matrix; `draw_items_mh()`: a list; `draw_invwishart()`: a matrix; `draw_mean_mvn()`: a vector;
#'   `draw_cor()`: a number with attribute `accepted`; `calc_person_fit()`: a data frame.
#' @references Marianti, S., Fox, J.-P., Avetisyan, M., Veldkamp, B. P., & Tijmstra, J. (2014).
#'   Testing for aberrant behavior in response time modeling. *Journal of Educational and
#'   Behavioral Statistics, 39*, 426-451.
#' @examples
#' set.seed(1)
#' N <- 300; K <- 5; theta <- rnorm(N); a <- rep(1.2, K); b <- seq(-1, 1, length.out = K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(theta, a), 2, a * b))), N)
#' om <- draw_omega(a, b, theta)
#' # b_i | omega: eta = a (theta - b) enters the PG pseudo-observations with x = -a
#' g <- calc_normal_part(matrix(-a, N, K, byrow = TRUE), (Y - 0.5) / om, om, rest = outer(theta, a))
#' draw_gauss(1 / 9 + g$prec, g$num)                       # prior N(0, 3^2)
#' # (a, b) jointly with the latent responses integrated out
#' draw_items_mh(Y, theta, a, b, prior_mean = c(1, 0), prior_prec = diag(c(1, 1 / 9)))
#' draw_mvn(diag(2) * 4, c(1, 2), lower = 0)
#' @name compose
NULL

#' @rdname compose
#' @export
calc_normal_part <- function(x, z, w = 1, rest = 0, by = c("item", "person")) {
  by <- match.arg(by); z <- as.matrix(z); N <- nrow(z); K <- ncol(z)
  chk <- function(v) {
    if (length(v) == 1 || identical(dim(v), c(N, K))) return(v)
    if (is.null(dim(v)) && length(v) == N * K) return(matrix(v, N, K))      # a vector when z has one column
    stop(sprintf("x, w and rest must be numbers or %d x %d matrices like z (got length %d); a value per item: matrix(v, N, K, byrow = TRUE)",
                 N, K, length(v)), call. = FALSE)
  }
  x <- chk(x); w <- chk(w); rest <- chk(rest)
  R <- z - rest; W <- w
  if (anyNA(R) || anyNA(W)) { obs <- !is.na(R) & !is.na(W); W <- ifelse(obs, W, 0); R <- ifelse(obs, R, 0) }
  WX <- W * x; if (length(WX) == 1) WX <- matrix(WX, N, K)
  s <- if (by == "item") colSums else rowSums
  normal_part(s(WX * x), s(WX * R))
}

#' @rdname compose
#' @export
draw_mvn <- function(prec, num, lower = -Inf, upper = Inf, max_tries = 100) {
  single <- is.matrix(prec) && length(dim(prec)) == 2
  Q <- if (single) array(prec, c(dim(prec), 1)) else prec
  p <- dim(Q)[1]; H <- matrix(num, nrow = p)
  if (dim(Q)[3] != ncol(H)) stop("prec has ", dim(Q)[3], " slices but num ", ncol(H), " columns", call. = FALSE)
  x <- .blk_mvn(Q, H, rep_len(lower, p), rep_len(upper, p), as.integer(max_tries))
  if (single) drop(x) else x
}

#' @rdname compose
#' @export
draw_items_mh <- function(Y, theta, a, b, prior_mean, prior_prec, link = c("logit", "probit"), form = c("difficulty", "intercept")) {
  link <- match.arg(link); form <- match.arg(form); Y <- as_mat(Y); K <- ncol(Y)
  m <- if (is.matrix(prior_mean)) prior_mean else matrix(prior_mean, K, 2, byrow = TRUE)
  r <- .blk_items_mh(Y, as.numeric(theta), as.numeric(a), as.numeric(b), m, as.matrix(prior_prec), as.integer(link == "probit"), as.integer(form == "intercept"))
  structure(list(a = as.vector(r$a), b = as.vector(r$b)), accepted = as.logical(r$accepted))
}

#' @rdname compose
#' @export
draw_invwishart <- function(df, scale) solve(stats::rWishart(1, df, solve(scale))[, , 1])

#' @rdname compose
#' @export
draw_mean_mvn <- function(X, Sigma, mean0, cov0) {
  X <- as.matrix(X); Si <- solve(Sigma); C0 <- solve(cov0)
  draw_mvn(nrow(X) * Si + C0, drop(Si %*% colSums(X) + C0 %*% mean0))
}

#' @rdname compose
#' @export
draw_cor <- function(x1, x2, cor, steps = 5, step = 0.05, method = c("rw", "mala")) {
  method <- match.arg(method)
  n <- length(x1); s11 <- sum(x1^2); s22 <- sum(x2^2); s12 <- sum(x1 * x2)
  lp <- function(r) -n / 2 * log(1 - r^2) - (s11 - 2 * r * s12 + s22) / (2 * (1 - r^2)) + log(1 - r^2)   # + Jacobian of tanh
  acc <- 0
  if (method == "mala") {                                # MALA on z = atanh(cor)
    lpz <- function(z) { r <- tanh(z); u <- 1 - r^2
      dq <- (-s12 * u + r * (s11 - 2 * r * s12 + s22)) / u^2
      list(lp = lp(r), grad = ((n - 2) * r / u - dq) * u) }
    z <- atanh(cor)
    for (k in seq_len(steps)) { z <- draw_mala(z, lpz, step); acc <- acc + attr(z, "accepted") }
    return(structure(tanh(as.numeric(z)), accepted = acc / steps))
  }
  for (k in seq_len(steps)) {
    rp <- tanh(atanh(cor) + step * stats::rnorm(1))
    if (log(stats::runif(1)) < lp(rp) - lp(cor)) { cor <- rp; acc <- acc + 1 }
  }
  structure(cor, accepted = acc / steps)
}

#' @rdname compose
#' @export
calc_person_fit <- function(L, fitted, sigma2) {
  L <- as_mat(L); R <- (L - as_mat(fitted))^2 / matrix(sigma2, nrow(L), ncol(L), byrow = TRUE)
  data.frame(lt = rowSums(R, na.rm = TRUE), df = rowSums(!is.na(R)))
}
