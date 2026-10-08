# MALA: a gradient-informed Metropolis-Hastings step for blocks without a conjugate full conditional.

#' Metropolis-adjusted Langevin steps
#'
#' A Metropolis–Hastings step whose proposal follows the gradient of the log density (MALA;
#' Roberts & Tweedie, 1996): \eqn{x^* = x + \frac{\epsilon^2}{2}\nabla\log\pi(x) + \epsilon z}{x* = x + eps^2 / 2 grad log pi(x) + eps z},
#' accepted with the usual MH probability (the proposal is not symmetric, so its density in both
#' directions enters the ratio). It replaces a random-walk Metropolis step where a block has no
#' conjugate full conditional: it moves in the direction of higher density and needs far fewer
#' rejections, especially in several dimensions (optimal acceptance about 0.574 against 0.234).
#'
#' * `draw_mala(x, logpost, step, joint)`: by default one step for many
#'   independent targets at once, one per row of `x` (or per element of a vector `x`, e.g. one per
#'   person); `logpost(x)` returns `list(lp = , grad = )`: the log density of every row (a vector)
#'   and its gradient (like `x`); `step` is one number or one per row. With `joint = TRUE`, `x`
#'   is one target and `lp` one number. Returns `x` with the attribute `accepted`.
#' * A parameter on a transformed scale (e.g. `u = log a`) needs the log Jacobian in `lp` (`+ u`)
#'   and in `grad` (`+ 1`); variables integrated out by the proposal (e.g. Pólya–Gamma omega) must be
#'   drawn anew before the next step uses them.
#' * `mcmc_adapt(step, accepted, iter, target)`: Robbins–Monro adaptation of the step size during
#'   burn-in (Andrieu & Thoms, 2008), \eqn{\log\epsilon \leftarrow \log\epsilon + t^{-0.6}(\alpha -
#'   \alpha^*)}{log eps <- log eps + t^-0.6 (acceptance - target)}; stop adapting after burn-in so that the chain is a proper Markov chain.
#'
#' [draw_cor()] uses it with `method = "mala"`.
#' @param x Current values: a vector (n scalar targets, or one target with `joint = TRUE`) or an
#'   n x d matrix (n independent targets of length d).
#' @param logpost Function of `x` returning `list(lp = , grad = )`.
#' @param step Step size(s) \eqn{\epsilon}{eps}.
#' @param joint `FALSE`: independent targets, one per row of `x`; `TRUE`: `x` is one target.
#' @param accepted Logical (or acceptance rate) of the last step(s).
#' @param iter Iteration number (for the decreasing adaptation rate).
#' @param target Target acceptance rate.
#' @return `draw_mala()`: the new `x` with attribute `accepted`; `mcmc_adapt()`: the new step(s).
#' @references Roberts, G. O., & Tweedie, R. L. (1996). Exponential convergence of Langevin
#'   distributions and their discrete approximations. *Bernoulli, 2*, 341-363.
#'
#'   Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC. *Statistics and Computing, 18*,
#'   343-373.
#' @examples
#' # abilities of a 2PL with known items: one MALA step per person, all persons at once
#' set.seed(1); N <- 200; a <- c(1, 1.5, 2); b <- c(-1, 0, 1)
#' theta0 <- rnorm(N); Y <- matrix(rbinom(N * 3, 1, plogis(sweep(outer(theta0, a), 2, a * b))), N)
#' lp <- function(th) { p <- plogis(sweep(outer(th, a), 2, a * b))
#'   list(lp = rowSums(Y * log(p) + (1 - Y) * log(1 - p)) - th^2 / 2,
#'        grad = drop((Y - p) %*% a) - th) }
#' th <- rep(0, N); eps <- rep(1, N)
#' for (t in 1:2000) {
#'   th <- draw_mala(th, lp, eps)
#'   if (t <= 500) eps <- mcmc_adapt(eps, attr(th, "accepted"), t)
#' }
#' mean(attr(th, "accepted"))
#' @name mala
NULL

#' @rdname mala
#' @export
draw_mala <- function(x, logpost, step, joint = FALSE) {
  shape <- x
  if (joint) x <- matrix(as.numeric(x), 1)              # one target of length(x)
  mat <- is.matrix(x); X <- if (mat) x else matrix(x, ncol = 1)
  n <- nrow(X); d <- ncol(X); eps <- rep_len(step, n)
  back <- function(Z) if (joint) { v <- shape; v[] <- as.numeric(Z); v } else if (mat) Z else drop(Z)
  ev <- function(z) {
    r <- logpost(back(z))
    if (!is.list(r) || is.null(r$lp) || is.null(r$grad))
      stop("draw_mala(): logpost(x) must return list(lp = log density, grad = gradient)", call. = FALSE)
    if (length(r$lp) != n)
      stop(sprintf("draw_mala(): logpost() returned %d log densities for %d target(s); %s", length(r$lp), n,
                   if (length(r$lp) == 1) "for one target of all elements use joint = TRUE" else "return one per row of x"), call. = FALSE)
    if (length(r$grad) != n * d) stop(sprintf("draw_mala(): grad has %d values, x has %d", length(r$grad), n * d), call. = FALSE)
    list(lp = r$lp, grad = matrix(r$grad, n, d))
  }
  cur <- ev(X)
  mu0 <- X + eps^2 / 2 * cur$grad
  Xp <- mu0 + eps * matrix(stats::rnorm(n * d), n, d)
  new <- ev(Xp)
  mu1 <- Xp + eps^2 / 2 * new$grad
  lq <- function(to, from_mu) -rowSums((to - from_mu)^2) / (2 * eps^2)       # log proposal density up to a constant
  logr <- new$lp - cur$lp + lq(X, mu1) - lq(Xp, mu0)
  acc <- is.finite(logr) & log(stats::runif(n)) < logr
  X[acc, ] <- Xp[acc, ]
  structure(back(X), accepted = acc)
}

#' @rdname mala
#' @export
mcmc_adapt <- function(step, accepted, iter, target = 0.574) exp(log(step) + iter^(-0.6) * (as.numeric(accepted) - target))
