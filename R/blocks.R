# Building blocks of the samplers, callable from R ---------------------------------------------

#' Building blocks for writing samplers in R
#'
#' Each `draw_*()` function draws one block from its full conditional and is named after the C++
#' function of the samplers that it calls (`irt_part()` draws nothing: it returns the accuracy part
#' of the conditional of ability). [algorithm()] names the function of every step of [gibbs()].
#'
#' The full conditionals used by [gibbs()], exposed so that new samplers (other blockings, cut or
#' semi-modular schemes, nested chains) can be put together in R while every draw still runs in C++.
#' All of them use R's random number generator, so `set.seed()` makes a composed sampler
#' reproducible, and a sampler written with them reproduces [gibbs()] step by step. Notation as in
#' the samplers: persons j = 1..N, items i = 1..K, `Y` the 0/1 responses, `L` the log response
#' times, logit P(y = 1) = a (theta - b), log T = lambda - zeta - rho theta + e.
#'
#' * `draw_gauss(prec, num)`: independent draws from N(num / prec, 1 / prec), elementwise (any
#'   conditionally normal scalar: ability, speed, rho_i, a shift).
#' * `draw_invgamma(shape, scale)`, `draw_tnorm(mean, sd, lower, upper)`: inverse gamma and
#'   truncated normal draws, elementwise.
#' * `draw_omega(a, b, theta)`: Pólya–Gamma variables PG(1, a_i (theta_j - b_i)), an N x K matrix.
#' * `draw_pg(eta)`: Pólya–Gamma variables PG(1, eta) for any matrix of linear predictors (models
#'   whose logit is not a (theta - b), e.g. several abilities or covariates in the item response).
#' * `draw_probit_z(Y, eta)`: for a probit response model, the latent normal responses of Albert & Chib
#'   (1993), N(eta, 1) truncated to the positive (y = 1) or negative (y = 0) half line; missing y
#'   (NA) are not truncated. Given them, the probit model is a normal linear model.
#' * `irt_part(Y, omega, a, b)`: the accuracy part of the conditional of
#'   ability given omega, a normal part (precision `prec` and linear term `num`, see [normal_parts]):
#'   `draw(irt_part(Y, omega, a, b) + normal_prior(0, 1))`.
#' * `draw_items_pg(Y, omega, theta, a, priors, one_pl, b_mean)`: b, then a, given omega (Gibbs);
#'   returns `list(a, b)`.
#' * `draw_items_joint(Y, theta, a, b, priors, b_mean)`: (a, b) jointly with omega integrated out
#'   (Metropolis–Hastings with a Fisher-scoring proposal in (a, -a b)); returns `list(a, b)` with the
#'   attribute `accepted` (one logical per item). It mixes much better than `draw_items_pg()` for very
#'   easy or very difficult items; draw omega anew after it.
#'
#' The item priors are those of [rtirt_priors()] (`a`: truncated normal on a > 0, `b`: normal);
#' `b_mean` gives b a prior mean per item (e.g. `W %*% delta` in an LLTM with item covariates `W`),
#' in place of the mean of `priors$b`.
#' * `draw_lambda(L, m, k1nu, W, mu, sd, lower)`: time intensities given the RT means
#'   `m` (N x K, e.g. -zeta - rho theta), the ALD shifts `k1nu` (zeros for normal times) and the
#'   precisions `W` (N x K); prior N(mu, sd^2) on (lower, Inf).
#' * `draw_sigma2(L, m, lambda, current, priors, nu = NULL, q = NULL)`: residual variances (normal) or
#'   ALD scales (with the mixing variables `nu` and the quantile `q`), from their `current` values
#'   (a half-t prior draws its inverse-gamma auxiliary variable from them).
#' * `draw_var(alpha, beta, current, prior)`: a variance with likelihood var^-alpha exp(-beta / var)
#'   (n residuals u ~ N(0, var): `alpha = n / 2`, `beta = sum(u^2) / 2`) and the prior `c(df, scale)` (half-t on its square root, two inverse-gamma steps) or `c(shape, scale)`
#'   (inverse gamma).
#' * `draw_nu(r, s, q)`: ALD mixing variables for residuals `r` and scales `s` (matrices).
#' * `draw_theta_collapsed(prec, num, W, y, rho, m0, g, tau)` and
#'   `draw_rt_items_collapsed(L, k1nu, W, theta, mz, tau, with_rho, mu, sd, rho_sd, lambda, rho)`:
#'   the collapsed steps of `gibbs(collapse = TRUE)` (speed integrated out).
#'
#' Missing responses (NA in `Y`) are skipped by `irt_part()`, `draw_items_pg()` and
#' `draw_items_joint()`. `draw_lambda()` and `draw_sigma2()` need complete log times (they stop on
#' NA); with missing times use [calc_normal_part()] with [draw()] and `draw_var()`. [run_sampler()] runs a sampler written from these blocks (chains, burn-in,
#' storage, summary) and [check_sampler()] checks it.
#'
#' These are low-level: arguments are checked only for their dimensions, and the caller is
#' responsible for a valid sampler (each block must be drawn from its conditional given the current
#' values of everything else, or from a marginal followed by the integrated variable).
#' `inst/examples/ex07_blocks.R` writes the `rtirt_cross()` sampler with them; `ex08_testlet.R`,
#' `ex09_mirt_rt.R` and `ex10_conditional_dependence.R` compose models from the literature.
#' @param prec,num,shape,scale,mean,sd,lower,upper Numeric vectors (see Details).
#' @param a,b,theta,lambda,rho,mz,m0,g,tau Numeric vectors.
#' @param Y,L,omega,m,k1nu,W,y,nu,r,s,eta Numeric matrices (persons x items).
#' @param priors Priors from [rtirt_priors()].
#' @param one_pl `TRUE`: keep a fixed.
#' @param b_mean Prior means of b, one per item (`NULL`: the mean of `priors$b`).
#' @param mu Prior mean (and `sd` the prior SD) of the time intensities; `lower` is their lower truncation.
#' @param current Current values (variances).
#' @param alpha,beta Exponents of the likelihood of a variance.
#' @param prior Prior of a variance: `c(df, scale)` or `c(shape, scale)`.
#' @param q Quantile level of the ALD.
#' @param with_rho Draw rho with lambda.
#' @param rho_sd Prior SD of rho.
#' @return Draws (vectors, matrices or lists, see Details).
#' @examples
#' set.seed(1)
#' draw_gauss(prec = c(1, 4), num = c(0, 2))          # N(0, 1) and N(0.5, 1/4)
#' theta <- rnorm(200); Y <- matrix(rbinom(600, 1, plogis(theta)), 200)
#' om <- draw_omega(a = rep(1, 3), b = rep(0, 3), theta = theta)
#' str(draw_items_pg(Y, om, theta, a = rep(1, 3), priors = rtirt_priors()))
#' @name blocks
NULL

pv_of <- function(priors) {
  if (!inherits(priors, "rtirt_priors")) stop("priors must come from rtirt_priors()", call. = FALSE)
  prior_vector(priors)
}
chk_b_mean <- function(b_mean, K) {
  if (is.null(b_mean)) return(numeric())
  if (length(b_mean) != K || anyNA(b_mean)) stop(sprintf("b_mean needs one prior mean per item (%d), got %d", K, length(b_mean)), call. = FALSE)
  as.numeric(b_mean)
}
as_mat <- function(x) { x <- as.matrix(x); storage.mode(x) <- "double"; x }

# dimension checks with readable messages; missing responses (NA) become y = 1/2 with omega = 0, so
# that the cell adds nothing to the Polya-Gamma terms (kappa = y - 1/2 = 0)
chk_irt <- function(Y, omega, a = NULL, b = NULL, theta = NULL, what) {
  Y <- as_mat(Y); N <- nrow(Y); K <- ncol(Y)
  err <- function(...) stop(sprintf(...), call. = FALSE)
  if (!is.null(a) && length(a) != K) err("%s(): a has %d values but Y has %d items (columns)", what, length(a), K)
  if (!is.null(b) && length(b) != K) err("%s(): b has %d values but Y has %d items (columns)", what, length(b), K)
  if (!is.null(theta) && length(theta) != N) err("%s(): theta has %d values but Y has %d persons (rows)", what, length(theta), N)
  if (!is.null(omega)) {
    omega <- as_mat(omega)
    if (!identical(dim(omega), dim(Y))) err("%s(): omega is %d x %d but Y is %d x %d (persons x items)", what, nrow(omega), ncol(omega), N, K)
    if (anyNA(Y)) { m <- is.na(Y); Y[m] <- 0.5; omega[m] <- 0 }
  }
  list(Y = Y, omega = omega)
}

#' @rdname blocks
#' @export
draw_gauss <- function(prec, num) {
  if (any(!(prec > 0), na.rm = TRUE)) stop("draw_gauss(): the precision must be positive", call. = FALSE)
  as.vector(.blk_gauss(as.numeric(prec), as.numeric(num)))
}
#' @rdname blocks
#' @export
draw_invgamma <- function(shape, scale) { n <- max(length(shape), length(scale)); as.vector(.blk_invgamma(rep_len(shape, n), rep_len(scale, n))) }
#' @rdname blocks
#' @export
draw_tnorm <- function(mean, sd, lower = -Inf, upper = Inf) { n <- max(length(mean), length(sd)); as.vector(.blk_tnorm(rep_len(mean, n), rep_len(sd, n), lower, upper)) }
#' @rdname blocks
#' @export
draw_omega <- function(a, b, theta) {
  if (length(a) != length(b))
    stop(sprintf("draw_omega(a, b, theta): a and b need one value per item (lengths %d and %d)%s", length(a), length(b),
                 if (length(theta) == length(b)) "; the order is draw_omega(a, b, theta): items first, persons last" else ""), call. = FALSE)
  .blk_omega(as.numeric(a), as.numeric(b), as.numeric(theta))
}
#' @rdname blocks
#' @export
draw_probit_z <- function(Y, eta) .blk_probit_z(as_mat(Y), as_mat(eta))
#' @rdname blocks
#' @export
draw_pg <- function(eta) .blk_pg(as_mat(eta))
#' @rdname blocks
#' @export
irt_part <- function(Y, omega, a, b) {
  d <- chk_irt(Y, omega, a = a, b = b, what = "irt_part")
  r <- .blk_irt_part(d$Y, d$omega, as.numeric(a), as.numeric(b)); normal_part(r$prec, r$num)
}
#' @rdname blocks
#' @export
draw_items_pg <- function(Y, omega, theta, a, priors = rtirt_priors(), one_pl = FALSE, b_mean = NULL) {
  d <- chk_irt(Y, omega, a = a, theta = theta, what = "draw_items_pg")
  r <- .blk_items_pg(d$Y, d$omega, as.numeric(theta), as.numeric(a), pv_of(priors), one_pl, chk_b_mean(b_mean, length(a)))
  lapply(r, as.vector)
}
#' @rdname blocks
#' @export
draw_items_joint <- function(Y, theta, a, b, priors = rtirt_priors(), b_mean = NULL) {
  chk_irt(Y, NULL, a = a, b = b, theta = theta, what = "draw_items_joint")
  r <- .blk_items_joint(as_mat(Y), as.numeric(theta), as.numeric(a), as.numeric(b), pv_of(priors), chk_b_mean(b_mean, length(a)))
  an <- as.vector(r$a); bn <- as.vector(r$b)
  structure(list(a = an, b = bn), accepted = an != a | bn != b)
}
#' @rdname blocks
#' @export
draw_lambda <- function(L, m, k1nu, W, mu, sd, lower = 0) {
  chk_no_na(L, "draw_lambda")
  as.vector(.blk_lambda(as_mat(L), as_mat(m), as_mat(k1nu), as_mat(W), mu, sd, lower))
}
# the RT blocks of the built-in samplers need complete log times
chk_no_na <- function(L, what) if (anyNA(L))
  stop(sprintf("%s(): the log times contain NA; for missing times use calc_normal_part() with draw() and draw_var()", what), call. = FALSE)
#' @rdname blocks
#' @export
draw_sigma2 <- function(L, m, lambda, current, priors = rtirt_priors(), nu = NULL, q = NULL) {
  if (!is.null(nu) && is.null(q)) stop("an ALD scale needs q", call. = FALSE)
  chk_no_na(L, "draw_sigma2")
  as.vector(.blk_sigma2(as_mat(L), as_mat(m), as.numeric(lambda), pv_of(priors), if (!is.null(nu)) as_mat(nu), if (is.null(q)) 0.5 else q,
                        rep_len(as.numeric(current), ncol(as.matrix(L)))))
}
#' @rdname blocks
#' @export
draw_var <- function(alpha, beta, current, prior = c(df = 3, scale = 1)) {
  fam <- if (identical(sort(names(prior)), c("df", "scale"))) 1 else if (identical(sort(names(prior)), c("scale", "shape"))) 0 else
    stop("prior must be c(df = , scale = ) or c(shape = , scale = )", call. = FALSE)
  p <- if (fam == 1) prior[c("df", "scale")] else prior[c("shape", "scale")]
  .blk_var(alpha, beta, fam, p[[1]], p[[2]], current)
}
#' @rdname blocks
#' @export
draw_nu <- function(r, s, q) .blk_nu(as_mat(r), as_mat(s), q)
#' @rdname blocks
#' @export
draw_theta_collapsed <- function(prec, num, W, y, rho, m0, g, tau)
  as.vector(.blk_theta_collapsed(as.numeric(prec), as.numeric(num), as_mat(W), as_mat(y), as.numeric(rho), as.numeric(m0), as.numeric(g), as.numeric(tau)))
#' @rdname blocks
#' @export
draw_rt_items_collapsed <- function(L, k1nu, W, theta, mz, tau, with_rho, mu, sd, rho_sd, lambda, rho, lower = 0) {
  r <- .blk_rt_items_collapsed(as_mat(L), as_mat(k1nu), as_mat(W), as.numeric(theta), as.numeric(mz), as.numeric(tau),
                               with_rho, mu, sd, rho_sd, as.numeric(lambda), as.numeric(rho), lower)
  lapply(r, as.vector)
}
