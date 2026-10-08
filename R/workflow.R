# The workflow layer of the building blocks: a state list, a step function, and helpers that do the
# bookkeeping (loops, chains, storage, adaptation, priors, checks, packing for EM). Persons j are the
# rows and items i the columns of the matrices.

# ---- normal parts ------------------------------------------------------------------------------

#' Normal parts of full conditionals
#'
#' Most full conditionals of the samplers are normal. A *normal part* stores one as its precision
#' \eqn{Q} and its linear term \eqn{h} (`prec` and `num`; mean \eqn{Q^{-1}h}{h / Q}, variance
#' \eqn{Q^{-1}}{1 / Q}), so that the contributions of the data and of the prior simply add up:
#'
#' * [irt_part()] and [calc_normal_part()] return the data part of independent scalars (one per
#'   person or item);
#' * `reg_part(X, y, w)` returns the data part of a vector of regression coefficients \eqn{g} in
#'   \eqn{y = Xg + e}, \eqn{e_j \sim N(0, 1/w_j)}{e_j ~ N(0, 1 / w_j)} (\eqn{Q = X'WX}{Q = X'WX},
#'   \eqn{h = X'Wy}{h = X'Wy}), e.g. a latent regression of the abilities on covariates;
#' * `normal_prior(mean, sd)` turns a normal prior into a part (\eqn{Q = 1/sd^2}, \eqn{h = mean/sd^2});
#' * `part1 + part2` adds parts (elementwise, or a prior added to the diagonal of a regression part);
#' * `draw(part, lower, upper)` draws from it (truncated to `[lower, upper]` when given;
#'   a regression part is drawn with [draw_mvn()]);
#' * `print(part)` shows the conditional means and SDs.
#'
#' ```
#' theta <- draw(irt_part(Y, omega, a, b) + normal_prior(mean = X %*% g, sd = 1))
#' g <- draw(reg_part(X, theta) + normal_prior(0, 10))
#' ```
#' @param mean,sd Prior mean and SD (recycled).
#' @param prec,num Precision (a vector, or a matrix for a regression part) and linear term.
#' @param X Design matrix (one row per observation).
#' @param y Observations.
#' @param w Precisions of the observations (one value or one per observation).
#' @param part,x A normal part.
#' @param e1,e2 Normal parts.
#' @param lower,upper Truncation bounds.
#' @param ... Unused.
#' @return `normal_prior()`, `normal_part()`, `reg_part()` and `+`: a normal part; `draw()`: a vector of draws.
#' @examples
#' set.seed(1)
#' N <- 200; a <- c(1, 1.5, 2); b <- c(-1, 0, 1); th <- rnorm(N)
#' Y <- matrix(rbinom(N * 3, 1, plogis(sweep(outer(th, a), 2, a * b))), N)
#' om <- draw_omega(a, b, th)
#' p <- irt_part(Y, om, a, b) + normal_prior(0, 1)        # data part + prior N(0, 1)
#' p
#' theta <- draw(p)
#' draw(normal_prior(1, 1), lower = 0)                     # N+(1, 1)
#' X <- cbind(1, rnorm(N))
#' draw(reg_part(X, theta) + normal_prior(0, 10))          # coefficients of theta on X
#' @name normal_parts
NULL

#' @rdname normal_parts
#' @export
normal_part <- function(prec, num) {
  num <- as.numeric(num)
  if (is.matrix(prec) && ncol(prec) == 1) prec <- as.numeric(prec)   # a column of precisions (from C++)
  if (is.matrix(prec)) {
    if (nrow(prec) != ncol(prec) || nrow(prec) != length(num))
      stop(sprintf("a regression part needs a p x p precision and p linear terms (got %d x %d and %d)", nrow(prec), ncol(prec), length(num)), call. = FALSE)
    storage.mode(prec) <- "double"
  } else {
    prec <- as.numeric(prec)
    if (length(prec) != length(num) && length(prec) != 1 && length(num) != 1)
      stop(sprintf("prec (length %d) and num (length %d) must have the same length", length(prec), length(num)), call. = FALSE)
  }
  structure(list(prec = prec, num = num), class = "normal_part")
}

#' @rdname normal_parts
#' @export
normal_prior <- function(mean = 0, sd = 1) {
  if (any(sd <= 0)) stop("normal_prior(): sd must be positive", call. = FALSE)
  normal_part(1 / as.numeric(sd)^2, as.numeric(mean) / as.numeric(sd)^2)
}

#' @rdname normal_parts
#' @export
reg_part <- function(X, y, w = 1) {
  X <- as.matrix(X); y <- as.numeric(y)
  if (nrow(X) != length(y)) stop(sprintf("reg_part(): X has %d rows but y has %d values", nrow(X), length(y)), call. = FALSE)
  if (!length(w) %in% c(1, length(y))) stop("reg_part(): w must be one value or one per observation", call. = FALSE)
  ok <- !is.na(y); X <- X[ok, , drop = FALSE]; y <- y[ok]; w <- rep_len(w, length(ok))[ok]
  normal_part(crossprod(X, w * X), crossprod(X, w * y))
}

#' @rdname normal_parts
#' @export
`+.normal_part` <- function(e1, e2) {
  if (!inherits(e1, "normal_part") || !inherits(e2, "normal_part"))
    stop("only normal parts can be added: use normal_prior(mean, sd) for a prior", call. = FALSE)
  if (is.matrix(e1$prec) || is.matrix(e2$prec)) {                # regression part: elementwise parts go on the diagonal
    p <- if (is.matrix(e1$prec)) nrow(e1$prec) else nrow(e2$prec)
    as_mv <- function(e) if (is.matrix(e$prec)) e else {
      if (!length(e$prec) %in% c(1, p) || !length(e$num) %in% c(1, p))
        stop(sprintf("a part of length %d cannot be added to a regression part of %d coefficients", max(length(e$prec), length(e$num)), p), call. = FALSE)
      list(prec = diag(rep_len(e$prec, p), p), num = rep_len(e$num, p))
    }
    e1 <- as_mv(e1); e2 <- as_mv(e2)
    if (nrow(e1$prec) != nrow(e2$prec)) stop("regression parts of different dimensions cannot be added", call. = FALSE)
    return(normal_part(e1$prec + e2$prec, e1$num + e2$num))
  }
  n <- max(length(e1$prec), length(e2$prec))
  if (!(length(e1$prec) %in% c(1, n) && length(e2$prec) %in% c(1, n)))
    stop(sprintf("normal parts of lengths %d and %d cannot be added", length(e1$prec), length(e2$prec)), call. = FALSE)
  normal_part(e1$prec + e2$prec, e1$num + e2$num)
}

#' @rdname normal_parts
#' @export
draw <- function(part, lower = -Inf, upper = Inf) {
  if (!inherits(part, "normal_part")) stop("draw() needs a normal part (irt_part(), calc_normal_part(), reg_part(), normal_prior(), ...)", call. = FALSE)
  if (is.matrix(part$prec)) {
    ev <- eigen(part$prec, symmetric = TRUE, only.values = TRUE)$values
    if (!all(ev > 0)) stop("draw(): the precision matrix must be positive definite (did you forget the prior?)", call. = FALSE)
    return(draw_mvn(part$prec, part$num, lower, upper))
  }
  n <- max(length(part$prec), length(part$num)); Q <- rep_len(part$prec, n); h <- rep_len(part$num, n)
  if (any(!(Q > 0))) stop("draw(): the precision must be positive (did you forget the prior?)", call. = FALSE)
  if (all(is.infinite(lower)) && all(is.infinite(upper))) draw_gauss(Q, h)
  else { lo <- rep_len(lower, n); hi <- rep_len(upper, n); vapply(seq_len(n), function(i) draw_tnorm(h[i] / Q[i], 1 / sqrt(Q[i]), lo[i], hi[i]), 0) }
}

#' @rdname normal_parts
#' @export
print.normal_part <- function(x, ...) {
  if (is.matrix(x$prec)) {
    V <- solve(x$prec); m <- drop(V %*% x$num); s <- sqrt(diag(V))
    cat(sprintf("regression part of %d coefficients: precision matrix and linear term; mean = solve(prec, num)\n", length(m)))
  } else {
    m <- x$num / x$prec; s <- 1 / sqrt(x$prec)
    cat(sprintf("normal part of %d value(s): precision and linear term; mean = num / prec, sd = 1 / sqrt(prec)\n", length(m)))
  }
  print(utils::head(data.frame(mean = m, sd = s), 6), digits = 3)
  if (length(m) > 6) cat("  ...\n")
  invisible(x)
}

# ---- run_sampler -----------------------------------------------------------------------------------

#' Run a sampler written from the building blocks
#'
#' The loop around a sampler: chains, burn-in, thinning, storage, running means of person
#' parameters and the summary table, as in JAGS. The sampler is one function `step(s, data)` that
#' takes the current state `s` (a named list of everything the sampler updates) and returns the new
#' state, as in a JAGS model block written out by hand.
#'
#' * `init(data)` returns the starting state (called once per chain, after `set.seed()`, so random
#'   starting values differ between chains).
#' * `step(s, data)` updates every block once and returns `s`. During burn-in `s$.burn` is `TRUE`
#'   and `s$.iter` is the iteration, for step-size adaptation ([mcmc_adapt()]).
#' * `monitor`: the names of the elements of `s` to store (vectors are stored as `a[1]`, `a[2]`, ...).
#' * `keep_mean`: the names of elements (e.g. person parameters) of which only the posterior mean and
#'   SD are kept.
#' @param data A list of data (passed to `init()` and `step()`).
#' @param init Function of `data` returning the starting state (or a list, used by every chain).
#' @param step Function of `(s, data)` returning the new state.
#' @param monitor Names of the elements of the state to store.
#' @param n_iter,n_burnin,n_thin Iterations per chain, burn-in and thinning (as [gibbs()]).
#' @param n_chain Number of chains.
#' @param seed Seed of the first chain (chain c uses `seed + c - 1`).
#' @param keep_mean Names of elements whose posterior mean and SD are kept.
#' @param object,x A `block_fit`.
#' @param digits Number of decimals.
#' @param ... Unused.
#' @return A `block_fit`: `draws` (a `coda::mcmc.list`), `mean` (posterior means and SDs of
#'   `keep_mean`), `last` (the last state of every chain) and `time`. `summary()` gives the table of
#'   posterior means, SDs, 95% intervals, R-hat and effective sample sizes; `as.matrix()` the draws.
#'   `print()` warns when an R-hat exceeds 1.05 or an effective sample size is below 100, and names
#'   the usual cures.
#' @seealso [check_sampler()] to check that the sampler draws from the right posterior, [run_em()].
#' @examples
#' set.seed(1)
#' N <- 300; K <- 6; a0 <- runif(K, 0.8, 2); b0 <- rnorm(K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(rnorm(N), a0), 2, a0 * b0))), N)
#' fit <- run_sampler(
#'   data = list(Y = Y),
#'   init = function(d) list(theta = rnorm(nrow(d$Y)), a = rep(1, ncol(d$Y)), b = rep(0, ncol(d$Y))),
#'   step = function(s, d) {
#'     s$omega <- draw_omega(s$a, s$b, s$theta)                     # Polya-Gamma
#'     s[c("a", "b")] <- draw_items_pg(d$Y, s$omega, s$theta, s$a)  # b, then a
#'     s$theta <- draw(irt_part(d$Y, s$omega, s$a, s$b) + normal_prior(0, 1))
#'     s
#'   },
#'   monitor = c("a", "b"), keep_mean = "theta", n_iter = 1000)
#' fit
#' @export
run_sampler <- function(data, init, step, monitor, n_iter = 2000, n_burnin = floor(n_iter / 2), n_chain = 2, n_thin = 1,
                        seed = 1, keep_mean = NULL) {
  if (!is.function(step)) stop("step must be a function(s, data) returning the new state", call. = FALSE)
  if (n_burnin >= n_iter) stop("n_burnin must be smaller than n_iter", call. = FALSE)
  t0 <- Sys.time()
  one <- function(ch) {
    set.seed(seed + ch - 1)
    s <- if (is.function(init)) init(data) else init
    if (!is.list(s)) stop("init must return the state: a named list", call. = FALSE)
    miss <- setdiff(c(monitor, keep_mean), names(s))
    if (length(miss)) stop("not in the starting state: ", paste(miss, collapse = ", "), call. = FALSE)
    out <- NULL; r <- 0; acc <- list()
    for (it in seq_len(n_iter)) {
      s$.iter <- it; s$.burn <- it <= n_burnin
      s <- step(s, data)
      if (!is.list(s)) stop(sprintf("step() must return the state list s (iteration %d)", it), call. = FALSE)
      if (it > n_burnin && (it - n_burnin) %% n_thin == 0) {
        v <- flatten_state(s, monitor)
        if (anyNA(v)) stop(sprintf("NA in %s at iteration %d", names(v)[is.na(v)][1], it), call. = FALSE)
        if (is.null(out)) out <- matrix(NA_real_, (n_iter - n_burnin) %/% n_thin, length(v), dimnames = list(NULL, names(v)))
        r <- r + 1; out[r, ] <- v
        for (k in keep_mean) { x <- as.numeric(s[[k]]); acc[[k]] <- if (is.null(acc[[k]])) list(s1 = x, s2 = x^2) else list(s1 = acc[[k]]$s1 + x, s2 = acc[[k]]$s2 + x^2) }
      }
    }
    s$.iter <- s$.burn <- NULL
    list(draws = coda::mcmc(out[seq_len(r), , drop = FALSE], start = n_burnin + n_thin, thin = n_thin), acc = acc, n = r, last = s)
  }
  ch <- lapply(seq_len(n_chain), one)
  mean <- lapply(stats::setNames(keep_mean, keep_mean), function(k) {
    n <- sum(vapply(ch, `[[`, 0, "n")); s1 <- Reduce(`+`, lapply(ch, function(c) c$acc[[k]]$s1)); s2 <- Reduce(`+`, lapply(ch, function(c) c$acc[[k]]$s2))
    data.frame(mean = s1 / n, sd = sqrt(pmax(s2 / n - (s1 / n)^2, 0)))
  })
  structure(list(draws = coda::mcmc.list(lapply(ch, `[[`, "draws")), mean = mean, last = lapply(ch, `[[`, "last"),
                 settings = list(n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin, n_chain = n_chain, seed = seed),
                 time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), class = "block_fit")
}

flatten_state <- function(s, names) {
  unlist(lapply(names, function(n) {
    x <- s[[n]]
    if (!is.numeric(x) && !is.logical(x)) stop(sprintf("monitor: %s is not numeric", n), call. = FALSE)
    nm <- if (length(x) == 1) n else if (is.matrix(x)) sprintf("%s[%d,%d]", n, row(x), col(x)) else sprintf("%s[%d]", n, seq_along(x))
    stats::setNames(as.numeric(x), nm)
  }))
}

#' @rdname run_sampler
#' @export
summary.block_fit <- function(object, ...) {
  M <- as.matrix(object$draws)
  rh <- if (coda::nchain(object$draws) > 1) tryCatch(coda::gelman.diag(object$draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1], error = function(e) rep(NA, ncol(M))) else rep(NA, ncol(M))
  data.frame(parameter = colnames(M), mean = colMeans(M), sd = apply(M, 2, stats::sd),
             q2.5 = apply(M, 2, stats::quantile, 0.025), q97.5 = apply(M, 2, stats::quantile, 0.975),
             rhat = unname(rh), ess = round(unname(coda::effectiveSize(object$draws))), row.names = NULL)
}

#' @rdname run_sampler
#' @export
print.block_fit <- function(x, digits = 3, ...) {
  s <- summary(x); st <- x$settings
  cat(sprintf("%s: %d chains x %d iterations (burn-in %d, thin %d), %.1f s\n",
              if (inherits(x, "block_model_fit")) "Gibbs for the declared model" else "Sampler from the building blocks",
              st$n_chain, st$n_iter, st$n_burnin, st$n_thin, x$time))
  print(format(utils::head(s, 20), digits = digits), row.names = FALSE)
  if (nrow(s) > 20) cat(sprintf("  ... %d more (summary(fit))\n", nrow(s) - 20))
  cat(sprintf("max R-hat %.3f, min ESS %d%s\n", max(s$rhat, -Inf, na.rm = TRUE), min(s$ess),
              if (length(x$mean)) sprintf("; posterior means kept for %s (fit$mean)", paste(names(x$mean), collapse = ", ")) else ""))
  mixing_note(s, declared = inherits(x, "block_model_fit"))
  invisible(x)
}

# a note when the chains have not mixed: worst parameters and the usual cures
mixing_note <- function(s, rhat = 1.05, ess = 100, declared = FALSE) {
  bad_r <- which(s$rhat > rhat); bad_e <- which(s$ess < ess)
  if (!length(bad_r) && !length(bad_e)) return(invisible())
  worst <- function(i, v, dec, fmt) { o <- utils::head(i[order(v[i], decreasing = dec)], 3); paste(sprintf(fmt, s$parameter[o], v[o]), collapse = ", ") }
  cat("Warning: the chains have not mixed:",
      paste(c(if (length(bad_r)) sprintf("R-hat > %g for %d parameter(s) (%s)", rhat, length(bad_r), worst(bad_r, s$rhat, TRUE, "%s %.2f")),
              if (length(bad_e)) sprintf("ESS < %d for %d parameter(s) (%s)", ess, length(bad_e), worst(bad_e, s$ess, FALSE, "%s %.0f"))), collapse = "; "), "\n")
  if (declared) cat("  Run longer, or change the steps: a direction identified only by the prior (e.g. theta + c, b + c) -> step_shift();",
                    "very easy or difficult items -> step_mala(a, b); chains in different modes -> data-based init.\n")
  else cat("  Run longer, or change the sampler: a direction identified only by the prior (e.g. theta + c, b + c, g0 + c)",
           "-> draw_shift(); very easy or difficult items -> draw_items_joint() before the Polya-Gamma steps.\n")
}

#' @rdname run_sampler
#' @export
as.matrix.block_fit <- function(x, ...) as.matrix(x$draws)

# ---- check_step: Geweke's test for a sampler written by the user (check_sampler(step = ...)) ----------
# The successive-conditional chain (step, then simulate the data) is compared with direct prior draws.

check_step <- function(step, prior, simulate, monitor, n_iter = 2e4, n_prior = 4000, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), seed = 1) {
  set.seed(seed)
  s <- prior(); D <- simulate(s)
  G <- NULL
  for (it in seq_len(n_iter)) {
    s <- step(s, D)
    v <- flatten_state(s, monitor)
    if (is.null(G)) G <- matrix(NA_real_, n_iter, length(v), dimnames = list(NULL, names(v)))
    G[it, ] <- v
    D <- simulate(s)
  }
  G <- G[-seq_len(round(n_iter / 10)), , drop = FALSE]
  R <- t(vapply(seq_len(n_prior), function(i) flatten_state(prior(), monitor), numeric(ncol(G))))
  colnames(R) <- colnames(G)
  geweke_table(G, R, probs, n_prior,
               sprintf("%d iterations (first 10%% dropped); P(parameter <= prior quantile) in the chain vs %g prior draws", n_iter, n_prior))
}

# ---- draw_shift ---------------------------------------------------------------------------------------

#' Exact shift of parameters that move together
#'
#' Some directions of a model are identified only by the priors: adding the same constant to the
#' abilities, the difficulties and the intercept of a latent regression leaves the likelihood of
#' a 2PL unchanged. A Gibbs sampler moves along such a direction very slowly. `draw_shift()` draws
#' the shift \eqn{c} of the move \eqn{x_k \leftarrow x_k + c\, d_k}{x_k <- x_k + c d_k} exactly from its conditional
#' (a generalized Gibbs step; Liu & Sabatti, 2000): only the normal priors of the shifted
#' parameters depend on \eqn{c}. Parameters listed in `values` without a prior are shifted with
#' a prior that moves along (e.g. abilities whose prior mean contains the shifted intercept).
#'
#' The caller states that the likelihood is unchanged by the shift. Common cases:
#'
#' | model | shift | `values`, `prior` |
#' |---|---|---|
#' | 2PL, theta ~ N(g0 + x'g, 1) | theta + c, b + c, g0 + c | `list(theta, b, g0)`, priors of `b` and `g0` |
#' | log T = lambda - zeta | lambda + c, zeta + c | `list(lambda, zeta)`, priors of both (speed: N(m_j, 1 / w_j)), `lower = list(lambda = 0)` for a truncated prior |
#' | intercept form a theta - d | theta + c, d + a c | `along = list(d = a)` |
#' @param values Named list of the current values to shift.
#' @param prior Named list of normal priors `c(mean, sd)` (or `list(mean = , sd = )` with one value
#'   per element) for the shifted parameters whose prior does not move along.
#' @param along Named list of direction multipliers (default 1).
#' @param lower Named list of lower bounds of truncated priors (e.g. positive loadings): the shift
#'   is drawn from its conditional truncated so that every bound holds.
#' @return `values` shifted, with attribute `shift`.
#' @references Liu, J. S., & Sabatti, C. (2000). Generalised Gibbs sampler and multigrid Monte Carlo
#'   for Bayesian computation. *Biometrika, 87*, 353-369.
#' @examples
#' s <- draw_shift(list(theta = rnorm(5), b = c(-1, 0, 1), g0 = 0.2),
#'                 prior = list(b = c(0, 3), g0 = c(0, 10)))
#' attr(s, "shift")
#' @export
draw_shift <- function(values, prior, along = list(), lower = list()) {
  if (!length(prior)) stop("draw_shift() needs the prior of at least one shifted parameter", call. = FALSE)
  bad <- setdiff(c(names(prior), names(along), names(lower)), names(values))
  if (length(bad)) stop("not in values: ", paste(bad, collapse = ", "), call. = FALSE)
  Q <- 0; h <- 0
  for (k in names(prior)) {
    x <- values[[k]]; d <- rep_len(along[[k]] %||% 1, length(x)); p <- prior[[k]]
    m <- rep_len(if (is.list(p)) p$mean else p[1], length(x)); sd <- rep_len(if (is.list(p)) p$sd else p[2], length(x))
    Q <- Q + sum(d^2 / sd^2); h <- h + sum(d * (m - x) / sd^2)
  }
  lo <- -Inf; hi <- Inf                                                # truncated priors: x + c d >= lower
  for (k in names(lower)) {
    x <- values[[k]]; d <- rep_len(along[[k]] %||% 1, length(x)); r <- (rep_len(lower[[k]], length(x)) - x) / d
    if (any(d > 0)) lo <- max(lo, r[d > 0]); if (any(d < 0)) hi <- min(hi, r[d < 0])
  }
  c0 <- if (is.finite(lo) || is.finite(hi)) draw_tnorm(h / Q, 1 / sqrt(Q), lo, hi) else draw_gauss(Q, h)
  for (k in names(values)) values[[k]] <- values[[k]] + c0 * (along[[k]] %||% 1)
  structure(values, shift = c0)
}

# ---- EM: run_em, estep_irt, mstep_items ---------------------------------------------------------------

#' EM algorithms from components, with the bookkeeping done
#'
#' The EM counterpart of [run_sampler()]: the algorithm is one function `step(s, data)` that does one
#' E- and M-step from the current state `s` (a named list of parameters) and returns the new state
#' with the log-likelihood **at the input state** in `s$.loglik`. `run_em()` packs the state into
#' a vector (log scale for the `positive` parameters), accelerates with SQUAREM, checks convergence
#' and computes standard errors from the numerical Hessian of `loglik()`.
#'
#' For posterior modes (priors in the M-step), `s$.loglik` and `loglik()` are the log-likelihood
#' plus the log prior: the objective that the M-step increases.
#'
#' `run_em()` checks the step, because a wrong M-step still converges, to the wrong values. It warns
#' when the objective decreases between iterations (an EM step never lowers it), when `s$.loglik`
#' differs from `loglik()` at the starting values, and when the estimates are not a maximum of
#' `loglik()` (score statistic \eqn{g'(-H)^{-1}g > 0.01}{g' (-H)^-1 g > 0.01}).
#'
#' For item response models the E- and M-steps are ready-made:
#'
#' * `estep_irt(Y, a, b, mean, sd, nodes)`: the Gauss-Hermite E-step of a 2PL with
#'   \eqn{\theta_j \sim N(\text{mean}_j, \text{sd}^2)}{theta_j ~ N(mean_j, sd^2)} (missing responses skipped): posterior weights
#'   of the nodes `W`, the marginal log-likelihood `loglik`, the EAPs and posterior SDs, and the
#'   Pólya–Gamma sums for the M-step.
#' * `mstep_items(e, a, b, priors, b_mean, one_pl)`: the item parameters from the E-step (Pólya–Gamma
#'   minorizer: never lowers the objective); maximum likelihood, or the posterior mode with the item
#'   priors of `priors = rtirt_priors(a = , b = )` and, optionally, a prior mean of b per item
#'   (`b_mean`), as [draw_items_pg()].
#'
#' The components behind them are in [em_components].
#' @param data A list of data.
#' @param init Function of `data` returning the starting state (or the state).
#' @param step Function of `(s, data)`: one EM step; returns the new state with `.loglik`.
#' @param loglik Function of `(s, data)` giving the marginal log-likelihood (plus the log prior for
#'   posterior modes), for the standard errors and the checks; `NULL`: neither.
#' @param positive Names of the parameters kept positive (estimated on the log scale).
#' @param se Compute standard errors.
#' @param max_iter,tol,accelerate As [em_squarem()].
#' @param score Optional function of `(s, data)` giving the gradient of `loglik()` as a named list
#'   (one vector per parameter, natural scale), e.g. from Fisher's identity; the Hessian for the
#'   standard errors is then its numerical derivative (2p calls instead of about 4p^2 calls of
#'   `loglik()`).
#' @param hessian Optional function of `(s, data)` giving the Hessian of `loglik()` on the natural
#'   scale (a matrix, parameters in the order of the state), e.g. from Louis' formula; used instead of
#'   numerical differentiation (`NULL` result: numerical).
#' @param Y Responses (0/1, NA missing), N x K.
#' @param a,b Item parameters.
#' @param mean,sd Mean (one value or one per person) and SD of the ability distribution.
#' @param nodes Nodes from [calc_nodes()].
#' @param e Result of `estep_irt()`.
#' @param priors Item priors from [rtirt_priors()] (posterior mode); `NULL`: maximum likelihood.
#' @param b_mean Prior means of b, one per item (`NULL`: the mean of `priors$b`).
#' @param one_pl Keep a fixed.
#' @param object,x A `block_em`.
#' @param digits Number of decimals.
#' @param ... Unused.
#' @return `run_em()`: a `block_em` with `estimates` (estimate and SE), `state`, `loglik`,
#'   `iterations`, `converged`, `trace` and `vcov` (working scale). `estep_irt()`: a list;
#'   `mstep_items()`: `list(a, b)`.
#' @examples
#' set.seed(1)
#' N <- 1000; K <- 8; a0 <- runif(K, 0.8, 2); b0 <- rnorm(K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(rnorm(N), a0), 2, a0 * b0))), N)
#' fit <- run_em(
#'   data = list(Y = Y),
#'   init = function(d) list(a = rep(1, ncol(d$Y)), b = rep(0, ncol(d$Y))),
#'   step = function(s, d) {
#'     e <- estep_irt(d$Y, s$a, s$b)                 # E-step
#'     s[c("a", "b")] <- mstep_items(e, s$a, s$b)    # M-step
#'     s$.loglik <- e$loglik
#'     s
#'   },
#'   loglik = function(s, d) estep_irt(d$Y, s$a, s$b)$loglik,
#'   positive = "a")
#' fit
#' @export
run_em <- function(data, init, step, loglik = NULL, positive = character(), se = !is.null(loglik), max_iter = 1000, tol = 1e-10,
                   accelerate = TRUE, score = NULL, hessian = NULL) {
  s0 <- if (is.function(init)) init(data) else init
  par <- names(s0)[!startsWith(names(s0), ".") & vapply(s0, is.numeric, TRUE)]
  bad <- setdiff(positive, par); if (length(bad)) stop("positive: not in the state: ", paste(bad, collapse = ", "), call. = FALSE)
  if (any(unlist(s0[positive]) <= 0)) stop("the starting values of the positive parameters must be positive", call. = FALSE)
  skel <- s0[par]
  pack <- function(s) unlist(lapply(par, function(n) if (n %in% positive) log(s[[n]]) else s[[n]]))
  unpack <- function(x, s = s0) {
    v <- utils::relist(x, skel)
    for (n in par) s[[n]] <- if (n %in% positive) exp(v[[n]]) else v[[n]]
    s
  }
  Fmap <- function(x) {
    s <- step(unpack(x), data)
    if (is.null(s$.loglik)) stop("step() must set s$.loglik (the log-likelihood at the input state)", call. = FALSE)
    list(x = pack(s), obj = s$.loglik)
  }
  if (!is.null(loglik)) {                                         # .loglik must be loglik() at the input state
    l1 <- Fmap(pack(s0))$obj; l2 <- loglik(unpack(pack(s0)), data)
    if (abs(l1 - l2) > 1e-6 * max(1, abs(l2)))
      warning(sprintf("s$.loglik from step() (%.4f) differs from loglik() (%.4f) at the starting values: .loglik must be the objective at the input state",
                      l1, l2), call. = FALSE)
  }
  sq <- em_squarem(pack(s0), Fmap, max_iter = max_iter, tol = tol, accelerate = accelerate)
  if (all(sq$x == pack(s0)))
    warning("step() did not change any parameter: is there an M-step?", call. = FALSE)
  down <- which(diff(sq$trace) < -1e-8 * pmax(1, abs(sq$trace[-1])))
  if (length(down))
    warning(sprintf("the objective decreased at %d iteration(s) (first: iteration %d, by %.3g): an EM step never lowers it; check the M-step",
                    length(down), down[1], sq$trace[down[1]] - sq$trace[down[1] + 1]), call. = FALSE)
  st <- unpack(sq$x); st$.loglik <- NULL
  x <- sq$x; nm <- unlist(lapply(par, function(n) if (length(skel[[n]]) == 1) n else sprintf("%s[%d]", n, seq_along(skel[[n]]))))
  V <- NULL; sev <- rep(NA_real_, length(x))
  if (se && !is.null(loglik) && sq$converged) {
    f <- function(z) loglik(unpack(z), data)
    gsc <- if (!is.null(score)) function(z) {                     # the score on the packed scale (log scale for positive parameters)
      st <- unpack(z); g <- score(st, data)
      unlist(lapply(par, function(n) if (n %in% positive) g[[n]] * st[[n]] else g[[n]]))
    }
    Hx <- if (!is.null(hessian)) hessian(unpack(x), data)
    H <- if (!is.null(Hx)) {                                       # natural scale -> packed scale: H_z = D H_x D + diag(g_z) for log-scale parameters
      pos <- rep(par %in% positive, lengths(skel)); dz <- ifelse(pos, exp(x), 1)
      gz <- if (!is.null(gsc)) gsc(x) else vapply(seq_along(x), function(k) { h <- 1e-5 * max(1, abs(x[k])); e <- replace(numeric(length(x)), k, h); (f(x + e) - f(x - e)) / (2 * h) }, 0)
      Hz <- outer(dz, dz) * Hx; diag(Hz) <- diag(Hz) + ifelse(pos, gz, 0); (Hz + t(Hz)) / 2
    } else if (is.null(gsc)) stats::optimHess(x, f) else calc_hessian(gsc, x, h = 1e-5)
    V <- tryCatch(solve(-H), error = function(e) NULL)
    if (!is.null(V) && !all(eigen(-(H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values > 0)) {
      warning("the estimates are not a maximum of loglik(): its Hessian is not negative definite (does step() have an M-step for every parameter?)",
              call. = FALSE)
      V <- NULL
    } else if (!is.null(V)) {
      g <- if (!is.null(gsc)) gsc(x) else
        vapply(seq_along(x), function(k) { h <- 1e-5 * max(1, abs(x[k])); e <- replace(numeric(length(x)), k, h); (f(x + e) - f(x - e)) / (2 * h) }, 0)
      sc <- sum(g * (V %*% g))
      if (sq$converged && sc > 0.01)
        warning(sprintf("the estimates are not a maximum of loglik() (score statistic %.3g): check the M-step (for posterior modes, loglik() and .loglik include the log prior)", sc),
                call. = FALSE)
      sx <- sqrt(pmax(diag(V), 0))
      pos <- rep(par %in% positive, lengths(skel))
      sev <- ifelse(pos, exp(x) * sx, sx)
    } else warning("the Hessian is singular: no standard errors", call. = FALSE)
  }
  pos <- rep(par %in% positive, lengths(skel))
  est <- data.frame(parameter = nm, est = ifelse(pos, exp(x), x), se = sev, row.names = NULL)
  if (!sq$converged) warning(sprintf("EM did not converge (%d iterations)", sq$iterations), call. = FALSE)
  structure(list(estimates = est, state = st, loglik = sq$obj, iterations = sq$iterations, converged = sq$converged,
                 trace = sq$trace, vcov = V, has_loglik = !is.null(loglik)), class = "block_em")
}

#' @rdname run_em
#' @export
print.block_em <- function(x, digits = 3, ...) {
  cat(sprintf("EM from components: %s after %d iterations, log-likelihood %.4f\n",
              if (x$converged) "converged" else "NOT converged", x$iterations, x$loglik))
  print(format(utils::head(x$estimates, 30), digits = digits), row.names = FALSE)
  if (nrow(x$estimates) > 30) cat(sprintf("  ... %d more (fit$estimates)\n", nrow(x$estimates) - 30))
  if (all(is.na(x$estimates$se)))
    cat(if (!isTRUE(x$has_loglik)) "se = NA: give run_em(..., loglik = ) for standard errors\n"
        else if (!x$converged) "se = NA: EM did not converge\n" else "se = NA: the Hessian of loglik() is singular or not negative definite\n")
  invisible(x)
}

#' @rdname run_em
#' @export
summary.block_em <- function(object, ...) object$estimates

#' @rdname run_em
#' @export
estep_irt <- function(Y, a, b, mean = 0, sd = 1, nodes = calc_nodes(41)) {
  Y <- as.matrix(Y); N <- nrow(Y); K <- ncol(Y)
  if (length(a) != K || length(b) != K) stop(sprintf("a and b need one value per item (column of Y): %d items, got %d and %d", K, length(a), length(b)), call. = FALSE)
  common <- length(mean) == 1
  if (!common && length(mean) != N) stop("mean must be one value or one per person (row of Y)", call. = FALSE)
  TH <- if (common) mean + sd * nodes$z else outer(mean, sd * nodes$z, `+`)
  nw <- calc_node_weights(calc_loglik_nodes(Y, a, b, TH), nodes$logw)
  S <- calc_item_sums(Y, a, b, TH, nw$W)
  eap <- if (common) drop(nw$W %*% TH) else rowSums(nw$W * TH)
  e2 <- if (common) drop(nw$W %*% TH^2) else rowSums(nw$W * TH^2)
  kap <- ifelse(is.na(Y), 0, Y - 0.5)
  list(W = nw$W, loglik = sum(nw$ll), ll = nw$ll, eap = eap, psd = sqrt(pmax(e2 - eap^2, 0)),
       S = S, T0 = colSums(kap), T1 = colSums(kap * eap), nodes = TH)
}

#' @rdname run_em
#' @export
mstep_items <- function(e, a, b, priors = NULL, b_mean = NULL, one_pl = FALSE) {
  if (is.null(priors)) {
    if (!is.null(b_mean)) stop("b_mean is a prior mean: give priors = rtirt_priors() as well (posterior mode)", call. = FALSE)
    if (one_pl) return(list(a = a, b = (a * e$S$S1 - e$T0) / (a * e$S$S0)))
    return(update_items_logistic(e$S, e$T0, e$T1, a, b))
  }
  pv_of(priors)
  pb <- list(if (is.null(b_mean)) priors$b[["mean"]] else chk_b_mean(b_mean, length(a)), priors$b[["sd"]], priors$b[["lower"]], priors$b[["upper"]])
  if (one_pl) { vb <- pb[[2]]^2; return(list(a = a, b = pmin(pmax((a * (a * e$S$S1 - e$T0) + pb[[1]] / vb) / (1 / vb + a^2 * e$S$S0), pb[[3]]), pb[[4]]))) }
  update_items_logistic(e$S, e$T0, e$T1, a, b, prior_a = unname(priors$a), prior_b = pb)
}
