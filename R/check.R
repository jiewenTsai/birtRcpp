# check_sampler(): does a Gibbs sampler draw from the right posterior? ------------------------------
#
# Two simulation checks of the C++ samplers of gibbs(), both based on the joint distribution of
# parameters and data under a proper, data-independent prior:
#   * Geweke (2004): alternate one sweep of the sampler given the data with a new data set drawn
#     given the current parameters (the successive-conditional simulator). The parameter chain then
#     has the prior as its stationary distribution; each parameter is compared with prior draws
#     through P(parameter <= prior quantile), with Monte Carlo SEs from the effective sample size.
#   * Simulation-based calibration (Cook, Gelman & Rubin, 2006; Talts et al., 2018): draw parameters
#     from the prior, simulate data, run the sampler; the rank of the true value among thinned
#     posterior draws is uniform when the sampler is correct (chi-square test over bins).

#' Check a Gibbs sampler by simulation
#'
#' Tests whether the sampler of [gibbs()] for a model draws from the correct posterior, with two
#' general methods that need only the prior and the likelihood:
#'
#' * `method = "geweke"` (Geweke, 2004): the successive-conditional simulator. One sweep of the
#'   sampler given the data alternates with new person values and a new data set simulated given the
#'   current parameters. (Simulating new data from the old person values is also exact, but the
#'   persons are then held by their own data and the chain crawls: a false alarm in practice.) If the sampler is correct, the chain of the parameters has the **prior** as
#'   its distribution. For every parameter and each probability in `probs`, the share of the chain
#'   below the prior quantile is compared with the probability (z-test with the effective sample
#'   size of the chain); p-values are Holm-adjusted over all comparisons. Cheap: one sweep per
#'   iteration on a small data set.
#' * `method = "sbc"`, simulation-based calibration (Cook, Gelman & Rubin, 2006; Talts et al.,
#'   2018; Modrák et al., 2025): `n_rep` times, draw the parameters from the prior, simulate
#'   data, run the sampler and record the rank of the true value among `n_draws` thinned posterior
#'   draws. The ranks are uniform for a correct sampler; a chi-square test over 10 bins per parameter
#'   (Holm-adjusted).
#'
#' Both use the model only as a template: its type, quantile, settings, number of items and
#' (the first `n_persons` rows of) its covariates. Data are simulated, so a small `n_persons` keeps
#' the posterior wide and the check fast. The priors must be proper and must not depend on the data
#' (`rtirt_priors_julia()` centres lambda on the data and cannot be used).
#'
#' `sampler_priors` lets the sampler use other priors than the simulation: a check of the check,
#' which must then report a failure.
#'
#' **Chain length matters for both.** The successive-conditional chain must travel over the whole
#' prior, and its parameters move by about one posterior width per iteration, so parameters with wide priors
#' (lambda ~ N(0, 10^2)) need long chains; the print method warns when an effective sample size is
#' small. In SBC the posterior draws of each fit must be close to independent after thinning: with
#' short fits a slowly mixing parameter (such as rho) piles its ranks up at one end, which looks
#' like an error of the sampler. Both effects disappear with longer chains, a real error does not.
#' **Samplers written by the user** from the building blocks ([run_sampler()]'s `step()`) are checked
#' with `step`, `prior`, `simulate` and `monitor` (Geweke test): `prior()` returns a state drawn from
#' the prior (all parameters and the person parameters), `simulate(s)` the data list given the
#' state, and `step(s, data)` is the sampler. Keep the model small (10-20 persons, a few items) and
#' the priors proper and not too wide, so that the chain mixes. The test finds wrong full
#' conditionals and forgotten or wrong prior terms, which summaries of a fit cannot show; its power
#' is that of a test on a small model (an error that moves the posterior by a fraction of a
#' posterior SD can pass), so passing is evidence, not proof.
#' @param model A model object (see [models]); sampled or not. Or a [block_model()]: its stacked steps
#'   (`steps`, see [gibbs()]) are checked by the Geweke test, with the parameters drawn from the
#'   declared priors and the data from the declared model (use a small data set, e.g. 20 persons and
#'   4 items). `NULL` with `step`.
#' @param step,prior,simulate A sampler written by the user: `step(s, data)` returns the new state,
#'   `prior()` a state drawn from the prior, `simulate(s)` the data list given the state.
#' @param method `"geweke"` or `"sbc"`.
#' @param n_persons Persons per simulated data set.
#' @param n_iter Geweke: iterations of the successive-conditional chain (the first 10% are dropped;
#'   default 2e5, about 15 s for three items). SBC: iterations per fit (the first half is burn-in;
#'   default 1e4).
#' @param n_rep SBC: number of simulated data sets.
#' @param n_draws SBC: thinned posterior draws per fit (ranks 0 to `n_draws`).
#' @param collapse,itemtype Sampler settings, as in [gibbs()] (default: those of the model).
#' @param priors Priors of the simulation and the sampler, from [rtirt_priors()].
#' @param sampler_priors Priors given to the sampler (default `priors`).
#' @param probs Geweke: probabilities of the prior quantiles that are compared.
#' @param n_prior Geweke: prior draws for the reference distribution.
#' @param seed Random seed.
#' @param steps For a [block_model()]: the steps (default as in [gibbs()]).
#' @param monitor For a [block_model()]: the parameters compared (default: all but the person parameters);
#'   with `step`: the names of the state elements compared (required).
#' @param x An `rtirt_check` object.
#' @param n Rows printed.
#' @param ... Unused.
#' @return A data frame of class `rtirt_check` with one row per comparison (Geweke: parameter,
#'   probability, share in the chain, effective sample size, z, p-value, Holm-adjusted p-value; SBC:
#'   parameter, chi-square statistic, p-value, Holm-adjusted p-value) and attributes `method`,
#'   `passed` (all Holm-adjusted p-values above 0.01), `draws` (Geweke chain) or `ranks` (SBC).
#' @references Cook, S. R., Gelman, A., & Rubin, D. B. (2006). Validation of software for Bayesian
#'   models using posterior quantiles. *Journal of Computational and Graphical Statistics, 15*,
#'   675-692.
#'
#'   Geweke, J. (2004). Getting it right: Joint distribution tests of posterior simulators.
#'   *Journal of the American Statistical Association, 99*, 799-804.
#'
#'   Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková, K., Gelman, A., &
#'   Vehtari, A. (2025). Simulation-based calibration checking for Bayesian computation: The choice
#'   of test quantities shapes sensitivity. *Bayesian Analysis, 20*(2).
#'
#'   Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian
#'   inference algorithms with simulation-based calibration. arXiv:1804.06788.
#' @examples
#' \donttest{
#' cond <- set_cond(n_subj = 100, n_item = 3)
#' m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
#' narrow <- rtirt_priors(lambda = c(mean = 0, sd = 1))           # any proper prior; narrow mixes fast
#' check_sampler(m, n_iter = 1e5, priors = narrow, seed = 1)        # Geweke
#' check_sampler(m, method = "sbc", n_rep = 100, priors = narrow, seed = 1)
#' # a wrong prior in the sampler is detected
#' wrong <- rtirt_priors(lambda = c(mean = 0, sd = 1), a = c(mean = 1, sd = 0.3))
#' check_sampler(m, n_iter = 2e4, priors = narrow, sampler_priors = wrong, seed = 1)
#' }
#' @export
check_sampler <- function(model = NULL, method = c("geweke", "sbc"), n_persons = 20, n_iter = NULL, n_rep = 200, n_draws = 99,
                          collapse = NULL, itemtype = NULL, priors = rtirt_priors(), sampler_priors = priors,
                          probs = c(0.1, 0.25, 0.5, 0.75, 0.9), n_prior = 1e5, seed = NULL, steps = NULL, monitor = NULL,
                          step = NULL, prior = NULL, simulate = NULL, ...) {
  method <- match.arg(method)
  if (!is.null(step)) {
    if (!is.function(prior) || !is.function(simulate) || !length(monitor))
      stop("a sampler written by the user is checked with step, prior, simulate and monitor", call. = FALSE)
    if (method != "geweke") stop("samplers written by the user are checked with method = \"geweke\"", call. = FALSE)
    return(check_step(step, prior, simulate, monitor, n_iter %||% 2e4, if (missing(n_prior)) 4000 else n_prior, probs, seed %||% 1))
  }
  if (is.null(model)) stop("give a model, or a sampler (step, prior, simulate, monitor)", call. = FALSE)
  if (inherits(model, "block_model")) {
    if (method != "geweke") stop("declared models are checked with method = \"geweke\"", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    return(check_model_steps(model, steps, n_iter %||% 2e4, probs, min(n_prior, 1e4), monitor))
  }
  if (inherits(model, "rtirt_qset")) model <- model$models[[1]]
  if (!inherits(model, "rtirt")) stop("model must be a model object, e.g. rtirt_cross(data)", call. = FALSE)
  for (p in list(priors, sampler_priors)) if (!inherits(p, "rtirt_priors")) stop("priors must come from rtirt_priors()", call. = FALSE)
  if (is.null(priors$lambda) || is.null(sampler_priors$lambda))
    stop("the check needs priors that do not depend on the data; rtirt_priors_julia() centres lambda on the log times", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  S <- ck_setup(model, n_persons, collapse %||% isTRUE(model$post$settings$collapse), itemtype, priors, sampler_priors)
  out <- if (method == "geweke") ck_geweke(S, n_iter %||% 2e5, probs, n_prior) else ck_sbc(S, n_iter %||% 1e4, n_rep, n_draws)
  attr(out, "setting") <- sprintf("%s%s, %d persons x %d items, collapse = %s", class(model)[1], if (S$qr) sprintf(" (q = %s)", S$q) else "",
                                  S$N, S$K, S$collapse)
  out
}

#' @rdname check_sampler
#' @export
print.rtirt_check <- function(x, n = 10, ...) {
  m <- attr(x, "method")
  if (is.null(m)) { print.data.frame(x, ...); return(invisible(x)) }       # a subset of the table
  cat(if (m == "geweke") "Geweke (2004) joint distribution test" else "Simulation-based calibration (Talts et al., 2018)", "\n")
  cat(" ", attr(x, "setting"), "\n ", attr(x, "detail"), "\n", sep = "")
  o <- x[order(x$p.holm, x$p.value), , drop = FALSE]
  cat(sprintf("  %d comparisons; smallest Holm-adjusted p = %.3g: %s\n\n", nrow(x), min(x$p.holm),
              if (attr(x, "passed")) "no evidence against the sampler" else "the sampler does NOT match the prior and likelihood"))
  print.data.frame(utils::head(o, n), row.names = FALSE, digits = 3)
  if (m == "geweke" && isTRUE(min(x$ess, na.rm = TRUE) < 400))
    cat(sprintf("\n  Note: the smallest effective sample size is %d; with fewer than about 400 the z-tests are unreliable,\n  so raise n_iter (or narrow the priors of the slow parameters, the check holds for any proper prior).\n", min(x$ess, na.rm = TRUE)))
  invisible(x)
}

# ---- setup: model template, design and the sampler call ----------------------------------------------
ck_setup <- function(model, N, collapse, itemtype, priors, sampler_priors) {
  info <- .model_info[[class(model)[1]]]; St <- model$settings; D <- model$data; eng <- info$engine
  one_pl <- identical(itemtype %||% St$itemtype, "1pl")
  X <- if (info$X) {
    if (nrow(D$X) < N) stop(sprintf("the model has %d persons; n_persons must not exceed it", nrow(D$X)), call. = FALSE)
    X0 <- D$X[seq_len(N), , drop = FALSE]
    if (isTRUE(St$intercept)) cbind(`(Intercept)` = 1, X0) else X0
  } else matrix(0, N, 0)
  qr <- isTRUE(model$qr)
  list(eng = eng, rt = info$rt, qr = qr, q = if (qr) model$q else 0.5, one_pl = one_pl, collapse = collapse,
       fixed = St$speed_var == "fixed" && eng %in% c("rtirt", "cross"), X = X, P = ncol(X), N = N, K = ncol(D$Y),
       items = D$items, pr = priors, pv = prior_vector(sampler_priors),
       names = par_names(eng, D$items, colnames(X) %||% character()))
}

ck_run <- function(S, Y, L, init, n_iter, n_burn) {
  aj <- S$collapse && !S$one_pl
  switch(S$eng,
    mlirt  = .gibbs_mlirt(Y, S$X, S$one_pl, S$collapse, aj, n_iter, n_burn, init, S$pv),
    rtirt  = .gibbs_rtirt(Y, L, S$X, S$fixed, S$one_pl, S$collapse, aj, n_iter, n_burn, init, S$pv),
    latent = .gibbs_latent(Y, L, S$X, S$qr, S$q, S$one_pl, S$collapse, aj, n_iter, n_burn, init, S$pv),
    cross  = .gibbs_cross(Y, L, S$qr, S$q, S$fixed, S$one_pl, S$collapse, aj, n_iter, n_burn, init, S$pv))
}

# draws of n variances (or ALD scales) from their prior: half-t on the square root, or inverse gamma
ck_rvar <- function(n, x, fam) if (fam == "half-t") (abs(stats::rt(n, x[["df"]])) * x[["scale"]])^2 else 1 / stats::rgamma(n, x[["shape"]], rate = x[["scale"]])

# M draws of the item and structural parameters from the prior, in the state layout of the samplers
ck_prior <- function(S, M) {
  pr <- S$pr; K <- S$K; P <- S$P
  tn <- function(n, m, s, lo, hi) matrix(draw_tnorm(rep(m, n), s, lo, hi), ncol = K)
  st <- list(a = if (S$one_pl) matrix(1, M, K) else tn(M * K, pr$a[1], pr$a[2], 0, Inf),
             b = tn(M * K, pr$b[["mean"]], pr$b[["sd"]], pr$b[["lower"]], pr$b[["upper"]]))
  if (S$rt) {
    st$lambda <- tn(M * K, pr$lambda[["mean"]], pr$lambda[["sd"]], pr$lambda[["lower"]], Inf)
    st$sigma2 <- matrix(ck_rvar(M * K, pr$sigma2t, pr$sigma2t_family), M)
  }
  nb <- function(k) matrix(stats::rnorm(M * k, 0, pr$beta[["sd"]]), M)
  vs <- function() ck_rvar(M, pr$var_speed, pr$var_speed_family)
  switch(S$eng,
    mlirt = { st$beta <- nb(P) },
    rtirt = { st$beta1 <- nb(P); st$beta2 <- nb(P)
              if (S$fixed) { st$c <- stats::runif(M, -1, 1); st$v <- 1 - st$c^2 } else { st$c <- stats::rnorm(M, 0, pr$cov[["sd"]]); st$v <- vs() } },
    latent = { st$beta <- nb(P + 1); st$s <- vs() },
    cross = { st$rho <- matrix(stats::rnorm(M * K, 0, pr$rho[["sd"]]), M); st$s <- if (S$fixed) rep(1, M) else vs() })
  st
}

# the parameters as stored by the samplers (one row per draw)
ck_stored <- function(S, st) {
  m <- with(st, switch(S$eng,
    mlirt = cbind(a, b, beta), cross = cbind(a, b, lambda, sigma2, rho, s), latent = cbind(a, b, lambda, sigma2, beta, s),
    rtirt = cbind(a, b, lambda, sigma2, beta1, beta2, c / sqrt(c^2 + v), c^2 + v)))
  colnames(m) <- S$names
  m
}

# one state of the sampler: parameters from the prior, persons given the parameters
ck_state <- function(S) ck_persons(S, lapply(ck_prior(S, 1), function(x) if (is.matrix(x)) x[1, ] else x))

# person values (ability, speed) drawn from their distribution given the parameters p
ck_persons <- function(S, p) {
  N <- S$N; X <- S$X; P <- S$P
  xb <- function(b) if (P) drop(X %*% b) else rep(0, N)
  k <- ald_k(S$q)
  switch(S$eng,
    mlirt = { p$theta <- xb(p$beta) + stats::rnorm(N) },
    cross = { p$theta <- stats::rnorm(N); p$zeta <- stats::rnorm(N, 0, sqrt(p$s)) },
    rtirt = { e1 <- stats::rnorm(N); p$theta <- xb(p$beta1) + e1; p$zeta <- xb(p$beta2) + p$c * e1 + stats::rnorm(N, 0, sqrt(p$v)) },
    latent = { p$theta <- stats::rnorm(N)
               u <- if (S$qr) { nu <- stats::rexp(N, 1 / p$s); k[1] * nu + sqrt(k[2] * p$s * nu) * stats::rnorm(N) } else stats::rnorm(N, 0, sqrt(p$s))
               p$zeta <- xb(p$beta[seq_len(P)]) + p$beta[P + 1] * p$theta + u })
  p
}
ald_k <- function(q) c((1 - 2 * q) / (q * (1 - q)), 2 / (q * (1 - q)))

# data given the parameters and the person values
ck_data <- function(S, st) {
  N <- S$N; K <- S$K
  eta <- outer(st$theta, st$a) - matrix(st$a * st$b, N, K, byrow = TRUE)
  Y <- matrix(stats::rbinom(N * K, 1, stats::plogis(eta)), N)
  L <- NULL
  if (S$rt) {
    mu <- matrix(st$lambda, N, K, byrow = TRUE) - st$zeta - if (S$eng == "cross") outer(st$theta, st$rho) else 0
    s2 <- matrix(st$sigma2, N, K, byrow = TRUE)
    e <- if (S$qr && S$eng == "cross") { k <- ald_k(S$q); nu <- matrix(stats::rexp(N * K, 1 / s2), N); k[1] * nu + sqrt(k[2] * s2 * nu) * stats::rnorm(N * K) }
         else sqrt(s2) * stats::rnorm(N * K)
    L <- mu + e
  }
  list(Y = Y, L = L %||% matrix(0, N, K))
}

# ---- Geweke's successive-conditional simulator ------------------------------------------------------
ck_geweke <- function(S, n_iter, probs, n_prior) {
  st <- ck_state(S); D <- ck_data(S, st)
  G <- matrix(NA_real_, n_iter, length(S$names), dimnames = list(NULL, S$names))
  for (it in seq_len(n_iter)) {
    r <- ck_run(S, D$Y, D$L, st, 1L, 0L)
    G[it, ] <- r$draws[1, ]
    st <- ck_persons(S, lapply(r$last, as.vector))     # persons and data anew given the parameters
    D <- ck_data(S, st)
  }
  G <- G[-seq_len(round(n_iter / 10)), , drop = FALSE]
  R <- ck_stored(S, ck_prior(S, n_prior))
  geweke_table(G, R, probs, n_prior, sprintf("%d iterations (first 10%% dropped); P(parameter <= prior quantile) in the chain vs %g prior draws",
                                             n_iter, n_prior))
}

# chain G against prior draws R: proportions below prior quantiles, ESS-based z-tests, Holm
geweke_table <- function(G, R, probs, n_prior, detail) {
  rows <- list()
  for (j in colnames(G)) {
    if (stats::sd(R[, j]) == 0) next                   # fixed by the model (1PL a, fixed speed variance)
    qs <- stats::quantile(R[, j], probs, names = FALSE)
    for (k in seq_along(probs)) {
      I <- as.numeric(G[, j] <= qs[k]); ph <- mean(I); pk <- mean(R[, j] <= qs[k])
      ess <- if (stats::var(I) > 0) max(1, unname(coda::effectiveSize(I))) else NA      # NA: the chain never crossed the quantile
      se <- sqrt(pk * (1 - pk) / (if (is.na(ess)) 1 else ess) + pk * (1 - pk) / n_prior)
      z <- (ph - pk) / se
      rows[[length(rows) + 1]] <- data.frame(parameter = j, prob = probs[k], prior = pk, sampler = ph, ess = round(ess), z = z,
                                             p.value = 2 * stats::pnorm(-abs(z)))
    }
  }
  out <- do.call(rbind, rows); out$p.holm <- stats::p.adjust(out$p.value, "holm")
  structure(out, class = c("rtirt_check", "data.frame"), method = "geweke", passed = all(out$p.holm > 0.01), draws = G, detail = detail)
}

# ---- simulation-based calibration --------------------------------------------------------------------
ck_sbc <- function(S, n_iter, n_rep, n_draws) {
  n_burn <- n_iter %/% 2
  keep <- round(seq(n_burn + 1, n_iter, length.out = n_draws))
  ranks <- matrix(NA_integer_, n_rep, length(S$names), dimnames = list(NULL, S$names))
  for (l in seq_len(n_rep)) {
    st <- ck_state(S); D <- ck_data(S, st)
    init <- initial_values(S$eng, list(Y = D$Y, log_t = D$L), S$X, S$fixed)
    r <- ck_run(S, D$Y, D$L, init, as.integer(n_iter), as.integer(n_burn))
    truth <- ck_stored(S, lapply(st, function(x) matrix(x, 1)))
    dr <- r$draws[keep, , drop = FALSE]
    ranks[l, ] <- colSums(dr < matrix(truth, n_draws, ncol(dr), byrow = TRUE))
  }
  B <- 10; rows <- list()
  for (j in colnames(ranks)) {
    if (length(unique(ranks[, j])) == 1) next
    cnt <- tabulate(pmin((ranks[, j] * B) %/% (n_draws + 1), B - 1) + 1, B)
    chi <- sum((cnt - n_rep / B)^2 / (n_rep / B))
    rows[[length(rows) + 1]] <- data.frame(parameter = j, chisq = chi, df = B - 1, p.value = stats::pchisq(chi, B - 1, lower.tail = FALSE))
  }
  out <- do.call(rbind, rows); out$p.holm <- stats::p.adjust(out$p.value, "holm")
  structure(out, class = c("rtirt_check", "data.frame"), method = "sbc", passed = all(out$p.holm > 0.01), ranks = ranks,
            detail = sprintf("%d simulated data sets; ranks among %d of %d post-burn-in draws; chi-square over %d bins", n_rep, n_draws, n_iter - n_burn, B))
}
