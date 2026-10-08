# EM for declared models. Gibbs and EM use the same full conditionals: a Gibbs step draws a
# parameter from its full conditional; an EM step (ECM; Meng & Rubin, 1993) replaces what is not
# observed by its expectation and moves the parameter to the mode. So the steps of gibbs()
# have EM counterparts: draw omega ~ PG <-> E[omega] = tanh(eta/2)/(2 eta) (a minorizer; Polson,
# Scott & Windle, 2013); draw the person parameters <-> integrate them by adaptive Gauss-Hermite
# quadrature (E-step); step_gibbs() <-> the mode of the expected full conditional (CM step);
# step_mala() <-> gradient ascent on Q (generalized EM). Persons j are rows, items i columns.

#' ECM for a declared model
#'
#' Fits a [block_model()] by EM with the same steps as [gibbs()]. Gibbs and EM share the
#' full conditionals of the model: a step that a sampler uses to **draw** a parameter is used by EM to
#' move the parameter to the **mode** of its full conditional averaged over the latent variables.
#'
#' * **E-step:** the person parameters (`latent`, by default every `[person]` parameter) are
#'   integrated out by adaptive Gauss–Hermite quadrature: a product grid of `nodes` points per
#'   latent variable, centred at each person's posterior mean and scaled by the posterior SD. The
#'   centres are fixed while EM runs (so the objective is one fixed function and EM never lowers
#'   it), then recentred, for up to `cycles` rounds. The Pólya–Gamma variables of a logistic part are
#'   replaced by their expectation \eqn{E[\omega] = \tanh(\eta/2)/(2\eta)}{E[omega] = tanh(eta/2)/(2 eta)}
#'   at each node (a minorizer: the step still never lowers the objective).
#' * **`step_gibbs()`** becomes a conditional maximization: the normal part of the parameter is summed
#'   over the model and over the nodes (weighted by the posterior of each person), and the parameter
#'   moves to its mode \eqn{m/p}{m/p} (clipped to the bounds); a variance moves to the mode of its
#'   inverse-gamma conditional (`prior = "ml"`: \eqn{\beta/\alpha}{beta/alpha}).
#' * **`step_mala()`** becomes a few gradient-ascent steps with step halving on the same weighted
#'   objective (generalized EM), on the log scale for bounded parameters (no Jacobian: a mode does not
#'   need one).
#' * **`step_shift()`** is skipped: the latent variables are integrated out.
#'
#' The objective is the log marginal likelihood (`prior = "ml"`; the priors of the latent variables
#' stay, they define the marginal likelihood) or the log marginal posterior (`prior = "map"`, the
#' default: also the priors of the other parameters). Standard errors come from its Hessian (Louis'
#' formula: expected complete-data Hessian plus the covariance of the complete-data score, from the
#' derivatives of the declared expressions); the
#' posterior means and SDs of the latent variables are in `fit$latent`. EM runs through [run_em()]
#' (SQUAREM, the checks of [run_em()]). Use [check_ecm()] to check the steps.
#'
#' @section Theory of each step:
#' * **ECM.** EM needs only the conditional expectations of the complete-data sufficient statistics
#'   (Dempster, Laird & Rubin, 1977, Sec. 2) and never lowers the objective when every M-step raises
#'   Q (their Theorem 1, generalized EM); the M-step may be a sequence of conditional maximizations
#'   (Meng & Rubin, 1993).
#' * **Adaptive quadrature.** Nodes centred and scaled at each person's posterior (Naylor & Smith,
#'   1982; Liu & Pierce, 1994; Bock & Gibbons, 2021, Sec. 6.4). Schilling & Bock (2005) recentre at
#'   every EM iteration; here the centres are kept fixed within a round and recentred between rounds,
#'   as Bock & Gibbons (2021, pp. 210–211) suggest to save time and as in the pseudo-adaptive rule of
#'   Rizopoulos (2012). With fixed centres the approximated objective is one fixed function, so the
#'   monotonicity of Dempster et al. (1977, Theorem 1) applies to it (our argument; no source states
#'   it).
#' * **`quadrature = "ba81"`.** Common nodes and expected counts per node ("artificial data"):
#'   Bock & Aitkin (1981); Bock & Gibbons (2021, Sec. 6.4, Eqs. 6.17–6.22).
#' * **Logistic parts.** \eqn{E[\omega] = \tanh(\eta/2)/(2\eta)}{E[omega] = tanh(eta/2)/(2 eta)} (Polson, Scott & Windle,
#'   2013, Sec. 2.2) makes the M-step an exact EM step for the Pólya–Gamma augmented model (Scott &
#'   Sun, 2013); it equals the Jaakkola–Jordan bound, a minorize–maximize step that never lowers the
#'   likelihood (Durante & Rigon, 2019, Appendix).
#' * **Exact integration (`analytic`).** A latent variable u that enters only normal parts linearly,
#'   with variances free of u, has a normal conditional \eqn{u \mid \cdot \sim N(h/p, 1/p)}{u | rest ~ N(h/p, 1/p)} (completing the square;
#'   Bishop, 2006, Sec. 2.3; for speed in the hierarchical RT model, van der Linden, 2007, Eq. 42),
#'   and the integral over u is the integrand at the mode times \eqn{\sqrt{2\pi/p}}{sqrt(2 pi / p)}, exact for a
#'   quadratic log integrand (Bishop, 2006, Eq. 4.135). The RT-IRT marginal likelihood papers we found
#'   integrate speed numerically together with ability (Glas & van der Linden, 2010; Molenaar &
#'   Bolsinova, 2017); integrating it exactly is our reduction, valid only under the condition above
#'   (checked from the declaration; not for heteroscedastic or skewed speed).
#' * **Exact variables in the M-step.** When every part of an exact variable is linear in all
#'   latent variables jointly, with variances free of them, these parts are quadratic in the latent
#'   variables, and the expectation of a quadratic depends only on the mean and covariance
#'   (\eqn{E[x'Ax] = \mathrm{tr}(A\Sigma) + \mu'A\mu}{E[x'Ax] = tr(A Sigma) + mu'A mu}, for any distribution). The M-step therefore uses
#'   each person's posterior mean and covariance of the latent variables through the 2L symmetric
#'   points \eqn{m \pm \sqrt{L}}{m -/+ sqrt(L)} x columns of the Cholesky factor, weight 1/(2L) each, which have exactly
#'   that mean and covariance (Julier & Uhlmann, 2004, Eq. 11). Logistic parts stay on the nodes.
#' * **Standard errors.** Louis (1982): observed information = expected complete-data information
#'   minus the covariance of the complete-data score. Under `prior = "map"` the log prior of the
#'   fixed parameters does not involve the latent variables, so its Hessian is added to that of the
#'   log marginal likelihood (Chalmers, 2018, notes that how priors enter Louis' form needs saying). The
#'   squared scores are quartic in an exact variable, so three Gauss–Hermite points per exact
#'   variable are used (n points are exact to degree 2n - 1; DLMF Sec. 3.5(v)).
#'
#' @param model A [block_model()].
#' @param steps Steps, as for [gibbs()] (the same list can be used for both); by default
#'   Gibbs (conditional maximization) where possible and MALA (gradient) otherwise.
#' @param latent Parameters integrated out in the E-step (one to three `[person]` parameters with a
#'   normal prior).
#' @param nodes Quadrature points per latent variable (default 21, 11 or 7 for one, two or three); one
#'   number, or one per latent variable (named, e.g. `c(theta = 11, zeta = 5)`).
#' @param prior `"map"` (posterior mode) or `"ml"` (maximum marginal likelihood).
#' @param init Optional starting values (named list).
#' @param cycles Maximum number of recentring rounds of the quadrature.
#' @param max_iter,tol Passed to [run_em()].
#' @param se Standard errors from the Hessian of the objective.
#' @param analytic Integrate latent variables that enter the model only through normal parts, linearly
#'   (e.g. speed in a log-normal RT model), exactly given the others instead of on the grid.
#' @param quadrature `"adaptive"` (default): nodes centred and scaled per person. `"ba81"`: the same
#'   nodes for every person (Bock & Aitkin, 1981; default 41, 21 or 11 per latent variable), centred
#'   at the marginal posterior. A part of the data that depends on the persons only through the
#'   variables on the grid is then summed over persons at each node ("artificial data"), so the E-step
#'   is a matrix product and the M-step works on nodes x items: faster for many persons. A latent
#'   variable whose posterior is narrow compared with its spread over persons needs many common nodes;
#'   the fit compares its objective with adaptive nodes at the estimates (`fit$adaptive_gap`) and
#'   warns when they differ by more than 0.5 (with `analytic = FALSE` speed is then on the common grid
#'   too, and its posterior is narrow when there are many items: use `analytic = TRUE`). Still ECM:
#'   only the numerical integration of the E-step changes, the conditional maximizations are the same.
#' @return A `block_em` (see [run_em()]) with `latent` (posterior means and SDs per person), `model`,
#'   `steps` and `cycles`; `algorithm(fit)` lists the E- and M-steps.
#' @references Bishop, C. M. (2006). *Pattern recognition and machine learning*. Springer.
#'
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm. *Psychometrika, 46*(4), 443–459. \doi{10.1007/BF02293801}
#'
#' Bock, R. D., & Gibbons, R. D. (2021). *Item response theory*. Wiley. \doi{10.1002/9781119716723}
#'
#' Chalmers, R. P. (2018). Numerical approximation of the observed information matrix with Oakes' identity. *British Journal of Mathematical and Statistical Psychology, 71*(3), 415–436. \doi{10.1111/bmsp.12127}
#'
#' Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society B, 39*(1), 1–22. \doi{10.1111/j.2517-6161.1977.tb01600.x}
#'
#' DLMF. *NIST Digital Library of Mathematical Functions*, Sec. 3.5(v), Gauss quadrature. https://dlmf.nist.gov/3.5
#'
#' Durante, D., & Rigon, T. (2019). Conditionally conjugate mean-field variational Bayes for logistic models. *Statistical Science, 34*(3), 472–485. \doi{10.1214/19-STS712}
#'
#' Glas, C. A. W., & van der Linden, W. J. (2010). Marginal likelihood inference for a model for item responses and response times. *British Journal of Mathematical and Statistical Psychology, 63*(3), 603–626. \doi{10.1348/000711009X481360}
#'
#' Julier, S. J., & Uhlmann, J. K. (2004). Unscented filtering and nonlinear estimation. *Proceedings of the IEEE, 92*(3), 401–422. \doi{10.1109/JPROC.2003.823141}
#'
#' Liu, Q., & Pierce, D. A. (1994). A note on Gauss–Hermite quadrature. *Biometrika, 81*(3), 624–629. \doi{10.1093/biomet/81.3.624}
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the EM algorithm. *Journal of the Royal Statistical Society B, 44*(2), 226–233. \doi{10.1111/j.2517-6161.1982.tb01203.x}
#'
#' Meng, X.-L., & Rubin, D. B. (1993). Maximum likelihood estimation via the ECM algorithm: A general framework. *Biometrika, 80*(2), 267–278. \doi{10.1093/biomet/80.2.267}
#'
#' Molenaar, D., & Bolsinova, M. (2017). A heteroscedastic generalized linear model with a non-normal speed factor for responses and response times. *British Journal of Mathematical and Statistical Psychology, 70*(2), 297–316. \doi{10.1111/bmsp.12087}
#'
#' Naylor, J. C., & Smith, A. F. M. (1982). Applications of a method for the efficient computation of posterior distributions. *Applied Statistics, 31*(3), 214–225. \doi{10.2307/2347995}
#'
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for logistic models using Pólya–Gamma latent variables. *Journal of the American Statistical Association, 108*(504), 1339–1349. \doi{10.1080/01621459.2013.829001}
#'
#' Rizopoulos, D. (2012). Fast fitting of joint models for longitudinal and event time data using a pseudo-adaptive Gaussian quadrature rule. *Computational Statistics & Data Analysis, 56*(3), 491–501. \doi{10.1016/j.csda.2011.09.007}
#'
#' Schilling, S., & Bock, R. D. (2005). High-dimensional maximum marginal likelihood item factor analysis by adaptive quadrature. *Psychometrika, 70*(3), 533–555. \doi{10.1007/s11336-003-1141-x}
#'
#' Scott, J. G., & Sun, L. (2013). Expectation-maximization for logistic regression. arXiv:1306.0040.
#'
#' van der Linden, W. J. (2007). A hierarchical framework for modeling speed and accuracy on test items. *Psychometrika, 72*(3), 287–308. \doi{10.1007/s11336-006-1478-z}
#' @examples
#' set.seed(1); N <- 300; K <- 6
#' th <- rnorm(N); a0 <- runif(K, 0.8, 2); b0 <- seq(-1.5, 1.5, length.out = K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(th, a0), 2, a0 * b0))), N)
#' m <- block_model(Y = Y,
#'   theta[person] ~ normal(0, 1),
#'   a[item] ~ normal(1, 1, lower = 0),
#'   b[item] ~ normal(0, 3),
#'   eta <- a * (theta - b),
#'   Y ~ bernoulli_logit(eta))
#' fit <- ecm(m, list(step_gibbs(a, b)), prior = "ml")
#' fit
#' algorithm(fit)
#' @export
ecm.block_model <- function(model, steps = NULL, latent = NULL, nodes = NULL, prior = c("map", "ml"), init = NULL,
                            cycles = 6, max_iter = 1000, tol = 1e-10, se = TRUE, analytic = TRUE, quadrature = c("adaptive", "ba81"), ...) {
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  bm_fit(model, steps, latent, nodes, match.arg(prior), init, cycles, max_iter, tol, se, analytic, match.arg(quadrature), "ecm")
}

#' Direct maximization for a declared model
#'
#' Maximizes the objective of [ecm.block_model()] (the marginal log-likelihood, or the log posterior
#' with `prior = "map"`) directly with [stats::nlminb()] (within the bounds of the parameters), all
#' parameters together, with the gradient from Fisher's identity and, for `method = "newton"`, the
#' Hessian from Louis' formula (Glas & van der Linden, 2010, Sec. 3.1, use Newton-Raphson or EM). No
#' closed-form step is needed, only the derivatives of the declared expressions, so this helps when a
#' parameter has no closed-form step. It reaches the same maximum as [ecm()] inside the bounds of the
#' parameters (ECM keeps a bounded parameter 1e-3 inside its bound; a parameter at its bound gets no
#' standard error, with a warning); ECM never lowers the objective and its iterations are cheap,
#' Newton needs few iterations.
#' @inheritParams ecm.block_model
#' @param method `"newton"` (Hessian by Louis' formula) or `"quasi-newton"`.
#' @return As [ecm.block_model()].
#' @references Glas, C. A. W., & van der Linden, W. J. (2010). Marginal likelihood inference for a
#'   model for item responses and response times. *British Journal of Mathematical and Statistical
#'   Psychology, 63*, 603-626.
#' @export
mml <- function(model, ...) UseMethod("mml")

#' @export
mml.default <- function(model, ...) stop("mml() maximizes declared models directly: model must come from block_model()", call. = FALSE)

#' @rdname mml
#' @export
mml.block_model <- function(model, latent = NULL, nodes = NULL, prior = c("map", "ml"), init = NULL, cycles = 6, max_iter = 1000,
                            tol = 1e-9, se = TRUE, analytic = TRUE, quadrature = c("adaptive", "ba81"), method = c("newton", "quasi-newton"), ...) {
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  bm_fit(model, NULL, latent, nodes, match.arg(prior), init, cycles, max_iter, tol, se, analytic, match.arg(quadrature), match.arg(method))
}

# the fit of a declared model: ECM steps (method = "ecm") or direct maximization ("newton", "quasi-newton")
bm_fit <- function(model, steps, latent, nodes, prior, init, cycles, max_iter, tol, se, analytic, quadrature, method) {
  E0 <- bm_setup(model, if (method == "ecm") steps, latent, nodes, prior, init, analytic, quadrature)
  s <- E0$s; C <- bm_centres(E0, s, NULL); ob <- -Inf; n_it <- 0
  for (cy in seq_len(cycles)) {
    fit <- if (method == "ecm") suppressWarnings(run_em(model, s, bm_step(E0, C), loglik = function(x, d) bm_estep(E0, x, C)$obj,
                                                         positive = E0$positive, se = FALSE, max_iter = max_iter, tol = tol))
           else bm_direct(E0, s, C, se = FALSE, max_iter = max_iter, tol = tol, newton = method == "newton")
    s <- fit$state[E0$fixed]; n_it <- n_it + fit$iterations
    e <- bm_estep(E0, s, C); C <- bm_centres(E0, s, e)
    done <- abs(e$obj - ob) < 1e-6 * max(1, abs(e$obj)); ob <- e$obj
    if (done) break
  }
  warns <- character()                                                 # held back: a maximum on a bound has a nonzero score
  fit <- withCallingHandlers(
    if (method == "ecm") run_em(model, s, bm_step(E0, C), loglik = function(x, d) bm_estep(E0, x, C)$obj, positive = E0$positive,
                                se = se, max_iter = max_iter, tol = tol, score = bm_score(E0, C), hessian = bm_information(E0, C))
    else bm_direct(E0, s, C, se = se, max_iter = max_iter, tol = tol, newton = method == "newton"),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  fit$iterations <- n_it + fit$iterations
  M <- model; bd <- unlist(lapply(E0$fixed, function(v) { P <- M$pars[[v]]; n <- length(fit$state[[v]])  # parameters at a bound of their prior
    if (P$dist != "normal") return(rep(FALSE, n)); x <- fit$state[[v]]
    (is.finite(P$lower) & x - P$lower <= 1.001e-3) | (is.finite(P$upper) & P$upper - x <= 1.001e-3) }))
  if (any(bd)) {
    fit$estimates$se[bd] <- NA
    warns <- c(grep("score statistic", warns, value = TRUE, invert = TRUE),
               sprintf("at the bound of its prior (no standard error; the maximum is on the boundary): %s", paste(fit$estimates$parameter[bd], collapse = ", ")))
  }
  for (w in warns) warning(w, call. = FALSE)
  if (E0$quadrature == "ba81") {                                       # common nodes can be too coarse: compare with adaptive nodes at the estimates
    Ea <- bm_setup(model, if (method == "ecm") steps, latent, NULL, E0$prior, init, analytic, "adaptive"); sa <- fit$state[E0$fixed]
    Ca <- bm_centres(Ea, sa, NULL); for (k in 1:3) Ca <- bm_centres(Ea, sa, bm_estep(Ea, sa, Ca))
    fit$adaptive_gap <- bm_estep(Ea, sa, Ca)$obj - fit$loglik
    if (abs(fit$adaptive_gap) > 0.5)                                    # e.g. a latent variable with a narrow posterior (speed with many items)
      warning(sprintf("quadrature = \"ba81\": at the estimates the objective with adaptive nodes differs by %.2f; the common nodes are too coarse for this posterior, so the %s and the estimates are not accurate. Use quadrature = \"adaptive\"%s",
                      fit$adaptive_gap, if (E0$prior == "map") "log posterior" else "log-likelihood",
                      if (!isTRUE(analytic)) ", or analytic = TRUE (latent variables that enter only normal parts linearly, such as speed, are then integrated exactly)" else " (or more nodes)"), call. = FALSE)
  }
  e <- bm_estep(E0, fit$state[E0$fixed], C)
  fit$latent <- lapply(stats::setNames(E0$latent, E0$latent), function(v) {
    if (E0$sigma && v %in% E0$ana) { X <- matrix(e$mu[[v]], model$N); X2 <- X^2 + matrix(e$vv[[v]], model$N) }
    else { X <- matrix(e$X[[v]], model$N); X2 <- X^2 }
    m <- rowSums(e$W * X); data.frame(mean = m, sd = sqrt(pmax(rowSums(e$W * X2) - m^2, 0)))
  })
  fit$objective <- if (E0$prior == "map") "log posterior" else "log-likelihood"
  fit$model <- model; fit$steps <- E0$steps; fit$cycles <- cy; fit$ecm <- list(latent = E0$latent, nodes = E0$Q, prior = E0$prior, C = C, analytic = analytic, quadrature = E0$quadrature, method = method)
  class(fit) <- c("block_model_ecm", class(fit))
  fit
}

# Direct maximization of the same objective (Glas & van der Linden, 2010, Sec. 3.1: Newton-Raphson or
# EM) by nlminb() (PORT, with the bounds of the parameters): the gradient from Fisher's identity
# (bm_score) and, for Newton steps, the Hessian from Louis' formula (bm_information), which also gives
# the standard errors.
bm_direct <- function(E, s, C, se = TRUE, max_iter = 1000, tol = 1e-9, newton = TRUE) {
  M <- E$M; v <- E$fixed; len <- vapply(v, function(n) length(s[[n]]), 0L); skel <- s[v]
  lo <- rep(vapply(v, function(n) { P <- M$pars[[n]]; if (P$dist == "normal") P$lower else 1e-8 }, 0), len)
  hi <- rep(vapply(v, function(n) { P <- M$pars[[n]]; if (P$dist == "normal") P$upper else Inf }, 0), len)
  st <- function(x) utils::relist(x, skel)
  sc <- bm_score(E, C); info <- bm_information(E, C)
  f <- function(x) { o <- bm_estep(E, st(x), C)$obj; if (is.finite(o)) -o else Inf }
  g <- function(x) -unlist(sc(st(x), NULL)[v])
  h <- if (newton) function(x) -info(st(x), NULL)
  x0 <- pmin(pmax(unlist(skel), lo + 1e-3 * is.finite(lo)), hi - 1e-3 * is.finite(hi))
  o <- stats::nlminb(x0, f, g, h, lower = lo, upper = hi,
                     control = list(eval.max = 2 * max_iter, iter.max = max_iter, rel.tol = tol))
  x <- o$par; state <- st(x)
  nm <- unlist(lapply(v, function(n) if (len[[n]] == 1) n else sprintf("%s[%d]", n, seq_len(len[[n]]))))
  V <- NULL; sev <- rep(NA_real_, length(x))
  if (se && o$convergence == 0) {
    H <- info(state, NULL)
    V <- tryCatch(solve(-H), error = function(e) NULL)
    if (!is.null(V) && all(eigen(-(H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values > 0)) sev <- sqrt(pmax(diag(V), 0))
    else { V <- NULL; warning("the Hessian of the objective is not negative definite (a parameter at its bound?): no standard errors", call. = FALSE) }
  }
  if (o$convergence != 0) warning(sprintf("nlminb() did not converge: %s", o$message), call. = FALSE)
  structure(list(estimates = data.frame(parameter = nm, est = x, se = sev, row.names = NULL), state = state, loglik = -o$objective,
                 iterations = o$iterations, converged = o$convergence == 0, trace = NULL, vcov = V, has_loglik = TRUE), class = "block_em")
}

# everything EM needs that does not change while it runs
bm_setup <- function(M, steps, latent, nodes, prior, init, analytic = TRUE, quadrature = "adaptive", sigma = TRUE, from_data = TRUE) {
  if (!inherits(M, "block_model")) stop("model must come from block_model()", call. = FALSE)
  pn <- names(M$pars)
  latent <- latent %||% pn[vapply(pn, function(v) M$pars[[v]]$dim == "person", TRUE)]
  bad <- setdiff(latent, pn); if (length(bad)) stop("latent: not parameters of the model: ", paste(bad, collapse = ", "), call. = FALSE)
  for (v in latent) {
    P <- M$pars[[v]]
    if (P$dim != "person" || P$dist != "normal" || is.finite(P$lower) || is.finite(P$upper))
      stop(sprintf("%s is integrated out, so it needs [person] and an untruncated normal() prior", v), call. = FALSE)
  }
  if (!length(latent)) stop("no latent [person] parameter to integrate out", call. = FALSE)
  if (length(latent) > 3) stop("at most three latent person parameters (the grid grows as nodes^L)", call. = FALSE)
  ord <- character(); todo <- latent                                    # latent variables in the order of their priors
  while (length(todo)) {
    ready <- todo[vapply(todo, function(v) !any(setdiff(intersect(c(all.vars(M$pars[[v]]$mean), all.vars(M$pars[[v]]$sd)), latent), ord) %in% todo), TRUE)]
    if (!length(ready)) stop("the priors of the latent variables refer to each other in a circle", call. = FALSE)
    ord <- c(ord, ready); todo <- setdiff(todo, ready)
  }
  fixed <- setdiff(pn, latent)
  steps <- bm_prepare_steps(M, steps, latent, fixed)
  ana <- character()                                                    # latent variables integrated exactly: normal and linear
  if (isTRUE(analytic)) for (u in rev(ord)) {
    fs <- factors_of(M, u)
    ok <- all(vapply(fs, function(f) f$dist == "normal" && !(u %in% f$sd_vars) && is_zero(tryCatch(stats::D(f$d1[[u]], u), error = function(e) 1)), TRUE)) &&
      !any(vapply(fs, function(f) length(intersect(setdiff(f$vars, u), ana)) > 0, TRUE))
    if (ok && length(setdiff(ord, c(ana, u)))) ana <- c(ana, u)          # keep at least one variable on the grid
  }
  grid <- setdiff(ord, ana)
  Q <- rep_len(nodes %||% (if (quadrature == "ba81") c(41, 21, 11) else c(21, 11, 7))[length(grid)], length(grid))
  if (!is.null(names(nodes))) Q <- unname(nodes[grid])
  gh <- lapply(Q, calc_nodes)
  Z <- as.matrix(expand.grid(lapply(gh, `[[`, "z"))); lw <- rowSums(as.matrix(expand.grid(lapply(gh, `[[`, "logw"))))
  S <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(ana))))       # two points per exact variable: m -/+ sd, exact for quadratic terms
  s <- init_model_state(M, list(), init, from_data)[fixed]
  positive <- fixed[vapply(fixed, function(v) { P <- M$pars[[v]]; P$dist != "normal" || identical(P$lower, 0) }, TRUE)]
  pf <- lapply(stats::setNames(pn, pn), function(v) M$factors[[which(vapply(M$factors, function(f) !f$observed && f$y == v, TRUE))]])
  ng <- nrow(Z); ns <- if (length(ana)) nrow(S) else 1
  fgrid <- Filter(function(f) (f$observed || f$y %in% latent) && !any(f$vars %in% ana), M$factors)
  fana <- lapply(stats::setNames(ana, ana), function(u) factors_of(M, u))
  # affine: every part of an exact variable is linear in all latent variables jointly, with a variance free of them. Then the exact
  # step reduces to sums over items per person (p, and h linear in the grid variables), and every expected log density the M-step
  # needs is quadratic in the latent variables, so it depends only on each person's posterior mean and covariance (2L points below).
  affine <- length(ana) > 0 && all(vapply(unlist(fana, recursive = FALSE), function(f) {
    lv <- intersect(f$vars, ord)
    !length(intersect(f$sd_vars, ord)) && all(vapply(lv, function(a) all(vapply(lv, function(b) is_zero(tryCatch(stats::D(f$d1[[a]], b), error = function(e) 1)), TRUE)), TRUE))
  }, TRUE))
  sigma <- affine && sigma                                             # sigma = FALSE: points per exact variable at every node (for tests)
  if (sigma) ns <- 1                                                    # no points per exact variable: the M-step uses the 2L points
  E <- list(M = M, Mg = bm_stackM(M, ng), Ms = bm_stackM(M, ng * ns), latent = ord, grid = grid, ana = ana, fixed = fixed, steps = steps,
            Q = Q, Z = Z, lw = lw, S = S, ng = ng, ns = ns, s = s, positive = positive, prior = prior, pf = pf, fgrid = fgrid, fana = fana,
            quadrature = quadrature, collapse = integer(), affine = affine, sigma = sigma, flat = character(),
            Msig = if (sigma) bm_stackM(M, 2L * length(ord)))
  if (quadrature == "ba81") E <- bm_ba81_setup(E)
  if (prior == "ml") for (v in fixed) if (!any(vapply(factors_of(M, v), function(f) bm_in_q(E, f), TRUE)))
    stop(sprintf("ecm(prior = \"ml\"): %s appears only in priors of other parameters, so the likelihood does not depend on it; use prior = \"map\"", v), call. = FALSE)
  E
}

# Common nodes (Bock & Aitkin, 1981): the same grid for every person. A part of the data that depends
# on the persons only through the grid variables then has the same value at a node for everyone, so
# the E-step is a matrix product and the M-step reads "artificial data" per node and item: for a
# logistic part the expected numbers of 1s and 0s, for a normal part the expected count, mean and
# spread (two pseudo-observations mean -/+ SD, weight count/2, give the same sums of y and y^2).
bm_ba81_setup <- function(E) {
  M <- E$M; pn <- names(M$pars)
  ok <- function(f) {
    if (!f$observed || f$shape != "NK" || !f$dist %in% c("normal", "bernoulli_logit") || !any(f$vars %in% E$grid)) return(FALSE)
    psy <- Filter(function(n) if (n %in% pn) M$pars[[n]]$dim == "person" else isTRUE(M$ddim[[n]] %in% c("person", "NK")), setdiff(f$syms, f$y))
    all(psy %in% E$grid)
  }
  ids <- vapply(Filter(ok, E$fgrid), `[[`, 0L, "id")
  Mc <- M; Mc$N <- 2L * E$ng; cs <- list()
  for (id in ids) { f <- M$factors[[id]]; mk <- fmask(f, M) * matrix(1, M$N, M$K); y0 <- M$data0[[f$y]]
    cs[[as.character(id)]] <- list(M = mk, Y = mk * y0, Y2 = mk * y0^2)
    Mc$data[[f$y]] <- Mc$data0[[f$y]] <- matrix(0, 2L * E$ng, M$K) }
  E$collapse <- ids; E$Mc0 <- Mc; E$cs <- cs
  E
}
# the grid variables at the nodes (twice: two pseudo-observations per node) and the other parameters
bm_nodes <- function(E, s, C) {
  out <- s[E$fixed]
  for (l in seq_along(E$grid)) { v <- E$grid[l]; out[[v]] <- rep(C[[v]]$mu[1] + C[[v]]$sd[1] * E$Z[, l], 2) }
  out
}
# log density of the collapsed parts at every person and node (N x nodes), by matrix products
bm_collapsed_ll <- function(E, sn) {
  LL <- 0; top <- seq_len(E$ng); Mc <- E$Mc0
  for (id in E$collapse) { fx <- Mc$factors[[id]]; cs <- E$cs[[as.character(id)]]
    at <- function(key) full_of(fval(fx, sn, Mc, key), fx, Mc)[top, , drop = FALSE]
    if (fx$dist == "bernoulli_logit") { eta <- at("eta"); LL <- LL + cs$Y %*% t(eta) - cs$M %*% t(softplus(eta)) }
    else { mu <- at("mean"); s2 <- at("sd2")
      LL <- LL - 0.5 * (cs$Y2 %*% t(1 / s2) - 2 * cs$Y %*% t(mu / s2) + cs$M %*% t(mu^2 / s2 + log(s2) + log(2 * pi))) }
  }
  LL
}
# the artificial data of the collapsed parts under the node weights Wg (N x nodes)
bm_collapse <- function(E, Wg) {
  Mc <- E$Mc0; tW <- t(Wg); ng <- E$ng; K <- E$M$K
  for (id in E$collapse) { f <- Mc$factors[[id]]; cs <- E$cs[[as.character(id)]]
    n <- tW %*% cs$M; s1 <- tW %*% cs$Y
    if (f$dist == "bernoulli_logit") {
      y <- rbind(matrix(1, ng, K), matrix(0, ng, K)); w <- rbind(s1, pmax(n - s1, 0)); Mc$factors[[id]]$kappa <- w * (y - 0.5)
    } else {
      yb <- ifelse(n > 0, s1 / n, 0); dd <- sqrt(pmax(ifelse(n > 0, (tW %*% cs$Y2) / n, 0) - yb^2, 0))
      y <- rbind(yb + dd, yb - dd); w <- rbind(n, n) / 2
    }
    Mc$data[[f$y]] <- Mc$data0[[f$y]] <- y; Mc$factors[[id]]$mask <- w
  }
  Mc
}

# is factor f part of the objective? With prior = "ml" the priors of the fixed parameters are not (the
# priors of the latent variables are: they define the marginal likelihood); nor are those of the parameters
# added by score_test(add = ) (E$flat)
bm_in_q <- function(E, f) !(!f$observed && f$y %in% E$fixed && (E$prior == "ml" || f$y %in% E$flat))

# the model with every person repeated reps times (node-major), so that all nodes are evaluated at once
bm_stackM <- function(M, reps) {
  rows <- rep(seq_len(M$N), reps); Ms <- M; Ms$N <- M$N * reps
  for (n in names(M$data)) for (slot in c("data", "data0"))
    Ms[[slot]][[n]] <- switch(M$ddim[[n]], NK = M[[slot]][[n]][rows, , drop = FALSE], person = M[[slot]][[n]][rows], M[[slot]][[n]])
  for (k in seq_along(M$factors)) { f <- M$factors[[k]]
    sub <- function(x) if (is.matrix(x)) x[rows, , drop = FALSE] else if (f$shape == "person") x[rows] else x
    if (!is.null(f$mask)) Ms$factors[[k]]$mask <- sub(f$mask)
    if (!is.null(f$kappa)) Ms$factors[[k]]$kappa <- sub(f$kappa) }
  Ms
}

bm_prepare_steps <- function(M, steps, latent, fixed) {
  if (is.null(steps)) {
    kinds <- vapply(fixed, function(v) gibbs_kind(M, v), "")
    g <- fixed[kinds %in% c("normal", "variance")]; o <- setdiff(fixed, g)
    steps <- c(if (length(g)) list(do.call(step_gibbs, lapply(g, as.name))),
               lapply(split(o, vapply(o, function(v) M$pars[[v]]$dim, "")), function(vs) do.call(step_mala, lapply(vs, as.name))))
  }
  if (inherits(steps, "block_step")) steps <- list(steps)
  out <- list()
  for (st in steps) {
    if (length(st$cut)) stop("a cut (step_gibbs(..., cut = )) is a sampling scheme; use gibbs()", call. = FALSE)
    if (st$type == "shift") next                                   # the latent variables are integrated out
    st$vars <- setdiff(st$vars, latent)
    if (!length(st$vars)) next
    bad <- setdiff(st$vars, fixed); if (length(bad)) stop("not a parameter of the model: ", paste(bad, collapse = ", "), call. = FALSE)
    if (st$type == "gibbs") for (v in st$vars) {
      kd <- gibbs_kind(M, v)
      if (!kd %in% c("normal", "variance")) stop(sprintf("no conditional mode in closed form for %s: %s; use step_mala(%s)", v, kd, v), call. = FALSE)
    }
    if (st$type == "mala") {
      d <- unique(vapply(st$vars, function(v) M$pars[[v]]$dim, ""))
      if (length(d) > 1) stop("step_mala() moves parameters of one dimension together", call. = FALSE)
      for (v in st$vars) if (M$pars[[v]]$dist != "normal" || is.finite(M$pars[[v]]$upper))
        stop(sprintf("step_mala(%s) needs a normal prior without an upper bound", v), call. = FALSE)
      st$dim <- d
    }
    out[[length(out) + 1]] <- st
  }
  miss <- setdiff(fixed, unlist(lapply(out, `[[`, "vars")))
  if (length(miss)) stop(sprintf("no step updates %s; add step_gibbs() or step_mala()", paste(miss, collapse = ", ")), call. = FALSE)
  out
}

# Q(s | e): the expected complete-data log density (data and latent priors, each part where the M-step evaluates it)
bm_qvalue <- function(E, s, e) {
  q <- bm_logprior(E, s)
  for (f in E$M$factors) if (f$observed || f$y %in% E$latent) {
    g <- bm_where(E, f, s, e); r <- bm_rows(g$f, g$s, g$M); q <- q + sum(if (is.null(g$w)) r else g$w * r)
  }
  q
}
# log density of one factor, one value per (stacked) person
bm_rows <- function(f, s, M) {
  fx <- M$factors[[f$id]]; m <- fmask(fx, M)
  v <- if (fx$dist == "bernoulli_logit") { eta <- fval(fx, s, M, "eta"); m * (fval(fx, s, M, "y") * eta - softplus(eta)) }
       else { r <- fval(fx, s, M, "r"); s2 <- fval(fx, s, M, "sd2"); m * (-r^2 / (2 * s2) - log(s2) / 2 - 0.5 * log(2 * pi)) }
  if (fx$shape == "NK") rowSums(v) else v
}

# the state on the stacked nodes: latent values at every node, the other parameters repeated for person parameters
bm_stack <- function(E, s, X, nq = E$ng * E$ns) {
  out <- s[E$fixed]
  for (v in E$fixed) if (E$M$pars[[v]]$dim == "person") out[[v]] <- rep(s[[v]], nq)
  c(out, X)
}
# a factor evaluated on the stacked nodes (when it uses a latent variable) or once
bm_where <- function(E, f, s, e, collapse = TRUE) {
  if (collapse && f$id %in% E$collapse) return(list(M = e$Mc, f = e$Mc$factors[[f$id]], s = c(s[E$fixed], e$Xc), w = NULL, lat = FALSE))
  if (E$sigma && any(f$vars %in% E$ana)) {                              # a part of an exact variable: at the 2L points of each person
    out <- s[E$fixed]; for (v in E$fixed) if (E$M$pars[[v]]$dim == "person") out[[v]] <- rep(s[[v]], e$sig$n)
    return(list(M = E$Msig, f = E$Msig$factors[[f$id]], s = c(out, e$sig$X), w = e$sig$w, lat = TRUE))
  }
  if (any(f$vars %in% E$latent)) list(M = E$Ms, f = E$Ms$factors[[f$id]], s = bm_stack(E, s, e$X), w = as.vector(e$W), lat = TRUE)
  else list(M = E$M, f = f, s = s, w = NULL, lat = FALSE)
}
bm_fold <- function(E, x, dim, lat) if (lat && dim == "person") rowSums(matrix(x, E$M$N)) else x
bm_logprior <- function(E, s) {
  if (E$prior == "ml") return(0)
  sum(vapply(setdiff(E$fixed, E$flat), function(v) bm_logprior_one(E, v, s), 0))
}
bm_logprior_one <- function(E, v, s) {
  P <- E$M$pars[[v]]; f <- E$pf[[v]]; x <- s[[v]]
  sum(switch(P$dist,
    normal = stats::dnorm(x, fval(f, s, E$M, "mean"), fval(f, s, E$M, "sd"), log = TRUE),
    half_t = log(2) + stats::dt(sqrt(x) / P$hyper[2], P$hyper[1], log = TRUE) - log(P$hyper[2]) - log(2 * sqrt(x)),
    inv_gamma = P$hyper[1] * log(P$hyper[2]) - lgamma(P$hyper[1]) - (P$hyper[1] + 1) * log(x) - P$hyper[2] / x))
}

# E-step on the frozen grid C: node values of the latent variables, weights, objective
bm_estep <- function(E, s, C) {
  N <- E$M$N; ng <- E$ng; ns <- E$ns; Xg <- list(); corr <- 0
  for (l in seq_along(E$grid)) {
    v <- E$grid[l]; z <- rep(E$Z[, l], each = N); sdv <- rep(C[[v]]$sd, ng)
    Xg[[v]] <- rep(C[[v]]$mu, ng) + sdv * z; corr <- corr + log(sdv) + z^2 / 2 + 0.5 * log(2 * pi)
  }
  sg <- bm_stack(E, s, Xg, ng); ll <- corr
  for (f in E$fgrid) if (!(f$id %in% E$collapse)) ll <- ll + bm_rows(f, sg, E$Mg)
  if (length(E$collapse)) { sn <- bm_nodes(E, s, C); ll <- ll + as.vector(bm_collapsed_ll(E, sn)) }
  mu <- list(); vv <- list()
  for (u in E$ana) {                                                  # exact: u | grid ~ N(h/p, 1/p); integral = f(m) sqrt(2 pi / p)
    if (E$affine) { a <- bm_exact_affine(E, u, s, Xg); mu[[u]] <- a$mu; vv[[u]] <- a$vv; ll <- ll + a$ll; next }
    su <- sg; su[[u]] <- rep(0, N * ng); p <- 0; h <- 0
    for (f in E$fana[[u]]) {
      fx <- E$Mg$factors[[f$id]]; m <- fmask(fx, E$Mg)
      c1 <- fval(fx, su, E$Mg, "d1", u); w <- m / fval(fx, su, E$Mg, "sd2"); r0 <- fval(fx, su, E$Mg, "r")
      p <- p + reduce_to(w * c1^2, "person", fx); h <- h - reduce_to(w * c1 * r0, "person", fx)
    }
    mu[[u]] <- h / p; vv[[u]] <- 1 / p; su[[u]] <- mu[[u]]
    for (f in E$fana[[u]]) ll <- ll + bm_rows(f, su, E$Mg)
    ll <- ll + 0.5 * log(2 * pi / p)
  }
  LL <- matrix(ll, N, ng) + matrix(E$lw, N, ng, byrow = TRUE)
  mx <- apply(LL, 1, max); lj <- mx + log(rowSums(exp(LL - mx))); Wg <- exp(LL - lj)
  if (E$sigma) {                                                       # grid variables at the nodes; exact ones through their moments
    out <- list(W = Wg, ll = lj, X = Xg, obj = sum(lj) + bm_logprior(E, s), mu = mu, vv = vv)
    out$sig <- bm_sigma(E, Wg, Xg, mu, vv)
  } else {
    X <- lapply(Xg, function(x) rep(x, ns))
    for (k in seq_along(E$ana)) { u <- E$ana[k]; X[[u]] <- as.vector(outer(mu[[u]], rep(1, ns)) + outer(sqrt(vv[[u]]), E$S[, k])) }
    sw <- E$Sw %||% rep(1 / ns, ns)                                     # weights of the points of the exact variables
    out <- list(W = matrix(rep(as.vector(Wg), ns) * rep(sw, each = N * ng), N), ll = lj, X = X, obj = sum(lj) + bm_logprior(E, s), mu = mu, vv = vv)
  }
  if (length(E$collapse)) { out$Mc <- bm_collapse(E, Wg); out$Xc <- sn[E$grid] }
  out
}

# The exact step for an affine part (completing the square). With w = mask / variance, c = dr/du and
# r0 = e0 + sum_l g_l x_l (r at u = 0; x_l the grid variables): p = sum_i w c^2 (one per person),
# h = -sum_i w c r0 = H0 + sum_l H_l x_l, and sum_i w r0^2 = R00 + 2 sum_l R0l x_l + sum_lk Rlk x_l x_k.
# At the conditional mode u = h/p, sum_i w r^2 = sum_i w r0^2 - h^2/p, so the integral over u is
# -1/2 (sum w r0^2 - h^2/p) - 1/2 sum log(2 pi variance) + 1/2 log(2 pi/p): sums over items once per person, not per node.
bm_exact_affine <- function(E, u, s, Xg) {
  M <- E$M; N <- M$N; ng <- E$ng
  s0 <- s[E$fixed]; for (v in E$latent) s0[[v]] <- numeric(N)          # every latent variable at 0
  gl <- intersect(E$grid, unique(unlist(lapply(E$fana[[u]], `[[`, "vars"))))
  p <- 0; H0 <- 0; R00 <- 0; cst <- 0; Hl <- list(); R0l <- list(); Rlk <- list()
  for (l in gl) { Hl[[l]] <- 0; R0l[[l]] <- 0; for (k in gl) Rlk[[paste(l, k)]] <- 0 }
  for (f in E$fana[[u]]) {
    red <- function(x) reduce_to(full_of(x, f, M), "person", f)
    m <- fmask(f, M) * full_of(1, f, M); s2 <- full_of(fval(f, s0, M, "sd2"), f, M); w <- m / s2
    c1 <- full_of(fval(f, s0, M, "d1", u), f, M); e0 <- full_of(fval(f, s0, M, "r"), f, M)
    g <- lapply(stats::setNames(gl, gl), function(l) if (l %in% f$vars) full_of(fval(f, s0, M, "d1", l), f, M) else 0)
    p <- p + red(w * c1^2); H0 <- H0 - red(w * c1 * e0); R00 <- R00 + red(w * e0^2)
    cst <- cst - 0.5 * red(m * log(2 * pi * s2))
    for (l in gl) { Hl[[l]] <- Hl[[l]] - red(w * c1 * g[[l]]); R0l[[l]] <- R0l[[l]] + red(w * e0 * g[[l]])
      for (k in gl) Rlk[[paste(l, k)]] <- Rlk[[paste(l, k)]] + red(w * g[[l]] * g[[k]]) }
  }
  rp <- function(x) rep_len(x, N * ng)
  h <- rp(H0); srr <- rp(R00)
  for (l in gl) { h <- h + rp(Hl[[l]]) * Xg[[l]]; srr <- srr + 2 * rp(R0l[[l]]) * Xg[[l]]
    for (k in gl) srr <- srr + rp(Rlk[[paste(l, k)]]) * Xg[[l]] * Xg[[k]] }
  P <- rp(p)
  list(mu = h / P, vv = 1 / P, ll = -0.5 * srr + h^2 / (2 * P) + rp(cst) + 0.5 * log(2 * pi / P))
}

# Each person's posterior mean m and covariance V of all latent variables (the grid ones from the node
# weights, the exact ones from their conditional normals), and 2L points m -/+ sqrt(L) * column of
# chol(V), weight 1/(2L) each: they have the same mean and covariance, so the expectation of any
# quadratic function of the latent variables is exact (unscented / third-degree cubature points).
bm_sigma <- function(E, Wg, Xg, mu, vv) {
  N <- E$M$N; lat <- E$latent; L <- length(lat)
  Z <- lapply(stats::setNames(lat, lat), function(v) matrix(if (v %in% E$ana) mu[[v]] else Xg[[v]], N))
  m <- lapply(Z, function(z) rowSums(Wg * z))
  C <- matrix(list(), L, L)
  for (a in seq_len(L)) for (b in seq_len(a)) {
    x <- rowSums(Wg * Z[[a]] * Z[[b]]) - m[[a]] * m[[b]]
    if (a == b && lat[a] %in% E$ana) x <- x + rowSums(Wg * matrix(vv[[lat[a]]], N))
    C[[a, b]] <- if (a == b) pmax(x, 0) else x
  }
  Lc <- matrix(list(0), L, L)                                          # Cholesky factor, one per person (vectorized over persons)
  for (k in seq_len(L)) {
    d <- C[[k, k]]; for (t in seq_len(k - 1)) d <- d - Lc[[k, t]]^2
    Lc[[k, k]] <- sqrt(pmax(d, 0))
    for (r in seq_len(L)[seq_len(L) > k]) { x <- C[[r, k]]; for (t in seq_len(k - 1)) x <- x - Lc[[r, t]] * Lc[[k, t]]
      Lc[[r, k]] <- ifelse(Lc[[k, k]] > 0, x / Lc[[k, k]], 0) }
  }
  X <- lapply(stats::setNames(seq_len(L), lat), function(r) unlist(lapply(seq_len(L), function(k) {
    d <- sqrt(L) * (if (r >= k) Lc[[r, k]] else 0); c(m[[r]] + d, m[[r]] - d) })))
  list(X = X, w = rep(1 / (2 * L), 2 * L * N), n = 2L * L)
}

# centres and scales of the grid: posterior means and SDs (the prior at the start)
bm_centres <- function(E, s, e) {
  N <- E$M$N; C <- list()
  if (is.null(e)) {
    sq <- s
    for (v in E$latent) { f <- E$pf[[v]]
      mu <- rep_len(fval(f, sq, E$M, "mean"), N); if (v %in% E$grid) C[[v]] <- list(mu = mu, sd = rep_len(fval(f, sq, E$M, "sd"), N)); sq[[v]] <- mu }
  } else for (v in E$grid) {
    Xv <- matrix(e$X[[v]], N); m <- rowSums(e$W * Xv)
    C[[v]] <- list(mu = m, sd = pmax(sqrt(pmax(rowSums(e$W * Xv^2) - m^2, 0)), 1e-3))
  }
  if (E$quadrature == "ba81") C <- lapply(C, function(c) { m <- mean(c$mu)     # one grid for everyone: the marginal mean and SD
    list(mu = rep(m, N), sd = rep(sqrt(mean(c$sd^2) + mean((c$mu - m)^2)), N)) })
  C
}

# expected normal part of a fixed parameter v over the model and the nodes
bm_normal_part <- function(E, v, s, e) {
  P <- E$M$pars[[v]]; prec <- 0; num <- 0
  for (f in factors_of(E$M, v)) {
    if (!bm_in_q(E, f)) next
    g <- bm_where(E, f, s, e); M <- g$M; fx <- g$f
    if (fx$dist != "normal" && !is.null(e$omc)) {                   # logistic part: one E[omega] per M-step for all its parameters
      om <- e$omc[[as.character(f$id)]]
      fast <- if (is.null(om)) ml_part_f(fx, v, g$s, M, "em", roww = g$w, keep_om = TRUE) else ml_part_f(fx, v, g$s, M, "gibbs", om = om, roww = g$w)
      if (!is.null(fast) && is.null(om) && length(fast) == 3) assign(as.character(f$id), fast[[3]], envir = e$omc)
    } else fast <- ml_part_f(fx, v, g$s, M, if (fx$dist == "normal") "normal" else "em", roww = g$w)
    if (!is.null(fast)) { prec <- prec + bm_fold(E, fast[[1]], P$dim, g$lat); num <- num + bm_fold(E, fast[[2]], P$dim, g$lat); next }
    m <- fmask(fx, M)
    s0 <- g$s; s0[[v]] <- 0 * s0[[v]]; c1 <- fval(fx, g$s, M, "d1", v)
    if (fx$dist == "normal") { wt <- m / fval(fx, g$s, M, "sd2"); r0 <- fval(fx, s0, M, "r"); Pp <- wt * c1^2; H <- -wt * c1 * r0 }
    else {
      om <- if (!is.null(e$omc)) e$omc[[as.character(f$id)]]               # the same E[omega] as the compiled path
      if (is.null(om)) { om <- m * pg_mean(fval(fx, g$s, M, "eta")); if (!is.null(e$omc)) assign(as.character(f$id), om, envir = e$omc) }
      Pp <- om * c1^2; H <- c1 * (fx$kappa - om * fval(fx, s0, M, "eta"))
    }
    if (g$lat) { Pp <- Pp * g$w; H <- H * g$w }
    prec <- prec + bm_fold(E, reduce_to(Pp, P$dim, fx), P$dim, g$lat); num <- num + bm_fold(E, reduce_to(H, P$dim, fx), P$dim, g$lat)
  }
  if (any(!(prec > 0))) stop(sprintf("the expected precision of %s is not positive (no information; give a prior and use prior = \"map\")", v), call. = FALSE)
  pmin(pmax(num / prec, P$lower + 1e-3 * is.finite(P$lower)), P$upper)   # keep a bounded parameter just inside, as update_items_logistic()
}

bm_variance <- function(E, v, s, e) {
  P <- E$M$pars[[v]]; n <- dim_len(E$M, P$dim); al <- numeric(n); be <- numeric(n)
  for (f in factors_of(E$M, v)) {
    if (!bm_in_q(E, f) || (!f$observed && f$y == v)) next
    g <- bm_where(E, f, s, e); M <- g$M; fx <- g$f
    m <- fmask(fx, M); r <- fval(fx, g$s, M, "r"); ww <- if (g$lat) m * g$w else m
    al <- al + bm_fold(E, reduce_to(ww * full_of(1, fx, M), P$dim, fx), P$dim, g$lat) / 2
    be <- be + bm_fold(E, reduce_to(ww * r^2, P$dim, fx), P$dim, g$lat) / 2
  }
  if (E$prior == "ml") return(be / al)
  if (P$dist == "inv_gamma") return((P$hyper[2] + be) / (P$hyper[1] + al + 1))
  vapply(seq_len(n), function(k) {                                    # half-t prior of the SD: one-dimensional maximization
    g <- function(lv) { x <- exp(lv); -al[k] * lv - be[k] / x + stats::dt(sqrt(x) / P$hyper[2], P$hyper[1], log = TRUE) - lv / 2 }
    exp(stats::optimize(g, c(-20, 20), maximum = TRUE)$maximum)
  }, 0)
}

# expected complete-data log posterior of the parameters vs (one dimension), per element, with gradient
bm_qtarget <- function(E, vs, s, e) {
  M <- E$M; d <- M$pars[[vs[1]]]$dim; lo <- vapply(M$pars[vs], `[[`, 0, "lower")
  facs <- Filter(function(f) length(intersect(f$vars, vs)) > 0 && bm_in_q(E, f), M$factors)
  function(x) {
    x <- matrix(x, ncol = length(vs)); s2 <- s
    for (k in seq_along(vs)) s2[[vs[k]]] <- if (is.finite(lo[k])) lo[k] + exp(x[, k]) else x[, k]
    lp <- 0; G <- matrix(0, nrow(x), length(vs))
    one <- function(f, sx, w, M, lat) {
      red <- function(x) bm_fold(E, reduce_to(x, d, f), d, lat)
      m <- fmask(f, M); if (!is.null(w)) m <- m * w
      if (f$dist == "normal") {
        r <- fval(f, sx, M, "r"); v2 <- fval(f, sx, M, "sd2")
        lp <<- lp + red(m * (-r^2 / (2 * v2) - log(v2) / 2))
        for (k in seq_along(vs)) if (vs[k] %in% f$vars) G[, k] <<- G[, k] + red(m * gnorm(f, sx, M, vs[k], r, v2))
      } else {
        eta <- fval(f, sx, M, "eta"); y <- fval(f, sx, M, "y")
        lp <<- lp + red(m * (y * eta - softplus(eta)))
        for (k in seq_along(vs)) if (vs[k] %in% f$vars)
          G[, k] <<- G[, k] + red(m * (y - stats::plogis(eta)) * fval(f, sx, M, "d1", vs[k]))
      }
    }
    for (f in facs) { g <- bm_where(E, f, s2, e); one(g$f, g$s, g$w, g$M, g$lat) }
    for (k in seq_along(vs)) if (is.finite(lo[k])) G[, k] <- G[, k] * exp(x[, k])     # chain rule only: a mode needs no Jacobian
    list(lp = lp, grad = G)
  }
}

bm_gradient <- function(E, vs, s, e, n_steps = 5) {
  lo <- vapply(E$M$pars[vs], `[[`, 0, "lower")
  x <- vapply(seq_along(vs), function(k) { z <- s[[vs[k]]]; if (is.finite(lo[k])) log(z - lo[k]) else z }, numeric(dim_len(E$M, E$M$pars[[vs[1]]]$dim)))
  x <- matrix(x, ncol = length(vs)); tg <- bm_qtarget(E, vs, s, e); cur <- tg(x); h <- rep(0.1, nrow(x))
  for (it in seq_len(n_steps)) {
    todo <- rep(TRUE, nrow(x))
    for (tr in 1:30) {
      xn <- x; xn[todo, ] <- x[todo, ] + h[todo] * cur$grad[todo, , drop = FALSE]
      nw <- tg(xn); ok <- todo & is.finite(nw$lp) & nw$lp >= cur$lp
      x[ok, ] <- xn[ok, ]; cur$lp[ok] <- nw$lp[ok]; cur$grad[ok, ] <- nw$grad[ok, , drop = FALSE]; h[ok] <- h[ok] * 1.5
      todo <- todo & !ok; h[todo] <- h[todo] / 2
      if (!any(todo)) break
    }
  }
  for (k in seq_along(vs)) s[[vs[k]]] <- if (is.finite(lo[k])) lo[k] + exp(x[, k]) else x[, k]
  s
}

# gradient of the objective by Fisher's identity: the expected complete-data score under the E-step
# weights (one E-step per call; run_em() differentiates it numerically for the Hessian)
bm_score <- function(E, C) function(s, d) {
  s <- s[E$fixed]; e <- bm_estep(E, s, C); out <- list()
  for (v in E$fixed) {
    P <- E$M$pars[[v]]; g <- numeric(dim_len(E$M, P$dim))
    for (f in factors_of(E$M, v)) {
      if (!bm_in_q(E, f) || (!f$observed && f$y == v && P$dist != "normal")) next
      h <- bm_where(E, f, s, e); M <- h$M; fx <- h$f; m <- fmask(fx, M); if (h$lat) m <- m * h$w
      x <- if (fx$dist == "normal") m * gnorm(fx, h$s, M, v, fval(fx, h$s, M, "r"), fval(fx, h$s, M, "sd2"))
           else m * (fval(fx, h$s, M, "y") - stats::plogis(fval(fx, h$s, M, "eta"))) * fval(fx, h$s, M, "d1", v)
      g <- g + bm_fold(E, reduce_to(x, P$dim, fx), P$dim, h$lat)
    }
    if (E$prior == "map" && P$dist != "normal") {                    # prior of a variance: numerical, it is cheap
      lp1 <- function(z) { s2 <- s; s2[[v]] <- z; bm_logprior_one(E, v, s2) }
      g <- g + vapply(seq_along(s[[v]]), function(k) { hh <- 1e-6 * s[[v]][k]; (lp1(replace(s[[v]], k, s[[v]][k] + hh)) - lp1(replace(s[[v]], k, s[[v]][k] - hh))) / (2 * hh) }, 0)
    }
    out[[v]] <- g
  }
  out
}

# Hessian of the objective by Louis' formula: E[complete-data Hessian] + Cov(complete-data score)
# under the E-step weights, from the first and second derivatives of the expressions (no numerical
# differentiation). Exact for the frozen grid; the exact latent variables get three Gauss-Hermite
# points (moments up to degree 5: the squared scores are quartic in them). A fixed parameter with one
# value per person (one not integrated out) enters only the parts of its own person: its complete-data
# score is 0 in the rows of the other persons, and the persons are independent given the fixed
# parameters, so the covariance of the scores is block diagonal by person. Its blocks are computed per
# person (sums over that person's nodes), not from a dense (persons x nodes) x persons score matrix.
bm_information <- function(E, C) {
  M <- E$M; N <- M$N
  if (length(E$ana)) {
    E$S <- as.matrix(expand.grid(rep(list(c(-sqrt(3), 0, sqrt(3))), length(E$ana))))
    E$Sw <- apply(as.matrix(expand.grid(rep(list(c(1, 4, 1) / 6), length(E$ana)))), 1, prod)
    E$ns <- nrow(E$S); E$Ms <- bm_stackM(M, E$ng * E$ns); E$sigma <- FALSE
  }
  facs <- Filter(function(f) f$dist %in% c("normal", "bernoulli_logit") && length(intersect(f$vars, E$fixed)) &&
                   bm_in_q(E, f), M$factors)
  per <- E$fixed[vapply(E$fixed, function(v) M$pars[[v]]$dim == "person", TRUE)]   # person-level parameters that are not integrated out
  function(s, d) {
    s <- s[E$fixed]; e <- bm_estep(E, s, C); nq <- E$ng * E$ns; pid <- rep(seq_len(N), nq)
    len <- vapply(E$fixed, function(v) length(s[[v]]), 0L); off <- stats::setNames(c(0L, cumsum(len))[seq_along(len)], E$fixed)
    at <- function(v, n) if (len[[v]] == 1) rep(off[[v]] + 1L, n) else off[[v]] + seq_len(len[[v]])
    ix <- function(v) off[[v]] + seq_len(len[[v]])
    gcol <- unlist(lapply(setdiff(E$fixed, per), ix)); gi <- integer(sum(len)); gi[gcol] <- seq_along(gcol)
    H <- matrix(0, sum(len), sum(len)); Sc <- matrix(0, N * nq, length(gcol)); Cp <- matrix(0, N * nq, length(per), dimnames = list(NULL, per))
    add <- function(i, j, x) { H[cbind(i, j)] <<- H[cbind(i, j)] + x }
    for (f in facs) {
      vs <- intersect(f$vars, E$fixed); g <- bm_where(E, f, s, e, collapse = FALSE); Mx <- g$M; fx <- g$f; sx <- g$s
      m <- fmask(fx, Mx); mw <- if (g$lat) m * g$w else m
      fold <- function(x) { x <- full_of(x, fx, Mx); switch(fx$shape, NK = colSums(x), person = sum(x), x) }
      byperson <- function(x) { x <- full_of(x, fx, Mx); if (g$lat) x <- rowsum(x, pid, reorder = FALSE); if (fx$shape == "NK") x else as.vector(x) }
      if (fx$dist == "normal") {
        r <- fval(fx, sx, Mx, "r"); s2 <- fval(fx, sx, Mx, "sd2")
        ru <- lapply(stats::setNames(vs, vs), function(u) fval(fx, sx, Mx, "d1", u))
        au <- lapply(stats::setNames(vs, vs), function(u) if (u %in% fx$sd_vars) bm_s2d(fx, sx, Mx, u) else 0)
        sc <- lapply(ru, function(x) NULL)
        for (u in vs) sc[[u]] <- -r * ru[[u]] / s2 + (r^2 / (2 * s2^2) - 1 / (2 * s2)) * au[[u]]
      } else {
        eta <- fval(fx, sx, Mx, "eta"); pr <- stats::plogis(eta); res <- fval(fx, sx, Mx, "y") - pr
        ru <- lapply(stats::setNames(vs, vs), function(u) fval(fx, sx, Mx, "d1", u))
        sc <- lapply(ru, function(x) res * x)
      }
      for (a in seq_along(vs)) for (b in seq_len(a)) {
        u <- vs[a]; v <- vs[b]; d2 <- bm_d2(fx, sx, Mx, u, v)
        h <- if (fx$dist == "bernoulli_logit") -pr * (1 - pr) * ru[[u]] * ru[[v]] + res * d2 else {
          auv <- if (u %in% fx$sd_vars && v %in% fx$sd_vars) bm_s2d(fx, sx, Mx, u, v) else 0
          -(ru[[u]] * ru[[v]] + r * d2) / s2 + r * (ru[[u]] * au[[v]] + ru[[v]] * au[[u]]) / s2^2 +
            r^2 * auv / (2 * s2^2) - r^2 * au[[u]] * au[[v]] / s2^3 - auv / (2 * s2) + au[[u]] * au[[v]] / (2 * s2^2)
        }
        if (!any(c(u, v) %in% per)) {                                  # item and model-level parameters
          x <- fold(mw * h)
          if (len[[u]] == 1 && len[[v]] == 1) x <- sum(x)
          iu <- at(u, length(x)); iv <- at(v, length(x))
        } else {                                                       # a person-level parameter: sums over each person's nodes
          xp <- byperson(mw * h); p <- if (u %in% per) u else v; q <- if (p == u) v else u
          if (q %in% per) { x <- if (is.matrix(xp)) rowSums(xp) else xp; iu <- ix(u); iv <- ix(v) }
          else if (M$pars[[q]]$dim == "item") { iu <- rep(ix(p), M$K); iv <- rep(ix(q), each = N); x <- as.vector(xp)
            if (p != u) { t0 <- iu; iu <- iv; iv <- t0 } }
          else { x <- if (is.matrix(xp)) rowSums(xp) else xp; iu <- if (p == u) ix(p) else rep(ix(q), N); iv <- if (p == u) rep(ix(q), N) else ix(p) }
        }
        add(iu, iv, x)
        if (u != v) add(iv, iu, x)
      }
      if (g$lat) for (u in vs) {                                        # complete-data scores on the nodes
        x <- full_of(m * sc[[u]], fx, Mx)
        if (u %in% per) Cp[, u] <- Cp[, u] + (if (fx$shape == "NK") rowSums(x) else x)
        else if (fx$shape == "NK" && len[[u]] > 1) Sc[, gi[at(u, M$K)]] <- Sc[, gi[at(u, M$K)]] + x
        else Sc[, gi[off[[u]] + 1L]] <- Sc[, gi[off[[u]] + 1L]] + (if (fx$shape == "NK") rowSums(x) else x)
      }
    }
    w <- as.vector(e$W); Sb <- rowsum(w * Sc, pid, reorder = FALSE)
    H[gcol, gcol] <- H[gcol, gcol] + crossprod(sqrt(w) * Sc) - crossprod(Sb)
    if (length(per)) {                                                 # covariance blocks of the person-level parameters, person by person
      Cb <- rowsum(w * Cp, pid, reorder = FALSE)
      for (k in seq_along(per)) {
        ip <- ix(per[k])
        GP <- rowsum(w * Cp[, k] * Sc, pid, reorder = FALSE) - Sb * Cb[, k]
        H[ip, gcol] <- H[ip, gcol] + GP; H[gcol, ip] <- H[gcol, ip] + t(GP)
        for (l in seq_len(k)) {
          iq <- ix(per[l]); x <- as.vector(rowsum(w * Cp[, k] * Cp[, l], pid, reorder = FALSE)) - Cb[, k] * Cb[, l]
          add(ip, iq, x); if (l != k) add(iq, ip, x)
        }
      }
    }
    if (E$prior == "map") for (v in E$fixed) if (M$pars[[v]]$dist != "normal") for (k in seq_len(len[[v]])) {
      lp <- function(z) { s2 <- s; s2[[v]][k] <- z; bm_logprior_one(E, v, s2) }        # prior of a variance: numerical, it is cheap
      x0 <- s[[v]][k]; hh <- 1e-4 * x0; i <- off[[v]] + k
      H[i, i] <- H[i, i] + (lp(x0 + hh) - 2 * lp(x0) + lp(x0 - hh)) / hh^2
    }
    if (!length(per) && all(vapply(Filter(function(f) f$observed, M$factors), function(f) any(f$vars %in% E$latent), TRUE))) {
      O <- matrix(0, sum(len), sum(len)); O[gcol, gcol] <- crossprod(Sb)   # the cross-product of the casewise scores (score_test() uses its rank)
      attr(H, "opg") <- O
    }
    H
  }
}
# second derivative of the main expression of f (eta, or r = y - mean) in u and v
bm_d2 <- function(f, s, M, u, v) {
  pk <- f$ir$d1[[u]]
  x <- if (!is.null(pk)) ml_eval(ml_pack(ml_deriv(pk$ir, v), M), s, M, f$shape) else eval(stats::D(f$d1[[u]], v), fenv(f, s, M))
  full_of(x, f, M)
}
# first (v = NULL) or second derivative of the variance of a normal factor
bm_s2d <- function(f, s, M, u, v = NULL) {
  pk <- f$ir$dsd2[[u]]
  if (!is.null(pk)) return(full_of(ml_eval(if (is.null(v)) pk else ml_pack(ml_deriv(pk$ir, v), M), s, M, f$shape), f, M))
  ev <- fenv(f, s, M); sd <- eval(f$sd, ev); du <- eval(f$dsd[[u]], ev)
  full_of(if (is.null(v)) 2 * sd * du else 2 * (eval(f$dsd[[v]], ev) * du + sd * eval(stats::D(f$dsd[[u]], v), ev)), f, M)
}

# one EM iteration for run_em(): E-step, then the steps in order with the same weights
bm_step <- function(E, C) function(s, d) {
  e <- bm_estep(E, s[E$fixed], C); s$.loglik <- e$obj
  e$omc <- new.env()                                                   # E[omega] of the logistic parts, shared within this M-step
  for (st in E$steps) {
    if (st$type == "gibbs") for (v in st$vars) s[[v]] <- if (E$M$pars[[v]]$dist == "normal") bm_normal_part(E, v, s[E$fixed], e) else bm_variance(E, v, s[E$fixed], e)
    if (st$type == "mala") s <- bm_gradient(E, st$vars, s, e)
  }
  s
}

#' @rdname ecm.block_model
#' @param x A fit of `ecm()`.
#' @export
print.block_model_ecm <- function(x, ...) {
  cat(sprintf("%s for the declared model: %s after %d iterations (%d rounds of the grid), %s %.4f\n",
              switch(x$ecm$method %||% "ecm", ecm = "ECM", newton = "Newton (nlminb)", `quasi-newton` = "Quasi-Newton (nlminb)"), if (x$converged) "converged" else "NOT converged",
              x$iterations, x$cycles, x$objective, x$loglik))
  print(format(utils::head(x$estimates, 30), digits = 3), row.names = FALSE)
  if (nrow(x$estimates) > 30) cat(sprintf("  ... %d more (fit$estimates)\n", nrow(x$estimates) - 30))
  if (all(is.na(x$estimates$se))) cat("se = NA: no standard errors (se = FALSE, not converged, or the Hessian is singular)\n")
  invisible(x)
}

# ---- check_ecm() --------------------------------------------------------------------------------------

#' Check an EM algorithm
#'
#' The EM counterpart of the Geweke check of [check_sampler()]. A Gibbs sampler is right when every
#' step leaves the posterior unchanged; an EM algorithm is right when
#'
#' 1. **every step never lowers the objective** (from many starting values, not only from the path
#'    that EM happened to take);
#' 2. **the E-step agrees with the objective** (Fisher's identity: the gradient of the log marginal
#'    likelihood equals the gradient of the expected complete-data log likelihood \eqn{Q(\cdot \mid x)}{Q(. | x)}
#'    at \eqn{x}{x}); and
#' 3. **ECM stops at a maximum** (a general optimizer started at the ECM estimate finds nothing better,
#'    and the gradient is zero there).
#'
#' For a [block_model()] all three are checked, with starting values drawn from the declared priors.
#' For a built-in normal model ([mlirt()], `rtirt_*()` without `quantile`) the three are checked for
#' the EM of [ecm()] with its own E- and M-step, from `n_start` starting values jittered around those
#' of [ecm()]: one whole M-step (which updates all blocks) never lowers the objective, the analytic
#' score (Fisher's identity) equals the numerical gradient of the objective, and the fit is a maximum.
#' For an EM written by hand, give `data`, `init`, `step` and `loglik` as for [run_em()]; checks 1
#' (from `n_start` jittered starts) and 3 are done.
#'
#' @param model A [block_model()], a built-in model (e.g. `rtirt_cross(data)`), or `NULL` for a
#'   hand-written EM.
#' @param steps,latent,nodes,prior,quadrature As in [ecm()] (for a built-in model only `prior` and
#'   `nodes`, default 41).
#' @param n_start Number of starting values.
#' @param data,init,step,loglik,positive A hand-written EM, as in [run_em()].
#' @param seed Random seed.
#' @param x An `em_check` object.
#' @param ... Unused.
#' @return A data frame with one row per check (`check`, `worst`, `limit`, `passed`), class `em_check`.
#' @examples
#' set.seed(1); N <- 200; K <- 5
#' th <- rnorm(N); b0 <- seq(-1, 1, length.out = K)
#' Y <- matrix(rbinom(N * K, 1, plogis(outer(th, b0, "-"))), N)
#' m <- block_model(Y = Y,
#'   theta[person] ~ normal(0, 1),
#'   a[item] ~ normal(1, 1, lower = 0),
#'   b[item] ~ normal(0, 3),
#'   Y ~ bernoulli_logit(a * (theta - b)))
#' check_ecm(m, list(step_gibbs(a, b)), n_start = 3, seed = 1)
#' # a built-in model: the EM of ecm()
#' cond <- set_cond(n_subj = 200, n_item = 5)
#' m_cross <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
#' check_ecm(m_cross, prior = "ml", n_start = 3, seed = 1)
#' @export
check_ecm <- function(model = NULL, steps = NULL, latent = NULL, nodes = NULL, prior = c("map", "ml"), n_start = 10,
                     data = NULL, init = NULL, step = NULL, loglik = NULL, positive = character(), seed = NULL, quadrature = c("adaptive", "ba81")) {
  if (!is.null(seed)) set.seed(seed)
  rows <- list()
  add <- function(check, worst, limit, detail) rows[[length(rows) + 1]] <<- data.frame(check = check, worst = worst, limit = limit,
                                                                                      passed = worst <= limit, detail = detail)
  if (inherits(model, "block_model")) {
    E <- bm_setup(model, steps, latent, nodes, match.arg(prior), NULL, quadrature = match.arg(quadrature))
    pk <- function(s) unlist(lapply(E$fixed, function(v) if (v %in% E$positive) log(s[[v]] - if (E$M$pars[[v]]$dist == "normal") E$M$pars[[v]]$lower else 0) else s[[v]]))
    up <- function(x, s) { i <- 0; for (v in E$fixed) { n <- length(s[[v]]); z <- x[i + seq_len(n)]; i <- i + n
      s[[v]] <- if (v %in% E$positive) (if (E$M$pars[[v]]$dist == "normal") E$M$pars[[v]]$lower else 0) + exp(z) else z }; s }
    qfun <- function(s2, e) bm_qvalue(E, s2, e)
    grad <- function(f, x) vapply(seq_along(x), function(k) { h <- 1e-5 * max(1, abs(x[k])); (f(replace(x, k, x[k] + h)) - f(replace(x, k, x[k] - h))) / (2 * h) }, 0)
    drop_q <- 0; drop_l <- 0; fisher <- 0; where <- ""
    for (r in seq_len(n_start)) {
      s <- simulate_model(E$M)[E$fixed]
      for (v in E$fixed) { P <- E$M$pars[[v]]                       # keep the starts away from the boundaries
        if (P$dist != "normal") s[[v]] <- pmin(pmax(s[[v]], 0.05), 20) else if (is.finite(P$lower)) s[[v]] <- pmax(s[[v]], P$lower + 0.2) }
      C <- bm_centres(E, s, NULL); e0 <- bm_estep(E, s, C); C <- bm_centres(E, s, e0); e0 <- bm_estep(E, s, C)
      q_prev <- qfun(s, e0); s1 <- s; k <- 0
      for (st in E$steps) {
        if (st$type == "gibbs") for (v in st$vars) {
          k <- k + 1
          s1[[v]] <- if (E$M$pars[[v]]$dist == "normal") bm_normal_part(E, v, s1, e0) else bm_variance(E, v, s1, e0)
          qn <- qfun(s1, e0); d <- (q_prev - qn) / max(1, abs(q_prev))
          if (d > drop_q) { drop_q <- d; where <- sprintf("step %d (%s)", k, v) }
          q_prev <- qn
        }
        if (st$type == "mala") { k <- k + 1; s1 <- bm_gradient(E, st$vars, s1, e0); qn <- qfun(s1, e0); d <- (q_prev - qn) / max(1, abs(q_prev))
          if (d > drop_q) { drop_q <- d; where <- sprintf("step %d (%s)", k, paste(st$vars, collapse = ", ")) }; q_prev <- qn }
      }
      drop_l <- max(drop_l, (e0$obj - bm_estep(E, s1, C)$obj) / max(1, abs(e0$obj)))
      x <- pk(s); g1 <- grad(function(z) bm_estep(E, up(z, s), C)$obj, x); g2 <- grad(function(z) qfun(up(z, s), e0), x)
      fisher <- max(fisher, max(abs(g1 - g2)) / max(1, max(abs(g1))))
    }
    add("each step raises Q(. | x) (relative drop)", drop_q, 1e-8, if (nzchar(where)) where else "")
    add("one iteration never lowers the objective (relative drop)", drop_l, 1e-8, "")
    add("Fisher identity: gradient of the objective = gradient of Q (relative difference)", fisher, 1e-4, "")
    fit <- suppressWarnings(ecm(model, steps, latent, nodes, E$prior, se = FALSE, quadrature = E$quadrature))
    s <- fit$state[E$fixed]; C <- fit$ecm$C; f <- function(z) bm_estep(E, up(z, s), C)$obj
    x <- pk(s); o <- stats::optim(x, f, method = "BFGS", control = list(fnscale = -1, maxit = 200))
    add("ECM stops at a maximum: gain of a general optimizer (relative)", (o$value - f(x)) / max(1, abs(f(x))), 1e-6, "")
    add("ECM stops at a maximum: largest |gradient|", max(abs(grad(f, x))), 1e-2, "")
  } else if (inherits(model, c("rtirt", "rtirt_qset"))) {
    em_check_builtin(model, match.arg(prior), nodes, n_start, add)
  } else {
    if (is.null(step) || is.null(loglik) || is.null(init)) stop("check_ecm(): give a block_model, a built-in model, or data, init, step and loglik as for run_em()", call. = FALSE)
    s0 <- if (is.function(init)) init(data) else init
    par <- names(s0)[!startsWith(names(s0), ".") & vapply(s0, is.numeric, TRUE)]
    drop_l <- 0; mism <- 0
    for (r in seq_len(n_start)) {
      s <- s0
      for (v in par) s[[v]] <- if (v %in% positive) s[[v]] * exp(stats::rnorm(length(s[[v]]), 0, 0.3)) else s[[v]] + stats::rnorm(length(s[[v]]), 0, 0.3)
      l0 <- loglik(s, data)
      for (it in 1:5) { s1 <- step(s, data); mism <- max(mism, abs(s1$.loglik - l0) / max(1, abs(l0))); s1$.loglik <- NULL
        l1 <- loglik(s1, data); drop_l <- max(drop_l, (l0 - l1) / max(1, abs(l0))); s <- s1; l0 <- l1 }
    }
    add("one iteration never lowers loglik() (relative drop)", drop_l, 1e-8, "")
    add(".loglik from step() = loglik() at the input state (relative difference)", mism, 1e-6, "")
    fit <- suppressWarnings(run_em(data, init, step, loglik = loglik, positive = positive, se = FALSE))
    pk <- function(s) unlist(lapply(par, function(v) if (v %in% positive) log(s[[v]]) else s[[v]]))
    skel <- fit$state[par]
    up_h <- function(x) { v <- utils::relist(x, skel); s <- fit$state; for (n in par) s[[n]] <- if (n %in% positive) exp(v[[n]]) else v[[n]]; s }
    f <- function(z) loglik(up_h(z), data); x <- pk(fit$state)
    o <- stats::optim(x, f, method = "BFGS", control = list(fnscale = -1, maxit = 200))
    g <- vapply(seq_along(x), function(k) { h <- 1e-5 * max(1, abs(x[k])); (f(replace(x, k, x[k] + h)) - f(replace(x, k, x[k] - h))) / (2 * h) }, 0)
    add("ECM stops at a maximum: gain of a general optimizer (relative)", (o$value - f(x)) / max(1, abs(f(x))), 1e-6, "")
    add("ECM stops at a maximum: largest |gradient|", max(abs(g)), 1e-2, "")
  }
  out <- do.call(rbind, rows)
  structure(out, class = c("em_check", "data.frame"), passed = all(out$passed), n_start = n_start)
}

#' @rdname check_ecm
#' @export
print.em_check <- function(x, ...) {
  cat(sprintf("%s (%d starting values): %s\n", attr(x, "what") %||% "ECM check", attr(x, "n_start"), if (attr(x, "passed")) "all checks passed" else "SOME CHECKS FAILED"))
  d <- as.data.frame(unclass(x)); d$worst <- signif(d$worst, 3)
  print(d, row.names = FALSE)
  invisible(x)
}
