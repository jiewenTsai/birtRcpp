# Marginal maximum likelihood by EM (Bock & Aitkin, 1981) for the normal models ----------
#
# E-step: ability on the Gauss-Hermite nodes of its N(mu_i, 1) distribution (the nodes of
# marginal_loglik()); given ability, speed is normal a posteriori (normal log RTs), so its
# moments are closed form and it is integrated analytically.
# M-step: the logistic part uses the Polya-Gamma bound (Polson, Scott & Windle, 2013): at the
# current eta, E[omega] = tanh(eta / 2) / (2 eta) gives a quadratic minorizer of
# y eta - log(1 + e^eta) (Jaakkola & Jordan, 2000), so each item is a 2 x 2 weighted least-
# squares problem and the marginal log-likelihood never decreases (an MM step inside EM).
# The RT and structural parameters have closed forms (a one-dimensional search for the
# correlation when the speed variance is fixed at 1). SQUAREM (Varadhan & Roland, 2008)
# accelerates the iterations, falling back to plain EM steps when it would lower the
# likelihood. Standard errors: observed information as the numerical Jacobian of the analytic
# score (Fisher identity; Louis, 1982) on an unconstrained scale, then the delta method.
#
# Internal parameters (natural scale): a, b, lambda, sigma2t; mlirt beta; rtirt beta1, beta2
# and zeta = x'beta2 + c (theta - x'beta1) + N(0, v); latent beta, gamma, s; cross rho, s.

#' Marginal maximum likelihood (ECM)
#'
#' Fits a built-in normal model by maximum marginal likelihood: the EM algorithm of Bock and
#' Aitkin (1981) with ability integrated by Gauss-Hermite quadrature and speed in closed form.
#' The item step uses the Polya-Gamma bound of the logistic likelihood (Polson, Scott &
#' Windle, 2013), so every step is closed form and the log-likelihood never decreases;
#' SQUAREM speeds up convergence. Standard errors come from the observed information (the
#' numerical derivative of the analytic score).
#'
#' **Quantile models.** With `quantile =` (models `"latent"` and `"cross"`, as in [gibbs()])
#' `ecm()` maximizes the marginal ALD working likelihood (Yu & Moyeed, 2001): the q-quantile of
#' the log times (`"cross"`) or of speed given ability and the covariates (`"latent"`) is linear
#' in the parameters. Ability is integrated on the nodes and speed in closed form, so the
#' working log-likelihood is smooth. The E-step and the item step are those of the normal models;
#' the response-time and structural parameters are updated by `nlminb()` with the analytic
#' gradient (two iterations per ECM step, which never lowers the working log-likelihood). The
#' fit starts at the ECM fit of the mean model, moved to the q-quantile. The standard errors are
#' sandwich (robust) standard errors, \eqn{H^{-1} J H^{-1}}{H^-1 J H^-1} with \eqn{J}{J} the
#' sum of the outer products of the scores of the persons: the ALD is not the distribution of the
#' data, so the observed information alone would be wrong. As in [gibbs()], `sigma2t` of
#' `"cross"` is the ALD scale of the log times and `var_speed` of `"latent"` the ALD scale of
#' speed. The working logLik, AIC and BIC compare only fits at the same quantile; [anova()]
#' refuses quantile fits. A vector of levels gives one fit per level (an `rtirt_qset`).
#'
#' The parameterization and identification are those of [gibbs()]: ability has mean
#' `x'beta` and variance 1, speed mean `x'beta` (no intercepts: the covariates are centred and
#' the item parameters carry the location). The fit starts the Gibbs sampler:
#'
#' ```
#' fit_ecm <- ecm(dat, model = "cross")
#' fit_gibbs <- gibbs(dat, model = "cross", init = fit_ecm)
#' ```
#'
#' **Small samples.** `print()` and `summary()` warn when the estimates suggest that the data
#' are too few for maximum likelihood: a difficulty with |b| > 6, a discrimination below 0.1
#' or above 10, or a standard error above 5 (or missing) for a, b or rho. These thresholds
#' flag only clearly broken fits (with ability N(0, 1), |b| > 6 means almost everyone or no
#' one solves the item, and an SE of 5 a 95% interval about 20 wide); `prior = "map"` or
#' [gibbs()] regularize such fits.
#'
#' @param data Data from [input_data()].
#' @param model The model by name (case does not matter):
#'   * `"irt"`: 2PL, no response times (with covariates: ability regressed on them).
#'   * `"null"`: responses and times; ability and speed correlated (van der Linden, 2007).
#'   * `"latreg"`: as `"null"`, with ability and speed regressed on the covariates.
#'   * `"latent"`: speed regressed on ability and the covariates.
#'   * `"cross"`: item-specific effects of ability on the log times (cross-relations).
#'
#'   A misspelled name gets the closest one in the error message. `"mixture"` (Wang & Xu,
#'   2015) is not available yet.
#' @param itemtype `"2pl"` or `"1pl"` (`"grm"` is not available yet); one value, or one per
#'   item (all equal for now).
#' @param speed_var `"free"` (default) or `"fixed"` (variance of speed 1; `"null"`,
#'   `"latreg"`, `"cross"`).
#' @param quantile `NULL` (the normal model), or one or more quantile levels in (0, 1) for
#'   `"latent"` and `"cross"` (see Quantile models).
#' @param cov Names of the covariates of `data` to use (default all); models without
#'   covariates ignore them.
#' @param prior `"ml"` (maximum likelihood) or `"map"`: `"map"` adds the item priors of
#'   `priors` (a, b, lambda, sigma2t and rho; by default those of [gibbs()], see
#'   [rtirt_priors()]): the posterior mode in the item parameters, for small samples. The
#'   structural parameters stay maximum likelihood.
#' @param priors Prior distributions for `prior = "map"`, from [rtirt_priors()].
#' @param se Compute standard errors.
#' @param control Further settings, from [ecm_control()], e.g.
#'   `ecm(dat, model = "cross", control = ecm_control(nodes = 61))`.
#' @param ... Not used: any other argument is an error before the fit starts. A setting of
#'   [ecm_control()] (`nodes`, `max_iter`, `tol`, `accelerate`, `verbose`) gets the
#'   `control = ecm_control(...)` form in the message, and a misspelled name the closest argument.
#' @return A fitted model with `ecm`: `estimates` (est, se, z, Wald 95% interval on the log scale
#'   for positive parameters and the Fisher-z scale for the correlation), `logLik`, `npar`,
#'   `AIC`, `BIC`, `iterations`, `converged`, `trace` (log-likelihood per iteration), `vcov`,
#'   `person` (EAP and posterior SD) and `secs`. Accessors: `print`, `summary`, `coef`,
#'   [estimates()], [scores()], [reliability()], [fit_indices()], `logLik`, `vcov`, `AIC`, `BIC`.
#' @references Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of
#'   item parameters: Application of an EM algorithm. *Psychometrika*, 46, 443-459.
#'
#'   Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for logistic models
#'   using Polya-Gamma latent variables. *JASA*, 108, 1339-1349.
#'
#'   White, H. (1982). Maximum likelihood estimation of misspecified models. *Econometrica*,
#'   50, 1-25.
#'
#'   Yu, K., & Moyeed, R. A. (2001). Bayesian quantile regression. *Statistics & Probability
#'   Letters*, 54, 437-447.
#'
#'   Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods for
#'   accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*,
#'   35, 335-353.
#' @examples
#' cond <- set_cond(n_subj = 300, n_item = 8)
#' dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' fit <- ecm(dat, model = "cross")
#' fit
#' head(estimates(fit))
#' \donttest{
#' fit_q <- ecm(dat, model = "cross", quantile = 0.9) # the 0.9-quantile of the log times
#' estimates(fit_q, pars = "^rho")
#' }
#' @export
ecm <- function(data, ...) UseMethod("ecm")

#' @export
ecm.default <- function(data, ...) {
  stop("data must come from input_data(), e.g. ecm(input_data(...), model = \"cross\")", call. = FALSE)
}

#' @export
ecm.rtirt_qset <- function(data, ...) {
  data$models <- lapply(data$models, ecm.rtirt, ...)
  data
}

#' Fit a model object by ECM (deprecated)
#'
#' The interface before `ecm(data, model = ...)`: [ecm()] on a model object from the
#' deprecated constructors ([models]).
#'
#' @param data A model object (see [models]).
#' @param max_iter,tol,nodes,accelerate,verbose See [ecm_control()].
#' @param se,prior,priors See [ecm()].
#' @param ... Unused.
#' @return The fitted model, as [ecm()].
#' @name ecm_model
#' @keywords internal
NULL

#' @rdname ecm_model
#' @export
ecm.rtirt <- function(data, max_iter = 1000, tol = 1e-10, nodes = 41, se = TRUE, prior = c("ml", "map"), priors = NULL, accelerate = TRUE,
                      verbose = FALSE, ...) {
  model <- data
  refuse_fit_as_data("ecm", model, names(list(...)))
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  t0 <- Sys.time()
  Q <- nodes
  map <- match.arg(prior) == "map"
  if (isTRUE(model$qr)) {
    if (map) stop("prior = \"map\" is for the normal models; quantile models are fitted by maximum (working) likelihood", call. = FALSE)
    return(em_qr_fit(model, max_iter, tol, Q, se, accelerate, verbose, t0))
  }
  priors <- priors %||% rtirt_priors()
  if (!inherits(priors, "rtirt_priors")) stop("priors must come from rtirt_priors()", call. = FALSE)
  E <- em_setup(model, Q, map, priors)
  x <- em_pack(em_start(E), E)
  F <- function(x) {
    p <- em_unpack(x, E)
    es <- em_estep(p, E)
    list(obj = em_obj(es, p, E), x = em_pack(em_mstep(p, E, es), E))
  }
  sq <- em_squarem(x, F,
    max_iter = max_iter, tol = tol, accelerate = accelerate,
    progress = if (verbose) function(it, obj) if (it %% 10 == 0) message(sprintf("ECM %d: logLik %.4f", it, obj))
  )
  x <- sq$x
  it <- sq$iterations
  conv <- sq$converged
  trace <- sq$trace
  p <- em_unpack(x, E)
  if (!conv) {
    # slow convergence usually means a parameter moving to the boundary of the space
    bd <- c(
      if (E$rt) sprintf("sigma2t[%s]", E$items[p$sigma2t < 1e-3]),
      if (!E$one_pl) sprintf("a[%s]", E$items[p$a > 10]),
      switch(E$eng,
        rtirt = if (!E$fixed && p$v < 1e-3) "var_speed (residual)" else if (abs(p$c) / sqrt(p$c^2 + p$v) > 0.99) "cor_ability_speed",
        latent = if (p$s < 1e-3) "var_speed",
        cross = if (!E$fixed && p$s < 1e-3) "var_speed",
        NULL
      )
    )
    warning(sprintf(
      "ecm() did not converge in %d iterations (raise max_iter)%s", max_iter,
      if (length(bd)) {
        sprintf(
          "; at the boundary: %s (a variance near 0 or a discrimination running off: the maximum is on the edge of the parameter space, so its standard errors are not valid; gibbs() or ecm(prior = \"map\") regularize)",
          paste(bd, collapse = ", ")
        )
      } else {
        ""
      }
    ), call. = FALSE)
  }
  es <- em_estep(p, E)
  ll <- sum(es$ll)
  npar <- length(x)
  N <- E$N
  V <- NULL
  if (se) {
    H <- em_hessian(x, E)
    V <- tryCatch(solve(-H), error = function(e) NULL)
    if (is.null(V) || any(!is.finite(diag(V))) || any(diag(V) < 0)) {
      warning("the observed information is not positive definite: no standard errors", call. = FALSE)
      V <- NULL
    }
  }
  mo <- em_moments(es)
  person <- list(theta = list(mean = mo$Et, sd = sqrt(pmax(mo$Et2 - mo$Et^2, 0))))
  if (E$rt) person$zeta <- list(mean = mo$Ez, sd = sqrt(pmax(mo$Ez2 - mo$Ez^2, 0)))
  tab <- em_table(x, V, E)
  vc <- attr(tab, "vcov")
  attr(tab, "vcov") <- NULL
  model$ecm <- list(
    estimates = tab, par = p, x = x, logLik = ll, npar = npar, nobs = N,
    AIC = if (map) NA_real_ else -2 * ll + 2 * npar, BIC = if (map) NA_real_ else -2 * ll + log(N) * npar, iterations = it, converged = conv,
    trace = trace, vcov = vc, person = person, map = map, priors = if (map) priors, nodes = Q,
    settings = model$settings, X = E$X, secs = as.numeric(Sys.time() - t0, units = "secs")
  )
  model
}

# model context: data, design, quadrature
em_setup <- function(object, Q, map, priors = rtirt_priors()) {
  type <- class(object)[1]
  info <- .model_info[[type]]
  S <- object$settings
  D <- object$data
  if (isTRUE(S$intercept)) stop("ecm() needs intercept = FALSE: without priors the intercepts of the latent regressions are not identified (the item parameters carry the location)", call. = FALSE)
  N <- nrow(D$Y)
  K <- ncol(D$Y)
  X <- if (info$X) D$X else matrix(0, N, 0)
  nd <- calc_nodes(Q)
  L <- D$log_t
  list(
    Y = D$Y, kappa = D$Y - 0.5, L = L, X = X, XX = crossprod(X), N = N, K = K, P = ncol(X), eng = info$engine, rt = info$rt,
    fixed = S$speed_var == "fixed" && info$engine %in% c("rtirt", "cross"), one_pl = S$itemtype == "1pl", map = map,
    z = nd$z, lw = nd$logw, items = D$items, xn = colnames(X) %||% character(),
    mu_l = if (info$rt) lambda_prior(priors, L)[1] else NA, sd_l = if (info$rt) lambda_prior(priors, L)[2] else NA, pr = priors
  )
}

em_start <- function(E) {
  pb <- pmin(pmax(colMeans(E$Y), 0.02), 0.98)
  p <- list(a = rep(1, E$K), b = -1.3 * stats::qlogis(pb))
  if (E$rt) {
    p$lambda <- colMeans(E$L)
    s2 <- var_floor(apply(E$L, 2, stats::var))
    vs <- max(stats::var(rowMeans(E$L)) - mean(s2) / E$K, 0.05)
    p$sigma2t <- pmax(s2 - vs, 0.05 * s2)
  }
  switch(E$eng,
    mlirt = p$beta <- rep(0, E$P),
    rtirt = {
      p$beta1 <- p$beta2 <- rep(0, E$P)
      p$c <- 0
      p$v <- if (E$fixed) 1 else vs
    },
    latent = {
      p$beta <- rep(0, E$P)
      p$gamma <- 0
      p$s <- vs
    },
    cross = {
      p$rho <- rep(0, E$K)
      p$s <- if (E$fixed) 1 else vs
    }
  )
  p
}

# natural-scale list <-> unconstrained vector (log for a and variances, atanh for a fixed-variance correlation)
em_pack <- function(p, E) {
  c(
    if (!E$one_pl) log(p$a), p$b, if (E$rt) c(p$lambda, log(p$sigma2t)),
    switch(E$eng,
      mlirt = p$beta,
      rtirt = c(p$beta1, p$beta2, if (E$fixed) atanh(p$c) else c(p$c, log(p$v))),
      latent = c(p$beta, p$gamma, log(p$s)),
      cross = c(p$rho, if (!E$fixed) log(p$s))
    )
  )
}

em_unpack <- function(x, E) {
  K <- E$K
  P <- E$P
  i <- 0
  take <- function(n) {
    v <- x[i + seq_len(n)]
    i <<- i + n
    v
  }
  p <- list(a = if (E$one_pl) rep(1, K) else exp(take(K)), b = take(K))
  if (E$rt) {
    p$lambda <- take(K)
    p$sigma2t <- exp(take(K))
  }
  switch(E$eng,
    mlirt = {
      p$beta <- take(P)
    },
    rtirt = {
      p$beta1 <- take(P)
      p$beta2 <- take(P)
      if (E$fixed) {
        p$c <- tanh(take(1))
        p$v <- 1 - p$c^2
      } else {
        p$c <- take(1)
        p$v <- exp(take(1))
      }
    },
    latent = {
      p$beta <- take(P)
      p$gamma <- take(1)
      p$s <- exp(take(1))
    },
    cross = {
      p$rho <- take(K)
      p$s <- if (E$fixed) 1 else exp(take(1))
    }
  )
  p
}

em_xb <- function(X, beta) if (length(beta)) drop(X %*% beta) else rep(0, nrow(X))
softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# E-step: nodes TH (N x Q), posterior weights W, posterior mean MZ (N x Q) and precision Pz of speed
em_estep <- function(p, E) {
  N <- E$N
  K <- E$K
  Q <- length(E$z)
  mu_t <- switch(E$eng,
    mlirt = em_xb(E$X, p$beta),
    rtirt = em_xb(E$X, p$beta1),
    rep(0, N)
  )
  TH <- outer(mu_t, E$z, "+")
  common <- all(mu_t == 0) # common nodes: Bock-Aitkin products; otherwise C++ (person-specific nodes)
  L <- calc_loglik_nodes(E$Y, p$a, p$b, if (common) E$z else TH) + matrix(E$lw, N, Q, byrow = TRUE)
  MZ <- NULL
  Pz <- NA
  if (E$rt) {
    C <- E$L - matrix(p$lambda, N, K, byrow = TRUE)
    w <- 1 / p$sigma2t
    rr <- if (E$eng == "cross") p$rho else rep(0, K)
    S1 <- drop(C %*% w)
    S2 <- drop(C^2 %*% w)
    Cr <- drop(C %*% (w * rr))
    wr <- sum(w * rr)
    wr2 <- sum(w * rr^2)
    M <- switch(E$eng,
      rtirt = em_xb(E$X, p$beta2) + p$c * (TH - mu_t),
      latent = em_xb(E$X, p$beta) + p$gamma * TH,
      0 * TH
    )
    # speed u ~ N(M, vz) seen through o_i = -(c_i + rho_i theta) = u + e_i: sums of w o and w o^2 at every node
    sp <- normal_given_sums(-(S1 + TH * wr), S2 + 2 * TH * Cr + TH^2 * wr2, w, M, if (E$eng == "rtirt") p$v else p$s)
    MZ <- sp$mean
    Pz <- 1 / sp$var
    L <- L + sp$logI
  }
  nw <- calc_node_weights(L)
  list(ll = nw$ll, TH = TH, W = nw$W, MZ = MZ, Pz = Pz, mu_t = mu_t, common = common)
}

em_moments <- function(es) {
  W <- es$W
  TH <- es$TH
  m <- list(Et = rowSums(W * TH), Et2 = rowSums(W * TH^2))
  if (!is.null(es$MZ)) {
    MZ <- es$MZ
    m$Ez <- rowSums(W * MZ)
    m$Ez2 <- rowSums(W * MZ^2) + 1 / es$Pz
    m$Etz <- rowSums(W * TH * MZ)
  }
  m
}

# item sums over persons and nodes: Polya-Gamma weights (MM step) or residuals y - p (score)
em_irt_sums <- function(p, E, es, score = FALSE) {
  calc_item_sums(E$Y, p$a, p$b, if (isTRUE(es$common)) E$z else es$TH, es$W, score)
}

em_obj <- function(es, p, E) sum(es$ll) + if (E$map) em_lprior(p, E) else 0

em_lprior <- function(p, E) {
  pr <- E$pr
  lp <- sum(stats::dnorm(p$b, pr$b[["mean"]], pr$b[["sd"]], log = TRUE)) +
    if (E$one_pl) 0 else sum(stats::dnorm(p$a, pr$a[1], pr$a[2], log = TRUE))
  if (E$rt) lp <- lp + sum(stats::dnorm(p$lambda, E$mu_l, E$sd_l, log = TRUE)) + sum(var_lprior(p$sigma2t, pr$sigma2t, pr$sigma2t_family)$lp)
  if (E$eng == "cross") lp <- lp + sum(stats::dnorm(p$rho, 0, pr$rho[["sd"]], log = TRUE))
  lp
}

em_mstep <- function(p, E, es) {
  mo <- em_moments(es)
  new <- em_mstep_items(p, E, es, mo)
  if (E$rt) {
    new <- em_mstep_rt(new, p, E, mo)
  }
  em_mstep_structural(new, E, mo)
}

# items: PG minorizer step (ML), or conditional modes with the Gibbs priors (MAP); eta = a theta - d
em_mstep_items <- function(p, E, es, mo) {
  new <- p
  sums <- em_irt_sums(p, E, es)
  T0 <- colSums(E$kappa)
  T1 <- colSums(E$kappa * mo$Et)
  if (!E$map && E$one_pl) {
    new$b <- (sums$S1 - T0) / sums$S0
    return(new)
  }
  prior_a <- if (E$map) unname(E$pr$a)
  prior_b <- if (E$map) unname(E$pr$b)
  update <- update_items_logistic(sums, T0, T1, p$a, p$b, prior_a = prior_a, prior_b = prior_b)
  new$b <- update$b
  if (!E$one_pl) {
    new$a <- update$a
  }
  new
}

# response times: log T + zeta = lambda - rho theta + e
em_mstep_rt <- function(new, p, E, mo) {
  N <- E$N
  K <- E$K
  Lt <- E$L
  if (E$eng == "cross") {
    lambda_rho <- em_cross_items(p, E, mo)
    new$lambda <- lambda_rho$lambda
    new$rho <- lambda_rho$rho
  } else if (!E$map) {
    new$lambda <- (colSums(Lt) + sum(mo$Ez)) / N
  } else {
    new$lambda <- ((colSums(Lt) + sum(mo$Ez)) / p$sigma2t + E$mu_l / E$sd_l^2) / (N / p$sigma2t + 1 / E$sd_l^2)
  }
  rho <- if (E$eng == "cross") new$rho else rep(0, K)
  SS <- em_rt_ss(Lt, new$lambda, rho, mo)
  new$sigma2t <- em_sigma2t(SS, E)
  new
}

# (lambda_i, rho_i) of the cross-relation model: a 2 x 2 weighted least-squares problem per item
em_cross_items <- function(p, E, mo) {
  N <- E$N
  K <- E$K
  Lt <- E$L
  w0 <- if (E$map) 1 / p$sigma2t else rep(1, K)
  n11 <- N * w0 + E$map / E$sd_l^2
  n12 <- -sum(mo$Et) * w0
  n22 <- sum(mo$Et2) * w0 + E$map / E$pr$rho[["sd"]]^2
  r1 <- (colSums(Lt) + sum(mo$Ez)) * w0 + E$map * E$mu_l / E$sd_l^2
  r2 <- -(drop(crossprod(Lt, mo$Et)) + sum(mo$Etz)) * w0
  determinant <- n11 * n22 - n12^2
  list(lambda = (n22 * r1 - n12 * r2) / determinant, rho = (n11 * r2 - n12 * r1) / determinant)
}

# residual variances: ML, or the conditional mode under the inverse-gamma or half-t prior (MAP)
em_sigma2t <- function(SS, E) {
  N <- E$N
  if (!E$map) {
    return(SS / N)
  }
  if (E$pr$sigma2t_family == "inverse-gamma") {
    return((E$pr$sigma2t[["scale"]] + SS / 2) / (E$pr$sigma2t[["shape"]] + N / 2 + 1))
  }
  log_posterior <- function(u, k) -N / 2 * u - SS[k] / 2 * exp(-u) + var_lprior(exp(u), E$pr$sigma2t, "half-t")$lp
  mode_k <- function(k) exp(stats::optimize(log_posterior, c(-15, 8), k = k, maximum = TRUE)$maximum)
  vapply(seq_len(E$K), mode_k, 0)
}

# structural part: latent regressions, speed variance
em_mstep_structural <- function(new, E, mo) {
  N <- E$N
  P <- E$P
  X <- E$X
  if (E$eng == "mlirt" && P) {
    new$beta <- drop(solve(E$XX, crossprod(X, mo$Et)))
  }
  if (E$eng %in% c("rtirt", "latent")) {
    new <- em_speed_regression(new, E, mo)
  }
  if (E$eng == "cross" && !E$fixed) {
    new$s <- max(sum(mo$Ez2) / N, 1e-8)
  }
  new
}

# regression of speed on (x, theta): coefficients and residual variance
em_speed_regression <- function(new, E, mo) {
  N <- E$N
  P <- E$P
  X <- E$X
  xt <- drop(crossprod(X, mo$Et))
  xz <- drop(crossprod(X, mo$Ez))
  G <- rbind(cbind(E$XX, xt), c(xt, sum(mo$Et2)))
  h <- c(xz, sum(mo$Etz))
  if (E$eng == "rtirt" && E$fixed) {
    # Var(zeta | x) = 1: c = correlation, v = 1 - c^2
    resid_ss <- function(c) {
      r <- xz - c * xt
      sum(mo$Ez2) - 2 * c * sum(mo$Etz) + c^2 * sum(mo$Et2) - if (P) sum(r * solve(E$XX, r)) else 0
    }
    profile <- function(c) -N / 2 * log(1 - c^2) - resid_ss(c) / (2 * (1 - c^2))
    c <- stats::optimize(profile, c(-0.9999, 0.9999), maximum = TRUE, tol = 1e-10)$maximum
    coefs <- c(if (P) solve(E$XX, xz - c * xt), c)
    v <- 1 - c^2
  } else {
    coefs <- drop(solve(G, h))
    v <- max((sum(mo$Ez2) - 2 * sum(coefs * h) + sum(coefs * (G %*% coefs))) / N, 1e-8)
  }
  if (E$eng == "rtirt") {
    new$beta1 <- if (P) drop(solve(E$XX, xt)) else numeric(0)
    new$c <- coefs[P + 1]
    new$v <- v
    new$beta2 <- coefs[seq_len(P)] + new$c * new$beta1
  } else {
    new$beta <- coefs[seq_len(P)]
    new$gamma <- coefs[P + 1]
    new$s <- v
  }
  new
}

# per item: sum_i E[(log T + zeta - lambda + rho theta)^2]
em_rt_ss <- function(Lt, lambda, rho, mo) {
  C <- sweep(Lt, 2, lambda)
  colSums(C^2) + 2 * colSums(C * mo$Ez) + 2 * rho * colSums(C * mo$Et) + sum(mo$Ez2) + 2 * rho * sum(mo$Etz) + rho^2 * sum(mo$Et2)
}

# score of the marginal log-likelihood (+ log prior) in the unconstrained parameters (Fisher identity)
em_score <- function(x, E) {
  p <- em_unpack(x, E)
  es <- em_estep(p, E)
  mo <- em_moments(es)
  N <- E$N
  K <- E$K
  X <- E$X
  ir <- em_irt_sums(p, E, es, score = TRUE)
  s_a <- ir$S1 - p$b * ir$S0
  s_b <- -p$a * ir$S0
  if (E$map) {
    pa <- E$pr$a
    s_a <- s_a - (p$a - pa[1]) / pa[2]^2
    s_b <- s_b - (p$b - E$pr$b[["mean"]]) / E$pr$b[["sd"]]^2
  }
  out <- c(if (!E$one_pl) p$a * s_a, s_b)
  if (E$rt) {
    rho <- if (E$eng == "cross") p$rho else rep(0, K)
    C <- sweep(E$L, 2, p$lambda)
    w <- 1 / p$sigma2t
    s_l <- (colSums(C) + sum(mo$Ez) + rho * sum(mo$Et)) * w
    SS <- em_rt_ss(E$L, p$lambda, rho, mo)
    s_s2 <- -N / 2 * w + SS / 2 * w^2
    if (E$map) {
      s_l <- s_l - (p$lambda - E$mu_l) / E$sd_l^2
      s_s2 <- s_s2 + var_lprior(p$sigma2t, E$pr$sigma2t, E$pr$sigma2t_family)$d
    }
    out <- c(out, s_l, p$sigma2t * s_s2)
  }
  str <- switch(E$eng,
    mlirt = drop(crossprod(X, mo$Et - em_xb(X, p$beta))),
    rtirt = {
      c <- p$c
      v <- p$v
      m1 <- em_xb(X, p$beta1)
      g <- em_xb(X, p$beta2) - c * m1
      Eu <- mo$Ez - c * mo$Et - g
      Eu2 <- mo$Ez2 - 2 * c * mo$Etz + c^2 * mo$Et2 - 2 * g * (mo$Ez - c * mo$Et) + g^2
      Eut <- mo$Etz - c * mo$Et2 - g * mo$Et
      s_b1 <- drop(crossprod(X, mo$Et - m1)) - c / v * drop(crossprod(X, Eu))
      s_b2 <- drop(crossprod(X, Eu)) / v
      s_c <- sum(Eut - m1 * Eu) / v
      s_v <- sum(-1 / (2 * v) + Eu2 / (2 * v^2))
      c(s_b1, s_b2, if (E$fixed) (1 - c^2) * (s_c - 2 * c * s_v) else c(s_c, v * s_v))
    },
    latent = {
      gm <- p$gamma
      s <- p$s
      g <- em_xb(X, p$beta)
      Eu <- mo$Ez - gm * mo$Et - g
      Eu2 <- mo$Ez2 - 2 * gm * mo$Etz + gm^2 * mo$Et2 - 2 * g * (mo$Ez - gm * mo$Et) + g^2
      Eut <- mo$Etz - gm * mo$Et2 - g * mo$Et
      c(drop(crossprod(X, Eu)) / s, sum(Eut) / s, s * sum(-1 / (2 * s) + Eu2 / (2 * s^2)))
    },
    cross = {
      C <- sweep(E$L, 2, p$lambda)
      s_r <- -(colSums(C * mo$Et) + sum(mo$Etz) + p$rho * sum(mo$Et2)) / p$sigma2t
      if (E$map) s_r <- s_r - p$rho / E$pr$rho[["sd"]]^2
      c(s_r, if (!E$fixed) p$s * sum(-1 / (2 * p$s) + mo$Ez2 / (2 * p$s^2)))
    }
  )
  c(out, str)
}

# observed information: central differences of the analytic score
em_hessian <- function(x, E) calc_hessian(function(z) em_score(z, E), x)

# reported parameters (the names of gibbs()) as a function of the unconstrained vector
em_report <- function(x, E) {
  p <- em_unpack(x, E)
  v <- c(
    p$a, p$b, if (E$rt) c(p$lambda, p$sigma2t),
    switch(E$eng,
      mlirt = p$beta,
      rtirt = c(p$beta1, p$beta2, p$c / sqrt(p$c^2 + p$v), p$c^2 + p$v),
      latent = c(p$beta, p$gamma, p$s),
      cross = c(p$rho, p$s)
    )
  )
  stats::setNames(v, par_names(E$eng, E$items, E$xn))
}

em_table <- function(x, V, E) {
  est <- em_report(x, E)
  n <- length(est)
  fixed <- rep(FALSE, n)
  nm <- names(est)
  if (E$one_pl) fixed[startsWith(nm, "a[")] <- TRUE
  if (E$fixed) fixed[nm == "var_speed"] <- TRUE
  se <- rep(NA_real_, n)
  Vr <- NULL
  if (!is.null(V)) {
    J <- vapply(seq_along(x), function(j) {
      h <- 1e-6 * max(1, abs(x[j]))
      e <- replace(numeric(length(x)), j, h)
      (em_report(x + e, E) - em_report(x - e, E)) / (2 * h)
    }, numeric(n))
    J <- matrix(J, n)
    Vr <- J %*% V %*% t(J)
    dimnames(Vr) <- list(nm, nm)
    se <- sqrt(pmax(diag(Vr), 0))
    se[fixed] <- NA
  }
  pos <- grepl("^(a|sigma2t)\\[", nm) | nm == "var_speed"
  cor <- nm == "cor_ability_speed"
  lo <- est - 1.959964 * se
  hi <- est + 1.959964 * se
  lo[pos] <- est[pos] * exp(-1.959964 * se[pos] / est[pos])
  hi[pos] <- est[pos] * exp(1.959964 * se[pos] / est[pos])
  if (any(cor)) {
    z <- atanh(est[cor])
    sz <- se[cor] / (1 - est[cor]^2)
    lo[cor] <- tanh(z - 1.959964 * sz)
    hi[cor] <- tanh(z + 1.959964 * sz)
  }
  zv <- est / se
  zv[!zero_null(nm)] <- NA # as sig: no Wald test where 0 is not a meaningful value
  tab <- data.frame(
    parameter = nm, est = unname(est), se = unname(se), z = unname(zv), q025 = unname(lo), q975 = unname(hi),
    sig = sig_flag(nm, lo, hi), row.names = NULL
  )
  attr(tab, "vcov") <- Vr
  tab
}

# starting values of a Gibbs chain from the ECM fit (jittered)
em_inits <- function(model, engine, X, fixed) {
  f <- model$ecm
  p <- f$par
  N <- nrow(model$data$Y)
  K <- ncol(model$data$Y)
  jit <- function(v, s = 0.05) v + stats::rnorm(length(v), 0, s)
  pad <- function(b) if (ncol(X) > length(b)) c(rep(0, ncol(X) - length(b)), b) else b # intercept column first
  init <- list(theta = jit(f$person$theta$mean, 0.1), a = p$a * exp(stats::rnorm(K, 0, 0.05)), b = jit(p$b))
  if (engine != "mlirt") {
    init <- c(init, list(
      zeta = jit(f$person$zeta$mean, 0.05), lambda = jit(p$lambda, 0.02),
      sigma2 = p$sigma2t * exp(stats::rnorm(K, 0, 0.05))
    ))
  }
  switch(engine,
    mlirt = c(init, list(beta = pad(jit(p$beta)))),
    rtirt = c(init, list(
      beta1 = pad(jit(p$beta1)), beta2 = pad(jit(p$beta2)),
      c = if (fixed) tanh(atanh(p$c) + stats::rnorm(1, 0, 0.05)) else jit(p$c), v = if (fixed) 1 else p$v
    )),
    latent = c(init, list(beta = c(pad(jit(p$beta)), jit(p$gamma)), s = p$s * exp(stats::rnorm(1, 0, 0.05)))),
    cross = c(init, list(rho = jit(p$rho), s = if (fixed) 1 else p$s * exp(stats::rnorm(1, 0, 0.05))))
  )
}

# check_ecm() for a built-in model, with the E- and M-step of ecm(): one M-step never lowers the objective,
# Fisher's identity (em_score() = numerical gradient of the objective), and the fit is a maximum
em_check_builtin <- function(model, prior, nodes, n_start, add) {
  if (inherits(model, "rtirt_qset")) {
    stop("check_ecm() checks one model: check one quantile level, e.g. check_ecm(fits$models[[1]])", call. = FALSE)
  }
  if (isTRUE(model$qr)) {
    return(em_qr_check(model, nodes, n_start, add))
  }
  Q <- nodes %||% 41
  E <- em_setup(model, Q, prior == "map")
  obj <- function(x) {
    p <- em_unpack(x, E)
    em_obj(em_estep(p, E), p, E)
  }
  grad <- function(f, x) {
    vapply(seq_along(x), function(k) {
      h <- 1e-5 * max(1, abs(x[k]))
      (f(replace(x, k, x[k] + h)) - f(replace(x, k, x[k] - h))) / (2 * h)
    }, 0)
  }
  xh <- suppressWarnings(ecm(model, nodes = Q, prior = prior, se = FALSE))$ecm$x
  x0 <- em_pack(em_start(E), E)
  drop_l <- 0
  fisher <- 0
  for (r in seq_len(n_start)) {
    x <- x0 + stats::rnorm(length(x0), 0, 0.3) # jittered on the unconstrained scale
    p <- em_unpack(x, E)
    es <- em_estep(p, E)
    o0 <- em_obj(es, p, E)
    drop_l <- max(drop_l, (o0 - obj(em_pack(em_mstep(p, E, es), E))) / max(1, abs(o0)))
    g <- grad(obj, x)
    fisher <- max(fisher, max(abs(g - em_score(x, E))) / max(1, max(abs(g))))
  }
  add("one iteration never lowers the objective (relative drop)", drop_l, 1e-8, "em_mstep(): all blocks")
  add("Fisher identity: gradient of the objective = gradient of Q (relative difference)", fisher, 1e-4, "em_score()")
  x <- xh
  o <- stats::optim(x, obj, method = "BFGS", control = list(fnscale = -1, maxit = 200))
  add("ECM stops at a maximum: gain of a general optimizer (relative)", (o$value - obj(x)) / max(1, abs(obj(x))), 1e-6, "")
  add("ECM stops at a maximum: largest |gradient|", max(abs(grad(obj, x))), 1e-2, "")
}

# ---- methods for ECM fits --------------------------------------------------------------------

has_em <- function(x) is.null(x$post) && !is.null(x$ecm)
need_em <- function(x) if (is.null(x$ecm)) stop("run ecm() first", call. = FALSE)

print_em <- function(x) {
  f <- x$ecm
  quantile <- isTRUE(x$qr)
  objective <- if (f$map) {
    "marginal posterior mode (item priors of gibbs())"
  } else if (quantile) {
    "maximum marginal ALD working likelihood"
  } else {
    "maximum marginal likelihood"
  }
  cat(sprintf(
    "  %s (ECM, %d Gauss-Hermite nodes): %d iterations%s, %.1f s\n",
    objective, f$nodes, f$iterations,
    if (f$converged) "" else " (not converged)", f$secs
  ))
  if (quantile) {
    cat(sprintf("  working logLik %s, npar %d, AIC %s, BIC %s (compare only fits at the same quantile)\n",
                fmt_fixed(f$logLik), f$npar, fmt_fixed(f$AIC), fmt_fixed(f$BIC)))
    if (!is.null(f$vcov)) cat("  Standard errors: sandwich (robust to the ALD being a working likelihood)\n")
  } else if (f$map) {
    cat(sprintf("  log posterior %s, npar %d (no AIC / BIC: the objective includes the priors)\n", fmt_fixed(f$logLik), f$npar))
  } else {
    cat(sprintf("  logLik %s, npar %d, AIC %s, BIC %s\n", fmt_fixed(f$logLik), f$npar, fmt_fixed(f$AIC), fmt_fixed(f$BIC)))
  }
  rel <- reliability(x)
  cat(sprintf("  Reliability: %s\n", paste(sprintf("%s %.3f", names(rel), rel), collapse = ", ")))
  extreme <- ecm_extreme_warning(x)
  if (length(extreme)) {
    cat(extreme, "\n", sep = "")
  }
  cat("Use summary(), coef(), scores(), estimates(), fit_indices(); gibbs(data, model = ..., init = fit) to sample from here.\n")
}

# Estimates that say the data are too few for maximum likelihood (the N = 40 fit with b = -34.1,
# SE 4105). Thresholds, chosen to flag only the clearly broken cases: |b| > 6 (with ability
# N(0, 1) almost everyone, or no one, solves the item), a < 0.1 or a > 10 (an item without
# information, or a step function), and a standard error above 5 (or none) for a, b or rho
# (a 95% interval wider than 20 on the logit or correlation scale). Empty when none.
ecm_extreme_warning <- function(x) {
  e <- x$ecm$estimates
  name <- e$parameter
  is_a <- startsWith(name, "a[") & !identical(x$settings$itemtype, "1pl")
  is_b <- startsWith(name, "b[")
  extreme <- (is_b & abs(e$est) > 6) | (is_a & (e$est < 0.1 | e$est > 10))
  has_se <- !is.null(x$ecm$vcov)
  with_se <- is_a | is_b | startsWith(name, "rho[")
  wild_se <- has_se & with_se & (is.na(e$se) | e$se > 5)
  flagged <- which(extreme | wild_se)
  if (!length(flagged)) {
    return(character())
  }
  se_text <- if (has_se) sprintf(", SE %s", format_number(e$se[flagged])) else ""
  shown <- sprintf("%s = %s%s", name[flagged], format_number(e$est[flagged]), se_text)
  more <- if (length(shown) > 3) sprintf(" and %d more", length(shown) - 3) else ""
  remedy <- if (isTRUE(x$ecm$map)) "gibbs() or fewer items" else "ecm(..., prior = \"map\") or gibbs()"
  sprintf("  Warning: extreme estimates or standard errors (%s%s): the sample may be too small for maximum likelihood; try %s",
          paste(utils::head(shown, 3), collapse = "; "), more, remedy)
}

format_number <- function(v) {
  fmt_fixed(v, na = "NA")
}

#' @export
logLik.rtirt <- function(object, ...) {
  need_em(object)
  structure(object$ecm$logLik, df = object$ecm$npar, nobs = object$ecm$nobs, class = "logLik")
}

#' Compare fits of the same data
#'
#' `anova(fit0, fit1, ...)` compares [ecm()] fits of the same responses and response times,
#' as `anova()` of `lm()`: the number of parameters, log-likelihood, AIC and BIC of each fit,
#' and the likelihood ratio test of each fit against the one in the row above
#' (`Chisq = 2 * (logLik - logLik above)`, `Df` the difference in the number of parameters).
#' The rows are ordered by the number of parameters. The test assumes that each fit is nested
#' in the next, e.g. `"null"` in `"latreg"` (the regression coefficients 0) or in `"cross"`
#' (equal rho), `itemtype = "1pl"` in `"2pl"`, `speed_var = "fixed"` in `"free"`; AIC and BIC
#' also compare fits that are not nested. A fit with more parameters but a lower
#' log-likelihood gives a warning (not nested, or not converged).
#'
#' For [gibbs()] fits the table has the criteria of [fit_indices()] (DIC, WAIC, LOOIC with its
#' standard error), smaller is better, and no test. Fits with `prior = "map"` (no likelihood),
#' quantile fits (an ALD working likelihood) and a mix of ECM and Gibbs fits are refused.
#' @param object,... Fits of [ecm()] (or of [gibbs()]) of the same data.
#' @return A data frame of class `rtirt_anova` (and `anova`), printed with a heading that names
#'   the fits.
#' @examples
#' cond <- set_cond(n_subj = 300, n_item = 8, n_feat = 1)
#' dat <- sim_data(cond, sim_para(cond, "latreg"), "latreg")
#' fit0 <- ecm(dat, model = "null")
#' fit1 <- ecm(dat, model = "latreg")
#' anova(fit0, fit1)
#' @export
anova.rtirt <- function(object, ...) {
  fits <- list(object, ...)
  labels <- vapply(as.list(substitute(list(object, ...)))[-1], function(e) paste(deparse(e, width.cutoff = 60L)[1], collapse = ""), "")
  if (length(fits) < 2) {
    stop("anova() compares two or more fits of the same data, e.g. anova(fit0, fit1)", call. = FALSE)
  }
  is_fit <- vapply(fits, function(f) inherits(f, "rtirt") && (has_em(f) || !is.null(f$post)), TRUE)
  if (!all(is_fit)) {
    stop(sprintf("anova() compares fits of ecm() or gibbs(); not a fit: %s", paste(labels[!is_fit], collapse = ", ")), call. = FALSE)
  }
  is_ecm <- vapply(fits, has_em, TRUE)
  if (any(is_ecm) && !all(is_ecm)) {
    stop("anova() compares ecm() fits with ecm() fits (likelihood ratio, AIC, BIC) or gibbs() fits with gibbs() fits (WAIC, LOOIC), not both",
         call. = FALSE)
  }
  check_same_data(fits, labels)
  quantile_fit <- vapply(fits, function(f) isTRUE(f$qr), TRUE)
  if (any(quantile_fit)) {
    stop(sprintf("%s: quantile fits have an ALD working likelihood, not comparable with other fits", paste(labels[quantile_fit], collapse = ", ")),
         call. = FALSE)
  }
  if (all(is_ecm)) {
    order_npar <- order(vapply(fits, function(f) as.numeric(f$ecm$npar), 0))
    fits <- fits[order_npar]
    labels <- labels[order_npar]
  }
  heading_models <- paste(sprintf("%s: %s", labels, vapply(fits, model_label, "")), collapse = "\n")
  if (all(is_ecm)) anova_ecm(fits, labels, heading_models) else anova_gibbs(fits, labels, heading_models)
}

# the fits of anova() must have the same responses and response times
check_same_data <- function(fits, labels) {
  first <- fits[[1]]$data
  has_rt <- vapply(fits, function(f) .model_info[[class(f)[1]]]$rt, TRUE)
  for (k in seq_along(fits)[-1]) {
    data <- fits[[k]]$data
    if (has_rt[1] != has_rt[k]) {
      stop(sprintf("%s and %s are not comparable: one has response times and the other not (model = \"irt\")", labels[1],
                   labels[k]), call. = FALSE)
    }
    if (!identical(dim(first$Y), dim(data$Y)) || !identical(first$Y, data$Y) || !identical(first$log_t, data$log_t)) {
      stop(sprintf("%s and %s are fits of different data; anova() compares fits of the same responses and times", labels[1],
                   labels[k]), call. = FALSE)
    }
  }
}

anova_ecm <- function(fits, labels, heading_models) {
  map <- vapply(fits, function(f) isTRUE(f$ecm$map), TRUE)
  if (any(map)) {
    stop(sprintf("%s: fitted with prior = \"map\", so there is no likelihood to compare; refit with prior = \"ml\"",
                 paste(labels[map], collapse = ", ")), call. = FALSE)
  }
  not_converged <- !vapply(fits, function(f) isTRUE(f$ecm$converged), TRUE)
  if (any(not_converged)) {
    warning(sprintf("%s did not converge; the comparison may be off", paste(labels[not_converged], collapse = ", ")), call. = FALSE)
  }
  npar <- vapply(fits, function(f) as.numeric(f$ecm$npar), 0)
  loglik <- vapply(fits, function(f) f$ecm$logLik, 0)
  chisq <- c(NA, 2 * diff(loglik))
  df <- c(NA, diff(npar))
  testable <- !is.na(df) & df > 0 & chisq > -1e-6
  p_value <- ifelse(testable, stats::pchisq(pmax(chisq, 0), df, lower.tail = FALSE), NA_real_)
  lower <- which(!is.na(chisq) & chisq < -1e-6 & df > 0)
  if (length(lower)) {
    warning(sprintf("%s has more parameters than %s but a lower log-likelihood: the fits are not nested (or did not converge); compare them by AIC or BIC",
                    labels[lower[1]], labels[lower[1] - 1]), call. = FALSE)
  }
  table <- data.frame(npar = npar, logLik = loglik, AIC = vapply(fits, function(f) f$ecm$AIC, 0),
                      BIC = vapply(fits, function(f) f$ecm$BIC, 0), Chisq = chisq, Df = df, `Pr(>Chisq)` = p_value,
                      row.names = labels, check.names = FALSE)
  heading <- c(paste0("Likelihood ratio tests of ecm() fits (each row against the row above; valid for nested fits)\n\n",
                      heading_models, "\n"))
  structure(table, heading = heading, class = c("rtirt_anova", "anova", "data.frame"))
}

#' @rdname anova.rtirt
#' @param x The result of `anova()`.
#' @param digits Decimals (log-likelihood, AIC, BIC, Chisq, p-value; `npar` and `Df` are counts).
#' @param signif.stars Significance stars after the p-value, as in [stats::printCoefmat()].
#' @export
print.rtirt_anova <- function(x, digits = 3, signif.stars = getOption("show.signif.stars", TRUE), ...) {
  cat(attr(x, "heading"), sep = "\n")
  table <- data.frame(model = rownames(x), stringsAsFactors = FALSE)
  for (column in names(x)) {
    table[[column]] <- x[[column]]
  }
  names(table)[1] <- ""
  p <- x[["Pr(>Chisq)"]]
  stars <- !is.null(p) && isTRUE(signif.stars) && any(!is.na(p))
  if (stars) {
    codes <- stats::symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " "))
    table[[" "]] <- ifelse(is.na(p), "", as.character(codes))
  }
  print_table(table, digits)
  if (stars) {
    cat("---\nSignif. codes:  ", attr(codes, "legend"), "\n", sep = "")
  }
  invisible(x)
}

anova_gibbs <- function(fits, labels, heading_models) {
  indices <- lapply(fits, function(f) fit_indices(f))
  get <- function(name) vapply(indices, function(x) as.numeric(x[[name]]), 0)
  table <- data.frame(DIC = get("DIC"), pD = get("pD"), WAIC = get("WAIC"), LOOIC = get("LOOIC"), SE_LOOIC = get("SE_LOOIC"),
                      row.names = labels, check.names = FALSE)
  heading <- paste0("Fit criteria of gibbs() fits (marginal likelihood; smaller is better)\n\n", heading_models, "\n")
  structure(table, heading = heading, class = c("rtirt_anova", "anova", "data.frame"))
}

#' @export
vcov.rtirt <- function(object, ...) {
  need_em(object)
  if (is.null(object$ecm$vcov)) stop("no standard errors: ecm(..., se = TRUE)", call. = FALSE)
  object$ecm$vcov
}

summary_em <- function(object, digits) {
  e <- object$ecm$estimates
  if (isTRUE(object$ecm$map)) { # item priors only; structural parameters by ML
    pr <- prior_strings(object)
    item <- grepl("^(a|b|lambda|sigma2t)\\[", e$parameter)
    e$prior <- ifelse(item, unname(pr[e$parameter]), "(ML)")
  }
  quantile <- isTRUE(object$qr)
  estimate <- if (isTRUE(object$ecm$map)) {
    "posterior mode under the item priors"
  } else if (quantile) {
    "ALD working-likelihood estimate"
  } else {
    "ML estimate"
  }
  p <- rtirt_partable(e, FALSE, quantile && class(object)[1] == "rtirt_cross", sprintf(
    "%s, %sSE, Wald z and 95%% interval (log scale for a, sigma2t, var_speed; Fisher-z for the correlation)",
    estimate, if (quantile) "sandwich " else ""
  ))
  w <- rtirt_item_wide(p, object$data$items)
  structure(
    list(
      header = c(rtirt_header(object), ecm_extreme_warning(object)),
      fit = rtable(fit_indices(object), if (isTRUE(object$qr)) {
        "marginal ALD working likelihood (persons integrated out; compare only fits at the same quantile)"
      } else {
        "marginal likelihood (persons integrated out)"
      }),
      items = rtable(w$est, "item parameters (ECM estimates)"),
      se = rtable(w$se, "standard errors of the item parameters"),
      parameters = p,
      reliability = rtable(as.data.frame(as.list(unclass(reliability(object)))), "reliability of the EAP scores (final E-step)")
    ),
    class = c("summary.rtirt", "rtirt_summary")
  )
}
