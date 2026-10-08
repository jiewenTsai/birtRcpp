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
#' Fits the normal models ([mlirt()], [rtirt_null()], [rtirt_latreg()], [rtirt_latent()] and
#' [rtirt_cross()] without `quantile`) by maximum marginal likelihood: the EM algorithm of
#' Bock and Aitkin (1981) with ability integrated by Gauss-Hermite quadrature and speed in
#' closed form. The item step uses the Polya-Gamma bound of the logistic likelihood (Polson,
#' Scott & Windle, 2013), so every step is closed form and the log-likelihood never decreases;
#' SQUAREM speeds up convergence. Standard errors come from the observed information (the
#' numerical derivative of the analytic score). Quantile (ALD) models use [gibbs()].
#'
#' The parameterization and identification are those of [gibbs()]: ability has mean
#' `x'beta` and variance 1, speed mean `x'beta` (no intercepts: the covariates are centred and
#' the item parameters carry the location). The same object can then be sampled, e.g. with
#' starting values from the ECM fit: `gibbs(fit, init = "ecm")`.
#'
#' @param model A model object (see [models]), or a [block_model()] (see [ecm.block_model()]).
#' @param max_iter Maximum number of ECM (SQUAREM) iterations.
#' @param tol Convergence: relative change of the log-likelihood below `tol` and largest
#'   parameter change (unconstrained scale) below `sqrt(tol)`.
#' @param nodes Gauss-Hermite nodes for ability.
#' @param se Compute standard errors.
#' @param prior `"ml"` (maximum likelihood) or `"map"`: `"map"` adds the item priors of `priors` (a, b, lambda, sigma2t and rho; by default
#'   those of [gibbs()], see [rtirt_priors()]): the posterior mode in the item parameters, for small
#'   samples. The structural parameters stay maximum likelihood.
#' @param priors Prior distributions for `prior = "map"`, from [rtirt_priors()].
#' @param accelerate Use SQUAREM.
#' @param verbose Print the log-likelihood every 10 iterations.
#' @param ... Unused.
#' @return The model with `ecm`: `estimates` (est, se, z, Wald 95% interval on the log scale
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
#'   Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods for
#'   accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*,
#'   35, 335-353.
#' @examples
#' cond <- set_cond(n_subj = 300, n_item = 8)
#' dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' fit <- ecm(rtirt_cross(dat))
#' fit
#' head(estimates(fit))
#' @export
ecm <- function(model, ...) UseMethod("ecm")

#' @export
ecm.default <- function(model, ...) stop("model must be a model object, e.g. rtirt_cross(data), or a block_model()", call. = FALSE)

#' @export
ecm.rtirt_qset <- function(model, ...) stop("ecm() fits the normal models; quantile models use gibbs()", call. = FALSE)

#' @rdname ecm
#' @export
ecm.rtirt <- function(model, max_iter = 1000, tol = 1e-10, nodes = 41, se = TRUE, prior = c("ml", "map"), priors = NULL, accelerate = TRUE,
                     verbose = FALSE, ...) {
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  if (isTRUE(model$qr)) stop("ecm() fits the normal models; quantile models use gibbs()", call. = FALSE)
  t0 <- Sys.time(); Q <- nodes; map <- match.arg(prior) == "map"
  priors <- priors %||% rtirt_priors()
  if (!inherits(priors, "rtirt_priors")) stop("priors must come from rtirt_priors()", call. = FALSE)
  E <- em_setup(model, Q, map, priors)
  x <- em_pack(em_start(E), E)
  F <- function(x) { p <- em_unpack(x, E); es <- em_estep(p, E); list(obj = em_obj(es, p, E), x = em_pack(em_mstep(p, E, es), E)) }
  sq <- em_squarem(x, F, max_iter = max_iter, tol = tol, accelerate = accelerate,
                   progress = if (verbose) function(it, obj) if (it %% 10 == 0) message(sprintf("ECM %d: logLik %.4f", it, obj)))
  x <- sq$x; it <- sq$iterations; conv <- sq$converged; trace <- sq$trace
  p <- em_unpack(x, E)
  if (!conv) {
    # slow convergence usually means a parameter moving to the boundary of the space
    bd <- c(if (E$rt) sprintf("sigma2t[%s]", E$items[p$sigma2t < 1e-3]),
            if (!E$one_pl) sprintf("a[%s]", E$items[p$a > 10]),
            switch(E$eng, rtirt = if (!E$fixed && p$v < 1e-3) "var_speed (residual)" else if (abs(p$c) / sqrt(p$c^2 + p$v) > 0.99) "cor_ability_speed",
                   latent = if (p$s < 1e-3) "var_speed", cross = if (!E$fixed && p$s < 1e-3) "var_speed", NULL))
    warning(sprintf("ecm() did not converge in %d iterations (raise max_iter)%s", max_iter,
                    if (length(bd)) sprintf("; at the boundary: %s (a variance near 0 or a discrimination running off: the maximum is on the edge of the parameter space, so its standard errors are not valid; gibbs() or ecm(prior = \"map\") regularize)",
                                            paste(bd, collapse = ", ")) else ""), call. = FALSE)
  }
  es <- em_estep(p, E)
  ll <- sum(es$ll); npar <- length(x); N <- E$N
  V <- NULL
  if (se) {
    H <- em_hessian(x, E)
    V <- tryCatch(solve(-H), error = function(e) NULL)
    if (is.null(V) || any(!is.finite(diag(V))) || any(diag(V) < 0)) {
      warning("the observed information is not positive definite: no standard errors", call. = FALSE); V <- NULL
    }
  }
  mo <- em_moments(es)
  person <- list(theta = list(mean = mo$Et, sd = sqrt(pmax(mo$Et2 - mo$Et^2, 0))))
  if (E$rt) person$zeta <- list(mean = mo$Ez, sd = sqrt(pmax(mo$Ez2 - mo$Ez^2, 0)))
  tab <- em_table(x, V, E); vc <- attr(tab, "vcov"); attr(tab, "vcov") <- NULL
  model$ecm <- list(estimates = tab, par = p, x = x, logLik = ll, npar = npar, nobs = N,
                    AIC = if (map) NA_real_ else -2 * ll + 2 * npar, BIC = if (map) NA_real_ else -2 * ll + log(N) * npar, iterations = it, converged = conv,
                    trace = trace, vcov = vc, person = person, map = map, priors = if (map) priors, nodes = Q,
                    settings = model$settings, X = E$X, secs = as.numeric(Sys.time() - t0, units = "secs"))
  model
}

# model context: data, design, quadrature
em_setup <- function(object, Q, map, priors = rtirt_priors()) {
  type <- class(object)[1]; info <- .model_info[[type]]; S <- object$settings; D <- object$data
  if (isTRUE(S$intercept)) stop("ecm() needs intercept = FALSE: without priors the intercepts of the latent regressions are not identified (the item parameters carry the location)", call. = FALSE)
  N <- nrow(D$Y); K <- ncol(D$Y)
  X <- if (info$X) D$X else matrix(0, N, 0)
  nd <- calc_nodes(Q)
  L <- D$log_t
  list(Y = D$Y, kappa = D$Y - 0.5, L = L, X = X, XX = crossprod(X), N = N, K = K, P = ncol(X), eng = info$engine, rt = info$rt,
       fixed = S$speed_var == "fixed" && info$engine %in% c("rtirt", "cross"), one_pl = S$itemtype == "1pl", map = map,
       z = nd$z, lw = nd$logw, items = D$items, xn = colnames(X) %||% character(),
       mu_l = if (info$rt) lambda_prior(priors, L)[1] else NA, sd_l = if (info$rt) lambda_prior(priors, L)[2] else NA, pr = priors)
}

em_start <- function(E) {
  pb <- pmin(pmax(colMeans(E$Y), 0.02), 0.98)
  p <- list(a = rep(1, E$K), b = -1.3 * stats::qlogis(pb))
  if (E$rt) {
    p$lambda <- colMeans(E$L); s2 <- var_floor(apply(E$L, 2, stats::var))
    vs <- max(stats::var(rowMeans(E$L)) - mean(s2) / E$K, 0.05)
    p$sigma2t <- pmax(s2 - vs, 0.05 * s2)
  }
  switch(E$eng,
    mlirt = p$beta <- rep(0, E$P),
    rtirt = { p$beta1 <- p$beta2 <- rep(0, E$P); p$c <- 0; p$v <- if (E$fixed) 1 else vs },
    latent = { p$beta <- rep(0, E$P); p$gamma <- 0; p$s <- vs },
    cross = { p$rho <- rep(0, E$K); p$s <- if (E$fixed) 1 else vs })
  p
}

# natural-scale list <-> unconstrained vector (log for a and variances, atanh for a fixed-variance correlation)
em_pack <- function(p, E)
  c(if (!E$one_pl) log(p$a), p$b, if (E$rt) c(p$lambda, log(p$sigma2t)),
    switch(E$eng, mlirt = p$beta,
           rtirt = c(p$beta1, p$beta2, if (E$fixed) atanh(p$c) else c(p$c, log(p$v))),
           latent = c(p$beta, p$gamma, log(p$s)),
           cross = c(p$rho, if (!E$fixed) log(p$s))))

em_unpack <- function(x, E) {
  K <- E$K; P <- E$P; i <- 0
  take <- function(n) { v <- x[i + seq_len(n)]; i <<- i + n; v }
  p <- list(a = if (E$one_pl) rep(1, K) else exp(take(K)), b = take(K))
  if (E$rt) { p$lambda <- take(K); p$sigma2t <- exp(take(K)) }
  switch(E$eng,
    mlirt = { p$beta <- take(P) },
    rtirt = { p$beta1 <- take(P); p$beta2 <- take(P)
              if (E$fixed) { p$c <- tanh(take(1)); p$v <- 1 - p$c^2 } else { p$c <- take(1); p$v <- exp(take(1)) } },
    latent = { p$beta <- take(P); p$gamma <- take(1); p$s <- exp(take(1)) },
    cross = { p$rho <- take(K); p$s <- if (E$fixed) 1 else exp(take(1)) })
  p
}

em_xb <- function(X, beta) if (length(beta)) drop(X %*% beta) else rep(0, nrow(X))
softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# E-step: nodes TH (N x Q), posterior weights W, posterior mean MZ (N x Q) and precision Pz of speed
em_estep <- function(p, E) {
  N <- E$N; K <- E$K; Q <- length(E$z)
  mu_t <- switch(E$eng, mlirt = em_xb(E$X, p$beta), rtirt = em_xb(E$X, p$beta1), rep(0, N))
  TH <- outer(mu_t, E$z, "+")
  common <- all(mu_t == 0)                           # common nodes: Bock-Aitkin products; otherwise C++ (person-specific nodes)
  L <- calc_loglik_nodes(E$Y, p$a, p$b, if (common) E$z else TH) + matrix(E$lw, N, Q, byrow = TRUE)
  MZ <- NULL; Pz <- NA
  if (E$rt) {
    C <- E$L - matrix(p$lambda, N, K, byrow = TRUE); w <- 1 / p$sigma2t
    rr <- if (E$eng == "cross") p$rho else rep(0, K)
    S1 <- drop(C %*% w); S2 <- drop(C^2 %*% w); Cr <- drop(C %*% (w * rr)); wr <- sum(w * rr); wr2 <- sum(w * rr^2)
    M <- switch(E$eng, rtirt = em_xb(E$X, p$beta2) + p$c * (TH - mu_t), latent = em_xb(E$X, p$beta) + p$gamma * TH, 0 * TH)
    # speed u ~ N(M, vz) seen through o_i = -(c_i + rho_i theta) = u + e_i: sums of w o and w o^2 at every node
    sp <- normal_given_sums(-(S1 + TH * wr), S2 + 2 * TH * Cr + TH^2 * wr2, w, M, if (E$eng == "rtirt") p$v else p$s)
    MZ <- sp$mean; Pz <- 1 / sp$var; L <- L + sp$logI
  }
  nw <- calc_node_weights(L)
  list(ll = nw$ll, TH = TH, W = nw$W, MZ = MZ, Pz = Pz, mu_t = mu_t, common = common)
}

em_moments <- function(es) {
  W <- es$W; TH <- es$TH
  m <- list(Et = rowSums(W * TH), Et2 = rowSums(W * TH^2))
  if (!is.null(es$MZ)) { MZ <- es$MZ
    m$Ez <- rowSums(W * MZ); m$Ez2 <- rowSums(W * MZ^2) + 1 / es$Pz; m$Etz <- rowSums(W * TH * MZ) }
  m
}

# item sums over persons and nodes: Polya-Gamma weights (MM step) or residuals y - p (score)
em_irt_sums <- function(p, E, es, score = FALSE)
  calc_item_sums(E$Y, p$a, p$b, if (isTRUE(es$common)) E$z else es$TH, es$W, score)

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
  N <- E$N; K <- E$K; P <- E$P; X <- E$X; mo <- em_moments(es); new <- p
  # items: PG minorizer, eta = a theta - d
  ac <- em_irt_sums(p, E, es); S0 <- ac$S0; S1 <- ac$S1; S2 <- ac$S2
  T0 <- colSums(E$kappa); T1 <- colSums(E$kappa * mo$Et)
  if (!E$map && E$one_pl) new$b <- (S1 - T0) / S0
  else {                                             # PG minorizer step (ML), or conditional modes with the Gibbs priors (MAP)
    up <- update_items_logistic(ac, T0, T1, p$a, p$b, prior_a = if (E$map) unname(E$pr$a), prior_b = if (E$map) unname(E$pr$b))
    new$b <- up$b; if (!E$one_pl) new$a <- up$a
  }
  # response times: log T + zeta = lambda - rho theta + e
  if (E$rt) {
    Lt <- E$L
    if (E$eng == "cross") {
      w0 <- if (E$map) 1 / p$sigma2t else rep(1, K)
      n11 <- N * w0 + E$map / E$sd_l^2; n12 <- -sum(mo$Et) * w0; n22 <- sum(mo$Et2) * w0 + E$map / E$pr$rho[["sd"]]^2
      r1 <- (colSums(Lt) + sum(mo$Ez)) * w0 + E$map * E$mu_l / E$sd_l^2
      r2 <- -(drop(crossprod(Lt, mo$Et)) + sum(mo$Etz)) * w0
      dt <- n11 * n22 - n12^2
      new$lambda <- (n22 * r1 - n12 * r2) / dt; new$rho <- (n11 * r2 - n12 * r1) / dt
    } else {
      new$lambda <- if (!E$map) (colSums(Lt) + sum(mo$Ez)) / N
                    else ((colSums(Lt) + sum(mo$Ez)) / p$sigma2t + E$mu_l / E$sd_l^2) / (N / p$sigma2t + 1 / E$sd_l^2)
    }
    SS <- em_rt_ss(Lt, new$lambda, if (E$eng == "cross") new$rho else rep(0, K), mo)
    new$sigma2t <- if (!E$map) SS / N
      else if (E$pr$sigma2t_family == "inverse-gamma") (E$pr$sigma2t[["scale"]] + SS / 2) / (E$pr$sigma2t[["shape"]] + N / 2 + 1)
      else vapply(seq_len(K), function(k) exp(stats::optimize(function(u) -N / 2 * u - SS[k] / 2 * exp(-u) +
                    var_lprior(exp(u), E$pr$sigma2t, "half-t")$lp, c(-15, 8), maximum = TRUE)$maximum), 0)
  }
  # structural part
  if (E$eng == "mlirt" && P) new$beta <- drop(solve(E$XX, crossprod(X, mo$Et)))
  if (E$eng %in% c("rtirt", "latent")) {               # regression of speed on (x, theta)
    xt <- drop(crossprod(X, mo$Et)); xz <- drop(crossprod(X, mo$Ez))
    G <- rbind(cbind(E$XX, xt), c(xt, sum(mo$Et2))); h <- c(xz, sum(mo$Etz))
    if (E$eng == "rtirt" && E$fixed) {                 # Var(zeta | x) = 1: c = correlation, v = 1 - c^2
      ss <- function(c) { r <- xz - c * xt; sum(mo$Ez2) - 2 * c * sum(mo$Etz) + c^2 * sum(mo$Et2) - if (P) sum(r * solve(E$XX, r)) else 0 }
      f <- function(c) -N / 2 * log(1 - c^2) - ss(c) / (2 * (1 - c^2))
      c <- stats::optimize(f, c(-0.9999, 0.9999), maximum = TRUE, tol = 1e-10)$maximum
      bb <- c(if (P) solve(E$XX, xz - c * xt), c); v <- 1 - c^2
    } else {
      bb <- drop(solve(G, h)); v <- max((sum(mo$Ez2) - 2 * sum(bb * h) + sum(bb * (G %*% bb))) / N, 1e-8)
    }
    if (E$eng == "rtirt") {
      new$beta1 <- if (P) drop(solve(E$XX, xt)) else numeric(0)
      new$c <- bb[P + 1]; new$v <- v; new$beta2 <- bb[seq_len(P)] + new$c * new$beta1
    } else { new$beta <- bb[seq_len(P)]; new$gamma <- bb[P + 1]; new$s <- v }
  }
  if (E$eng == "cross" && !E$fixed) new$s <- max(sum(mo$Ez2) / N, 1e-8)
  new
}

# per item: sum_i E[(log T + zeta - lambda + rho theta)^2]
em_rt_ss <- function(Lt, lambda, rho, mo) {
  C <- sweep(Lt, 2, lambda)
  colSums(C^2) + 2 * colSums(C * mo$Ez) + 2 * rho * colSums(C * mo$Et) + sum(mo$Ez2) + 2 * rho * sum(mo$Etz) + rho^2 * sum(mo$Et2)
}

# score of the marginal log-likelihood (+ log prior) in the unconstrained parameters (Fisher identity)
em_score <- function(x, E) {
  p <- em_unpack(x, E); es <- em_estep(p, E); mo <- em_moments(es); N <- E$N; K <- E$K; X <- E$X
  ir <- em_irt_sums(p, E, es, score = TRUE)
  s_a <- ir$S1 - p$b * ir$S0; s_b <- -p$a * ir$S0
  if (E$map) {
    pa <- E$pr$a
    s_a <- s_a - (p$a - pa[1]) / pa[2]^2
    s_b <- s_b - (p$b - E$pr$b[["mean"]]) / E$pr$b[["sd"]]^2
  }
  out <- c(if (!E$one_pl) p$a * s_a, s_b)
  if (E$rt) {
    rho <- if (E$eng == "cross") p$rho else rep(0, K)
    C <- sweep(E$L, 2, p$lambda); w <- 1 / p$sigma2t
    s_l <- (colSums(C) + sum(mo$Ez) + rho * sum(mo$Et)) * w
    SS <- em_rt_ss(E$L, p$lambda, rho, mo)
    s_s2 <- -N / 2 * w + SS / 2 * w^2
    if (E$map) { s_l <- s_l - (p$lambda - E$mu_l) / E$sd_l^2; s_s2 <- s_s2 + var_lprior(p$sigma2t, E$pr$sigma2t, E$pr$sigma2t_family)$d }
    out <- c(out, s_l, p$sigma2t * s_s2)
  }
  str <- switch(E$eng,
    mlirt = drop(crossprod(X, mo$Et - em_xb(X, p$beta))),
    rtirt = {
      c <- p$c; v <- p$v; m1 <- em_xb(X, p$beta1); g <- em_xb(X, p$beta2) - c * m1
      Eu <- mo$Ez - c * mo$Et - g
      Eu2 <- mo$Ez2 - 2 * c * mo$Etz + c^2 * mo$Et2 - 2 * g * (mo$Ez - c * mo$Et) + g^2
      Eut <- mo$Etz - c * mo$Et2 - g * mo$Et
      s_b1 <- drop(crossprod(X, mo$Et - m1)) - c / v * drop(crossprod(X, Eu))
      s_b2 <- drop(crossprod(X, Eu)) / v
      s_c <- sum(Eut - m1 * Eu) / v; s_v <- sum(-1 / (2 * v) + Eu2 / (2 * v^2))
      c(s_b1, s_b2, if (E$fixed) (1 - c^2) * (s_c - 2 * c * s_v) else c(s_c, v * s_v))
    },
    latent = {
      gm <- p$gamma; s <- p$s; g <- em_xb(X, p$beta)
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
    })
  c(out, str)
}

# observed information: central differences of the analytic score
em_hessian <- function(x, E) calc_hessian(function(z) em_score(z, E), x)

# reported parameters (the names of gibbs()) as a function of the unconstrained vector
em_report <- function(x, E) {
  p <- em_unpack(x, E)
  v <- c(p$a, p$b, if (E$rt) c(p$lambda, p$sigma2t),
         switch(E$eng, mlirt = p$beta,
                rtirt = c(p$beta1, p$beta2, p$c / sqrt(p$c^2 + p$v), p$c^2 + p$v),
                latent = c(p$beta, p$gamma, p$s), cross = c(p$rho, p$s)))
  stats::setNames(v, par_names(E$eng, E$items, E$xn))
}

em_table <- function(x, V, E) {
  est <- em_report(x, E); n <- length(est)
  fixed <- rep(FALSE, n); nm <- names(est)
  if (E$one_pl) fixed[startsWith(nm, "a[")] <- TRUE
  if (E$fixed) fixed[nm == "var_speed"] <- TRUE
  se <- rep(NA_real_, n); Vr <- NULL
  if (!is.null(V)) {
    J <- vapply(seq_along(x), function(j) {
      h <- 1e-6 * max(1, abs(x[j])); e <- replace(numeric(length(x)), j, h)
      (em_report(x + e, E) - em_report(x - e, E)) / (2 * h)
    }, numeric(n))
    J <- matrix(J, n); Vr <- J %*% V %*% t(J); dimnames(Vr) <- list(nm, nm)
    se <- sqrt(pmax(diag(Vr), 0)); se[fixed] <- NA
  }
  pos <- grepl("^(a|sigma2t)\\[", nm) | nm == "var_speed"; cor <- nm == "cor_ability_speed"
  lo <- est - 1.959964 * se; hi <- est + 1.959964 * se
  lo[pos] <- est[pos] * exp(-1.959964 * se[pos] / est[pos]); hi[pos] <- est[pos] * exp(1.959964 * se[pos] / est[pos])
  if (any(cor)) { z <- atanh(est[cor]); sz <- se[cor] / (1 - est[cor]^2); lo[cor] <- tanh(z - 1.959964 * sz); hi[cor] <- tanh(z + 1.959964 * sz) }
  zv <- est / se
  tab <- data.frame(parameter = nm, est = unname(est), se = unname(se), z = unname(zv), q025 = unname(lo), q975 = unname(hi),
                    sig = ifelse(!is.na(lo) & (lo > 0 | hi < 0), "*", ""), row.names = NULL)
  attr(tab, "vcov") <- Vr
  tab
}

# starting values of a Gibbs chain from the ECM fit (jittered)
em_inits <- function(model, engine, X, fixed) {
  f <- model$ecm; p <- f$par; N <- nrow(model$data$Y); K <- ncol(model$data$Y)
  jit <- function(v, s = 0.05) v + stats::rnorm(length(v), 0, s)
  pad <- function(b) if (ncol(X) > length(b)) c(rep(0, ncol(X) - length(b)), b) else b      # intercept column first
  init <- list(theta = jit(f$person$theta$mean, 0.1), a = p$a * exp(stats::rnorm(K, 0, 0.05)), b = jit(p$b))
  if (engine != "mlirt") init <- c(init, list(zeta = jit(f$person$zeta$mean, 0.05), lambda = jit(p$lambda, 0.02),
                                               sigma2 = p$sigma2t * exp(stats::rnorm(K, 0, 0.05))))
  switch(engine,
    mlirt = c(init, list(beta = pad(jit(p$beta)))),
    rtirt = c(init, list(beta1 = pad(jit(p$beta1)), beta2 = pad(jit(p$beta2)),
                         c = if (fixed) tanh(atanh(p$c) + stats::rnorm(1, 0, 0.05)) else jit(p$c), v = if (fixed) 1 else p$v)),
    latent = c(init, list(beta = c(pad(jit(p$beta)), jit(p$gamma)), s = p$s * exp(stats::rnorm(1, 0, 0.05)))),
    cross = c(init, list(rho = jit(p$rho), s = if (fixed) 1 else p$s * exp(stats::rnorm(1, 0, 0.05)))))
}

# check_ecm() for a built-in model, with the E- and M-step of ecm(): one M-step never lowers the objective,
# Fisher's identity (em_score() = numerical gradient of the objective), and the fit is a maximum
em_check_builtin <- function(model, prior, nodes, n_start, add) {
  if (inherits(model, "rtirt_qset") || isTRUE(model$qr)) stop("check_ecm(): quantile models have no EM; check their sampler with check_sampler()", call. = FALSE)
  Q <- nodes %||% 41; E <- em_setup(model, Q, prior == "map")
  obj <- function(x) { p <- em_unpack(x, E); em_obj(em_estep(p, E), p, E) }
  grad <- function(f, x) vapply(seq_along(x), function(k) { h <- 1e-5 * max(1, abs(x[k])); (f(replace(x, k, x[k] + h)) - f(replace(x, k, x[k] - h))) / (2 * h) }, 0)
  xh <- suppressWarnings(ecm(model, nodes = Q, prior = prior, se = FALSE))$ecm$x
  x0 <- em_pack(em_start(E), E); drop_l <- 0; fisher <- 0
  for (r in seq_len(n_start)) {
    x <- x0 + stats::rnorm(length(x0), 0, 0.3)                      # jittered on the unconstrained scale
    p <- em_unpack(x, E); es <- em_estep(p, E); o0 <- em_obj(es, p, E)
    drop_l <- max(drop_l, (o0 - obj(em_pack(em_mstep(p, E, es), E))) / max(1, abs(o0)))
    g <- grad(obj, x); fisher <- max(fisher, max(abs(g - em_score(x, E))) / max(1, max(abs(g))))
  }
  add("one iteration never lowers the objective (relative drop)", drop_l, 1e-8, "em_mstep(): all blocks")
  add("Fisher identity: gradient of the objective = gradient of Q (relative difference)", fisher, 1e-4, "em_score()")
  x <- xh; o <- stats::optim(x, obj, method = "BFGS", control = list(fnscale = -1, maxit = 200))
  add("ECM stops at a maximum: gain of a general optimizer (relative)", (o$value - obj(x)) / max(1, abs(obj(x))), 1e-6, "")
  add("ECM stops at a maximum: largest |gradient|", max(abs(grad(obj, x))), 1e-2, "")
}

# ---- methods for ECM fits --------------------------------------------------------------------

has_em <- function(x) is.null(x$post) && !is.null(x$ecm)
need_em <- function(x) if (is.null(x$ecm)) stop("run ecm() first", call. = FALSE)

print_em <- function(x) {
  f <- x$ecm
  cat(sprintf("  %s (ECM, %d Gauss-Hermite nodes): %d iterations%s, %.1f s\n",
              if (f$map) "marginal posterior mode (item priors of gibbs())" else "maximum marginal likelihood", f$nodes, f$iterations,
              if (f$converged) "" else " (not converged)", f$secs))
  if (f$map) cat(sprintf("  log posterior %.3f, npar %d (no AIC / BIC: the objective includes the priors)\n", f$logLik, f$npar))
  else cat(sprintf("  logLik %.3f, npar %d, AIC %.1f, BIC %.1f\n", f$logLik, f$npar, f$AIC, f$BIC))
  rel <- reliability(x)
  cat(sprintf("  Reliability: %s\n", paste(sprintf("%s %.3f", names(rel), rel), collapse = ", ")))
  cat("Use summary(), coef(), scores(), estimates(), fit_indices(); gibbs(fit, init = \"ecm\") to sample from here.\n")
}

#' @export
logLik.rtirt <- function(object, ...) {
  need_em(object)
  structure(object$ecm$logLik, df = object$ecm$npar, nobs = object$ecm$nobs, class = "logLik")
}

#' @export
vcov.rtirt <- function(object, ...) {
  need_em(object)
  if (is.null(object$ecm$vcov)) stop("no standard errors: ecm(..., se = TRUE)", call. = FALSE)
  object$ecm$vcov
}

summary_em <- function(object, digits) {
  e <- object$ecm$estimates
  if (isTRUE(object$ecm$map)) {                      # item priors only; structural parameters by ML
    pr <- prior_strings(object); item <- grepl("^(a|b|lambda|sigma2t)\\[", e$parameter)
    e$prior <- ifelse(item, unname(pr[e$parameter]), "(ML)")
  }
  p <- rtirt_partable(e, FALSE, FALSE, sprintf("%s, SE, Wald z and 95%% interval (log scale for a, sigma2t, var_speed; Fisher-z for the correlation)",
                                               if (isTRUE(object$ecm$map)) "posterior mode under the item priors" else "ML estimate"))
  w <- rtirt_item_wide(p, object$data$items)
  structure(list(header = rtirt_header(object),
                 fit = rtable(fit_indices(object), "marginal likelihood (persons integrated out)"),
                 items = rtable(w$est, "item parameters (ECM estimates)"),
                 se = rtable(w$se, "standard errors of the item parameters"),
                 parameters = p,
                 reliability = rtable(as.data.frame(as.list(reliability(object))), "reliability of the EAP scores (final E-step)")),
            class = c("summary.rtirt", "rtirt_summary"))
}
