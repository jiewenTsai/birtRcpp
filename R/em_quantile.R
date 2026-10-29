# ECM for the quantile (ALD) models of "latent" and "cross" ------------------------------
#
# The ALD is a working likelihood (Yu & Moyeed, 2001): the q-quantile is linear when its
# check loss is minimized. Ability is integrated on the Gauss-Hermite nodes and speed in closed
# form given ability (the integrals of marginal_loglik()), so the marginal working
# log-likelihood is smooth in every parameter.
# E-step: the posterior weights of the nodes.
# CM-step 1: the 2PL items, the Polya-Gamma step of the normal models (em_mstep_items()).
# CM-step 2: the RT and structural parameters maximize sum_jk W_jk L_rt(j, k), the weights
#   fixed, by nlminb() with the analytic gradient (posterior expectations of speed: Fisher
#   identity). Each step raises the objective, so the log-likelihood never decreases.
# Standard errors: sandwich H^-1 (sum_j s_j s_j') H^-1, from the score s_j of each person and
#   H the numerical Jacobian of the total score (the information identity fails under a
#   working likelihood). Derivation: designs/quantile-ecm.md of the lab folder.
# Internal parameters as in em.R: for "cross" the slot sigma2t holds the ALD scales of the log
# times; for "latent" the slot s holds the ALD scale of speed (the names of gibbs()).

em_qr_fit <- function(model, max_iter, tol, Q, se, accelerate, verbose, t0) {
  E <- em_setup(model, Q, map = FALSE)
  E$q <- model$q
  x <- em_pack(em_qr_start(model, E, Q, max_iter, tol), E)
  ab <- em_qr_ab_index(E)
  F <- function(x) {
    p <- em_unpack(x, E)
    es <- em_qr_estep(p, E)
    new <- em_mstep_items(p, E, es, em_moments(es))
    x_new <- em_pack(new, E)
    x_new[-ab] <- em_qr_cm_rt(x_new, es$W, E)
    list(obj = sum(es$ll), x = x_new)
  }
  progress <- if (verbose) function(it, obj) if (it %% 10 == 0) message(sprintf("ECM %d: working logLik %.4f", it, obj))
  sq <- em_squarem(x, F, max_iter = max_iter, tol = tol, accelerate = accelerate, progress = progress)
  x <- sq$x
  if (!sq$converged) {
    warning(sprintf("ecm() did not converge in %d iterations (raise max_iter)", max_iter), call. = FALSE)
  }
  p <- em_unpack(x, E)
  es <- em_qr_estep(p, E, moments = TRUE)
  V <- NULL
  if (se) {
    V <- em_qr_sandwich(x, E)
  }
  tab <- em_table(x, V, E)
  vc <- attr(tab, "vcov")
  attr(tab, "vcov") <- NULL
  ll <- sum(es$ll)
  npar <- length(x)
  N <- E$N
  person <- list(
    theta = list(mean = rowSums(es$W * es$TH), sd = sqrt(pmax(rowSums(es$W * es$TH^2) - rowSums(es$W * es$TH)^2, 0))),
    zeta = list(mean = es$Ez, sd = sqrt(pmax(es$Ez2 - es$Ez^2, 0)))
  )
  model$ecm <- list(
    estimates = tab, par = p, x = x, logLik = ll, npar = npar, nobs = N,
    AIC = -2 * ll + 2 * npar, BIC = -2 * ll + log(N) * npar, iterations = sq$iterations, converged = sq$converged,
    trace = sq$trace, vcov = vc, person = person, map = FALSE, priors = NULL, nodes = Q, se_type = "sandwich",
    settings = model$settings, X = E$X, secs = as.numeric(Sys.time() - t0, units = "secs")
  )
  model
}

# positions of log a and b in the unconstrained vector (the rest: RT and structural)
em_qr_ab_index <- function(E) seq_len(if (E$one_pl) E$K else 2 * E$K)

# start: ECM of the mean model, its RT location moved to the q-quantile and its scale matched
# to the variance of the ALD, Var = scale^2 (1 - 2q + 2q^2) / (q (1 - q))^2
em_qr_start <- function(model, E, Q, max_iter, tol) {
  mean_model <- model
  mean_model$qr <- FALSE
  mean_model$q <- NA_real_
  fit <- suppressWarnings(ecm.rtirt(mean_model, max_iter = max_iter, tol = tol, nodes = Q, se = FALSE))
  p <- fit$ecm$par
  q <- E$q
  to_scale <- q * (1 - q) / sqrt(1 - 2 * q + 2 * q^2)
  if (E$eng == "cross") {
    sd_t <- sqrt(p$sigma2t)
    p$lambda <- p$lambda + sd_t * stats::qnorm(q)
    p$sigma2t <- sd_t * to_scale
  } else {
    sd_u <- sqrt(p$s)
    p$lambda <- p$lambda - sd_u * stats::qnorm(q)
    p$s <- sd_u * to_scale
  }
  p
}

# E-step: weights of the nodes from the 2PL and the RT part; with moments = TRUE also the
# posterior mean and second moment of speed
em_qr_estep <- function(p, E, moments = FALSE) {
  N <- E$N
  Q <- length(E$z)
  TH <- matrix(E$z, N, Q, byrow = TRUE)
  L <- calc_loglik_nodes(E$Y, p$a, p$b, E$z) + matrix(E$lw, N, Q, byrow = TRUE)
  rt <- em_qr_rt(p, E, TH, moments)
  nw <- calc_node_weights(L + rt$L)
  es <- list(ll = nw$ll, TH = TH, W = nw$W, common = TRUE)
  if (moments) {
    es$Ez <- rowSums(nw$W * rt$Ez)
    es$Ez2 <- rowSums(nw$W * rt$Ez2)
  }
  es
}

# CM-step 2: the RT and structural parameters (unconstrained) given the node weights W; the
# value and the gradient come from one pass (cached for nlminb, which asks for both at a point)
em_qr_cm_rt <- function(x, W, E, iter_max = 2) {
  ab <- em_qr_ab_index(E)
  last <- list(y = NULL)
  at <- function(y) {
    if (!identical(y, last$y)) {
      p <- em_unpack(replace(x, -ab, y), E)
      parts <- em_qr_rt_score(p, E, W, value = TRUE)
      last <<- list(y = y, value = -parts$value, gradient = -colSums(parts$score))
    }
    last
  }
  start <- x[-ab]
  fit <- stats::nlminb(start, function(y) at(y)$value, function(y) at(y)$gradient,
                       control = list(iter.max = iter_max, rel.tol = 1e-12))
  if (is.finite(fit$objective) && fit$objective <= at(start)$value) fit$par else start
}

# RT part at the nodes: L (N x Q), and with moments = TRUE the posterior mean Ez and second
# moment Ez2 of speed at each node
em_qr_rt <- function(p, E, TH, moments = FALSE) {
  if (E$eng == "cross") {
    B <- em_qr_cross_rows(p, E, TH)
    v <- p$s
    if (!moments) {
      return(list(L = matrix(.log_ald_rt(B, p$sigma2t, v, E$q), E$N)))
    }
    g <- .ald_rt_grad(B, p$sigma2t, v, E$q)
    return(list(L = matrix(g$ll, E$N), Ez = matrix(g$ez, E$N), Ez2 = matrix(g$ez2, E$N)))
  }
  lat <- em_qr_latent(p, E, TH)
  out <- list(L = lat$const + lat$F)
  if (moments) {
    out$Ez <- lat$Ez
    out$Ez2 <- lat$Ez2
  }
  out
}

# "cross": the rows b_ijk = log t_ij - lambda_i + rho_i theta_k, persons within nodes
em_qr_cross_rows <- function(p, E, TH) {
  C <- E$L - matrix(p$lambda, E$N, E$K, byrow = TRUE)
  rows <- lapply(seq_len(ncol(TH)), function(k) C + outer(TH[, k], p$rho))
  do.call(rbind, rows)
}

# "latent": the normal RT likelihood is a normal kernel in speed, N(zeta; zh_j, 1 / A), with
# constant const_j; F = log of its convolution with the ALD of speed (mean m_jk, scale s), and
# the derivatives of F (designs/quantile-ecm.md, section 4)
em_qr_latent <- function(p, E, TH) {
  N <- E$N
  Q <- ncol(TH)
  K <- E$K
  q <- E$q
  s <- p$s
  w <- 1 / p$sigma2t
  C <- E$L - matrix(p$lambda, N, K, byrow = TRUE)
  A <- sum(w)
  S1 <- drop(C %*% w)
  S2 <- drop(C^2 %*% w)
  zh <- -S1 / A
  tau2 <- 1 / A
  tau <- sqrt(tau2)
  const <- -0.5 * K * log(2 * pi) + 0.5 * sum(log(w)) - 0.5 * (S2 - S1^2 / A) + 0.5 * log(2 * pi / A)
  m <- em_xb(E$X, p$beta) + p$gamma * TH
  d <- zh - m
  a <- q / s
  b <- (1 - q) / s
  u_pos <- (d - a * tau2) / tau
  u_neg <- (-d - b * tau2) / tau
  l_pos <- -a * d + a^2 * tau2 / 2 + stats::pnorm(u_pos, log.p = TRUE)
  l_neg <- b * d + b^2 * tau2 / 2 + stats::pnorm(u_neg, log.p = TRUE)
  top <- pmax(l_pos, l_neg)
  lse <- top + log(exp(l_pos - top) + exp(l_neg - top))
  pi_pos <- exp(l_pos - lse)
  pi_neg <- exp(l_neg - lse)
  mills <- function(u) exp(stats::dnorm(u, log = TRUE) - stats::pnorm(u, log.p = TRUE))
  h_pos <- mills(u_pos)
  h_neg <- mills(u_neg)
  F_d <- pi_pos * (-a + h_pos / tau) + pi_neg * (b - h_neg / tau)
  F_s <- -1 / s - a / s * pi_pos * (-d + a * tau2 - tau * h_pos) - b / s * pi_neg * (d + b * tau2 - tau * h_neg)
  F_tau2 <- pi_pos * (a^2 / 2 - h_pos * (d / tau2 + a) / (2 * tau)) + pi_neg * (b^2 / 2 + h_neg * (d / tau2 - b) / (2 * tau))
  # speed a posteriori: TN(zh - a tau2, tau2; (m, Inf)) w.p. pi_pos, TN(zh + b tau2, tau2; (-Inf, m)) w.p. pi_neg
  mean_pos <- zh - a * tau2 + tau * h_pos     # upper truncation at m: (m - mean) / tau = -u_pos
  mean_neg <- zh + b * tau2 - tau * h_neg
  var_pos <- tau2 * (1 - h_pos * (h_pos + u_pos))
  var_neg <- tau2 * (1 - h_neg * (h_neg + u_neg))
  Ez <- pi_pos * mean_pos + pi_neg * mean_neg
  Ez2 <- pi_pos * (var_pos + mean_pos^2) + pi_neg * (var_neg + mean_neg^2)
  list(
    const = const, F = log(q * (1 - q) / s) + lse, F_d = F_d, F_s = F_s, F_tau2 = F_tau2,
    zh = zh, A = A, C = C, w = w, Ez = Ez, Ez2 = Ez2
  )
}

# per-person score of the RT and structural parameters (unconstrained scale), N x n_rt, with
# the node weights W: sum_k W_jk dL_rt(j, k) / dx; with value = TRUE a list of the score and
# sum_jk W_jk L_rt(j, k)
em_qr_rt_score <- function(p, E, W, value = FALSE) {
  with_value <- function(score, L) if (value) list(score = score, value = sum(W * L)) else score
  N <- E$N
  K <- E$K
  Q <- length(E$z)
  TH <- matrix(E$z, N, Q, byrow = TRUE)
  if (E$eng == "cross") {
    g <- .ald_rt_grad(em_qr_cross_rows(p, E, TH), p$sigma2t, p$s, E$q)
    wv <- as.vector(W)
    person <- rep(seq_len(N), times = Q)
    by_person <- function(M) rowsum(M, person, reorder = FALSE)
    s_lambda <- -by_person(wv * g$gb)
    s_rho <- by_person(wv * as.vector(TH) * g$gb)
    s_scale <- by_person(wv * g$gs) * matrix(p$sigma2t, N, K, byrow = TRUE)
    s_v <- if (E$fixed) NULL else by_person(wv * g$gv) * p$s
    out <- cbind(s_lambda, s_scale, s_rho, s_v)
    dimnames(out) <- NULL
    return(with_value(out, g$ll))
  }
  lat <- em_qr_latent(p, E, TH)
  fz <- rowSums(W * lat$F_d)
  ft <- rowSums(W * lat$F_tau2)
  fm <- -W * lat$F_d
  w <- lat$w
  A <- lat$A
  e <- lat$C + lat$zh
  s_lambda <- sweep(e + fz / A, 2, w, "*")
  s_w <- matrix(0.5 / w - 0.5 / A, N, K, byrow = TRUE) - e^2 / 2 - e * fz / A - ft / A^2
  s_log_sigma2 <- -sweep(s_w, 2, w, "*")
  s_beta <- if (E$P) E$X * rowSums(fm) else NULL
  s_gamma <- rowSums(fm * TH)
  s_log_s <- rowSums(W * lat$F_s) * p$s
  with_value(cbind(s_lambda, s_log_sigma2, s_beta, s_gamma, s_log_s, deparse.level = 0), lat$const + lat$F)
}

# per-person score of the 2PL items (log a, b), N x n_ab
em_qr_item_score <- function(p, E, W, TH) {
  cols <- lapply(seq_len(E$K), function(i) {
    P <- stats::plogis(p$a[i] * (TH - p$b[i]))
    R <- (E$Y[, i] - P) * W
    list(log_a = p$a[i] * rowSums(R * (TH - p$b[i])), b = -p$a[i] * rowSums(R))
  })
  s_b <- vapply(cols, `[[`, numeric(E$N), "b")
  if (E$one_pl) {
    return(s_b)
  }
  cbind(vapply(cols, `[[`, numeric(E$N), "log_a"), s_b)
}

# score of every person (unconstrained scale), N x length(x)
em_qr_person_score <- function(x, E) {
  p <- em_unpack(x, E)
  es <- em_qr_estep(p, E)
  matrix(cbind(em_qr_item_score(p, E, es$W, es$TH), em_qr_rt_score(p, E, es$W)), E$N)
}

em_qr_score <- function(x, E) colSums(em_qr_person_score(x, E))

# sandwich covariance of the unconstrained parameters, or NULL with a warning
em_qr_sandwich <- function(x, E) {
  H <- calc_hessian(function(z) em_qr_score(z, E), x)
  S <- em_qr_person_score(x, E)
  bread <- tryCatch(solve(-H), error = function(e) NULL)
  V <- if (is.null(bread)) NULL else bread %*% crossprod(S) %*% bread
  if (is.null(V) || any(!is.finite(diag(V))) || any(diag(V) < 0)) {
    warning("the observed information is not positive definite: no standard errors", call. = FALSE)
    return(NULL)
  }
  V
}

# check_ecm() for a quantile model: one iteration never lowers the working log-likelihood, the
# analytic score is its gradient, and the fit is a maximum
em_qr_check <- function(model, nodes, n_start, add) {
  Q <- nodes %||% 41
  E <- em_setup(model, Q, map = FALSE)
  E$q <- model$q
  ab <- em_qr_ab_index(E)
  loglik <- function(x) sum(em_qr_estep(em_unpack(x, E), E)$ll)
  gradient <- function(x) {
    vapply(seq_along(x), function(k) {
      h <- 1e-5 * max(1, abs(x[k]))
      (loglik(replace(x, k, x[k] + h)) - loglik(replace(x, k, x[k] - h))) / (2 * h)
    }, 0)
  }
  fitted <- has_em(model) && identical(model$ecm$nodes, Q)
  xh <- if (fitted) model$ecm$x else suppressWarnings(ecm.rtirt(model, nodes = Q, se = FALSE))$ecm$x
  drop_l <- 0
  fisher <- 0
  for (r in seq_len(n_start)) {
    x <- xh + stats::rnorm(length(xh), 0, 0.3) # jittered on the unconstrained scale
    p <- em_unpack(x, E)
    es <- em_qr_estep(p, E)
    x_new <- em_pack(em_mstep_items(p, E, es, em_moments(es)), E)
    x_new[-ab] <- em_qr_cm_rt(x_new, es$W, E)
    l0 <- sum(es$ll)
    drop_l <- max(drop_l, (l0 - loglik(x_new)) / max(1, abs(l0)))
    g <- gradient(x)
    fisher <- max(fisher, max(abs(g - em_qr_score(x, E))) / max(1, max(abs(g))))
  }
  add("one iteration never lowers the working log-likelihood (relative drop)", drop_l, 1e-8, "items step, then nlminb() on the RT part")
  add("analytic score = numerical gradient of the working log-likelihood (relative difference)", fisher, 1e-4, "em_qr_score()")
  add("ECM stops at a maximum: largest |gradient|", max(abs(em_qr_score(xh, E))), 1e-2, "")
}
