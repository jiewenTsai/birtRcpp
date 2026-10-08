# Marginal log-likelihood: the person parameters integrated out ---------------------------
#
# Per person j and posterior draw d, log p(y_j, log t_j | item and structural parameters),
# integrating ability theta by Gauss-Hermite quadrature on its N(mu_j, 1) distribution and
# speed zeta
#   - analytically for normal log response times: given theta, speed is normal a posteriori
#     (completing the square, normal_given_sums(); the same integral as in ecm());
#   - in closed form for rtirt_latent quantile models: the product of the normal RT
#     likelihoods is a normal density in zeta, and the convolution of a normal with an ALD
#     has a closed form;
#   - in closed form for rtirt_cross quantile models: with ALD residuals u_k = b_k + zeta the
#     log-likelihood is piecewise linear in zeta (kinks at -b_k), so the integral against
#     N(0, var_speed) is a sum of K + 1 normal-CDF differences.
# This is the likelihood DIC, WAIC and PSIS-LOO need to compare models with different latent
# structures (Merkle, Furr & Rabe-Hesketh, 2019); the complete-data DIC conditions on the
# person parameters.

# draws x parameters matrix after burn-in, thinned to at most n_draws rows
post_draws <- function(object, n_draws) {
  M <- as.matrix(post_mcmc(object))
  if (nrow(M) > n_draws) M <- M[unique(round(seq(1, nrow(M), length.out = n_draws))), , drop = FALSE]
  M
}

draw_pars <- function(row, items, xn) {
  g <- function(p, n) { v <- row[sprintf("%s[%s]", p, n)]; if (anyNA(v)) NULL else unname(v) }
  s <- function(p) if (p %in% names(row)) unname(row[[p]]) else NULL
  list(a = g("a", items) %||% rep(1, length(items)), b = g("b", items), lambda = g("lambda", items),
       sigma2t = g("sigma2t", items), rho = g("rho", items), beta = if (length(xn)) g("beta", xn),
       beta_theta = if (length(xn)) g("beta_ability", xn), beta_speed = if (length(xn)) g("beta_speed", xn),
       b_theta = s("b_ability"), var_speed = s("var_speed"), cor = s("cor_ability_speed"))
}

# log of the convolution  int N(z; zh, tau2) ALD(z; m, s, q) dz   (ALD scale s: rho_q(u) / s)
log_norm_ald <- function(zh, tau2, m, s, q) {
  tau <- sqrt(tau2); d <- zh - m; a <- q / s; b <- (1 - q) / s
  lp <- -a * d + a^2 * tau2 / 2 + stats::pnorm((d - a * tau2) / tau, log.p = TRUE)
  ln <- b * d + b^2 * tau2 / 2 + stats::pnorm((-d - b * tau2) / tau, log.p = TRUE)
  mx <- pmax(lp, ln)
  log(q * (1 - q) / s) + mx + log(exp(lp - mx) + exp(ln - mx))
}

marg_ll_draw <- function(object, p, Q = 41) {
  D <- object$data; X <- object$post$settings$X; info <- .model_info[[class(object)[1]]]
  eng <- info$engine; qr <- isTRUE(object$qr); q <- object$post$settings$q
  Y <- D$Y; N <- nrow(Y); K <- ncol(Y)
  nd <- calc_nodes(Q)
  xb <- function(beta) if (length(beta) && ncol(X)) drop(X %*% beta) else rep(0, N)
  mu_t <- switch(eng, mlirt = xb(p$beta), rtirt = xb(p$beta_theta), rep(0, N))
  TH <- outer(mu_t, nd$z, "+")                                              # ability at the nodes, N x Q
  L <- calc_loglik_nodes(Y, p$a, p$b, if (all(mu_t == 0)) nd$z else TH)
  if (eng != "mlirt") {
    C <- D$log_t - matrix(p$lambda, N, K, byrow = TRUE)                     # c = log t - lambda
    w <- 1 / p$sigma2t; rr <- if (eng == "cross") p$rho else rep(0, K); vs <- p$var_speed
    # normal RT part: speed u ~ N(m, v) seen through o_i = -(c_i + rho_i theta) = u + e_i, sums of w o and w o^2
    sp <- function(m, v) normal_given_sums(-(drop(C %*% w) + TH * sum(w * rr)),
                                           drop(C^2 %*% w) + 2 * TH * drop(C %*% (w * rr)) + TH^2 * sum(w * rr^2), w, m, v)$logI
    L <- L + switch(eng,
      rtirt = { cc <- p$cor * sqrt(vs); sp(xb(p$beta_speed) + cc * (TH - mu_t), vs - cc^2) },
      latent = {
        m <- xb(p$beta) + p$b_theta * TH
        if (!qr) sp(m, vs)
        else {          # prod_i N(c_i + zeta; 0, sigma2_i) = exp(cst - (S2 - S1^2 / A) / 2) sqrt(2 pi / A) N(zeta; -S1 / A, 1 / A)
          A <- sum(w); S1 <- drop(C %*% w); S2 <- drop(C^2 %*% w)
          -0.5 * K * log(2 * pi) - 0.5 * sum(log(p$sigma2t)) - 0.5 * (S2 - S1^2 / A) + 0.5 * log(2 * pi / A) + log_norm_ald(-S1 / A, 1 / A, m, vs, q)
        }
      },
      cross = if (!qr) sp(0, vs) else vapply(seq_len(Q), function(k) .log_ald_rt(C + outer(TH[, k], p$rho), p$sigma2t, vs, q), numeric(N)))
  }
  calc_node_weights(L, nd$logw)$ll
}

#' Marginal log-likelihood and fit indices
#'
#' `marginal_loglik()` returns the log-likelihood of each person with the person parameters
#' (ability, speed) integrated out, for posterior draws of the item and structural
#' parameters. `fit_indices()` uses it for DIC, WAIC and PSIS-LOO (with the 'loo' package);
#' these, unlike the complete-data DIC of [dic()], compare models with different latent
#' structures (Merkle, Furr & Rabe-Hesketh, 2019). Ability is integrated by Gauss-Hermite
#' quadrature and speed in closed form (normal log times; normal-ALD convolution for
#' `rtirt_latent` quantile models; piecewise-linear ALD log-likelihood for `rtirt_cross`). Quantile models use the ALD as a
#' working likelihood, so their criteria are not comparable with mean models or across levels.
#' @param object A sampled model.
#' @param n_draws Number of posterior draws (after burn-in, evenly spaced; default 200, 100 for
#'   `rtirt_cross` quantile models).
#' @param Q Gauss-Hermite nodes for ability (default 41; 61 for `rtirt_cross` quantile models,
#'   whose integrand in ability has sharper features: about 3e-4 per person from the exact value).
#' @return A draws x persons matrix.
#' @export
marginal_loglik <- function(object, n_draws = NULL, Q = NULL) {
  need_post(object)
  cq <- isTRUE(object$qr) && .model_info[[class(object)[1]]]$engine == "cross"
  if (is.null(Q)) Q <- if (cq) 61 else 41
  if (is.null(n_draws)) n_draws <- if (cq) 100 else 200
  M <- post_draws(object, n_draws); xn <- colnames(object$post$settings$X) %||% character()
  items <- object$data$items
  t(apply(M, 1, function(row) marg_ll_draw(object, draw_pars(row, items, xn), Q)))
}
