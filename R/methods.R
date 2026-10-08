# Methods: summary / coef / estimates / scores / reliability / convergence / fit_indices / plot
# (the accessors of the birt package), plus precis / dic / compare_para of the Julia package.

model_label <- function(x) paste0(class(x)[1], if (isTRUE(x$qr)) sprintf(" (quantile %s)", format(x$q)) else "")

#' @export
print.rtirt <- function(x, ...) {
  D <- x$data; C <- x$cond
  cat(model_label(x), "\n")
  cat(sprintf("  %d persons x %d items%s\n", nrow(D$Y), ncol(D$Y),
              if (.model_info[[class(x)[1]]]$X) sprintf("; covariates: %s", paste(colnames(D$X), collapse = ", ")) else ""))
  if (has_em(x)) { print_em(x); return(invisible(x)) }
  if (is.null(x$post)) { cat("  not sampled yet: fit <- gibbs(model), or ecm(model) for maximum likelihood\n"); return(invisible(x)) }
  cat(sprintf("  %d chains x %d iterations (burn-in %d, thin %d)%s; %.1f s\n", C$n_chain, C$n_iter, C$n_burnin, C$n_thin,
              if (isTRUE(x$post$settings$collapse)) ", collapsed" else "", x$post$secs))
  rel <- reliability(x)
  cat(sprintf("  Reliability: %s\n", paste(sprintf("%s %.3f", names(rel), rel), collapse = ", ")))
  if (!is.null(x$post$fit)) { fi <- x$post$fit$table
    cat(sprintf("  Marginal: DIC %.1f (pD %.1f), WAIC %.1f%s\n", fi$DIC, fi$pD, fi$WAIC, if (is.na(fi$LOOIC)) "" else sprintf(", LOOIC %.1f", fi$LOOIC))) }
  rh_warn(x)
  cat("Use summary(), coef(), scores(), estimates(), convergence(), fit_indices(), plot().\n")
  invisible(x)
}

rh_warn <- function(x) {
  mr <- suppressWarnings(max(estimates(x)$rhat, na.rm = TRUE))
  if (is.finite(mr) && mr > 1.05) cat(sprintf("  Warning: max R-hat %.2f > 1.05; see convergence()\n", mr))
}

need_post <- function(x) if (is.null(x$post)) stop("run gibbs() first", call. = FALSE)
post_mcmc <- function(x) {                         # draws after burn-in, as an mcmc.list
  b <- x$post$burnin; th <- x$cond$n_thin
  coda::mcmc.list(lapply(x$post$draws, function(m) coda::mcmc(m[seq_len(nrow(m)) * th > b, , drop = FALSE], thin = th)))
}

#' Parameter estimates
#'
#' `estimates()` returns the posterior summary of the item and structural parameters as a
#' data frame: posterior mean (`est`), `sd`, 2.5% and 97.5% quantiles (`q025`, `q975`), R-hat and effective sample size (after
#' burn-in, over the independent chains); `sig = "*"` when the 95% interval excludes 0.
#' `precis()` (the Julia name) prints the same table. `coef()` returns the posterior means
#' as a named vector. For a model fitted by [ecm()] the table has the maximum likelihood
#' estimates with standard errors (`est`, `se`, `z`, Wald 95% interval `q025`, `q975`).
#' @param object A sampled model, or one fitted by [ecm()].
#' @param pars Optional regular expression selecting parameters.
#' @param method `"auto"` (the Gibbs results when the model was sampled, else the ECM fit),
#'   `"gibbs"` or `"ecm"`.
#' @param diagnostics `"classic"`: R-hat of Gelman & Rubin (1992) and the effective sample size of
#'   coda. `"rank"`: the rank-normalized split R-hat and the bulk and tail effective sample sizes of
#'   Vehtari, Gelman, Simpson, Carpenter & Bürkner (2021) (package posterior), which also catch
#'   chains with different variances or poorly mixing tails; the table then has `ess` (bulk) and
#'   `ess_tail`.
#' @param digits Decimals to print.
#' @param ... Unused.
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using
#'   multiple sequences. *Statistical Science, 7*, 457-472.
#'
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization,
#'   folding, and localization: An improved R-hat for assessing convergence of MCMC.
#'   *Bayesian Analysis, 16*, 667-718.
#' @export
estimates <- function(object, ...) UseMethod("estimates")

#' @rdname estimates
#' @export
estimates.rtirt <- function(object, pars = NULL, method = c("auto", "gibbs", "ecm"), diagnostics = c("classic", "rank"), ...) {
  method <- match.arg(method); diagnostics <- match.arg(diagnostics)
  if (method == "ecm" || (method == "auto" && has_em(object))) {
    need_em(object); e <- object$ecm$estimates
    return(if (is.null(pars)) e else e[grepl(pars, e$parameter), , drop = FALSE])
  }
  need_post(object)
  if (!is.null(object$post$estimates) && is.null(pars) && diagnostics == "classic") return(object$post$estimates)
  mc <- post_mcmc(object); M <- as.matrix(mc)
  if (!is.null(pars)) { keep <- grepl(pars, colnames(M)); M <- M[, keep, drop = FALSE]; mc <- mc[, keep, drop = FALSE] }
  sdv <- apply(M, 2, stats::sd); ok <- sdv > 0
  rh <- rep(NA_real_, ncol(M))
  if (coda::nchain(mc) > 1 && any(ok))
    rh[ok] <- tryCatch(coda::gelman.diag(mc[, ok, drop = FALSE], autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
                       error = function(e) NA_real_)
  q <- apply(M, 2, stats::quantile, c(0.025, 0.975))
  if (diagnostics == "rank") {
    if (!requireNamespace("posterior", quietly = TRUE)) stop("diagnostics = \"rank\" needs the 'posterior' package", call. = FALSE)
    arr <- function(j) sapply(mc, function(ch) as.numeric(ch[, j]))          # iterations x chains
    d <- t(vapply(seq_len(ncol(M)), function(j) if (!ok[j]) rep(NA_real_, 3) else { x <- arr(j)
      c(posterior::rhat(x), posterior::ess_bulk(x), posterior::ess_tail(x)) }, numeric(3)))
    return(data.frame(parameter = colnames(M), est = colMeans(M), sd = sdv, q025 = q[1, ], q975 = q[2, ], rhat = d[, 1],
                      ess = as.integer(round(d[, 2])), ess_tail = as.integer(round(d[, 3])), sig = ifelse(q[1, ] > 0 | q[2, ] < 0, "*", ""),
                      row.names = NULL))
  }
  data.frame(parameter = colnames(M), est = colMeans(M), sd = sdv, q025 = q[1, ], q975 = q[2, ],
             rhat = rh, ess = as.integer(round(tryCatch(coda::effectiveSize(mc), error = function(e) rep(NA_real_, ncol(M))))), sig = ifelse(q[1, ] > 0 | q[2, ] < 0, "*", ""),
             row.names = NULL)
}

#' @rdname estimates
#' @export
precis <- function(object, ...) UseMethod("precis")

#' @rdname estimates
#' @export
precis.rtirt <- function(object, pars = NULL, digits = 3, ...) {
  tab <- estimates(object, pars)
  print_table(tab, digits)
  invisible(tab)
}

#' @rdname estimates
#' @export
coef.rtirt <- function(object, ...) {
  e <- estimates(object)
  stats::setNames(e$est, e$parameter)
}

# prior of every parameter of estimates() (as in blavaan's summary), from the samplers in src/
prior_strings <- function(object) {
  info <- .model_info[[class(object)[1]]]; eng <- info$engine; S <- object$settings
  qr <- isTRUE(object$qr); fixed <- identical(S$speed_var, "fixed")
  nm <- par_names(eng, object$data$items, colnames(object$post$settings$X %||% object$data$X) %||% character())
  pr <- stats::setNames(rep("", length(nm)), nm)
  set <- function(pat, v) pr[grepl(pat, nm)] <<- v
  tx <- prior_text(fit_priors(object), object$data$log_t)
  set("^a\\[", if (identical(S$itemtype, "1pl")) "fixed 1" else tx$a)
  set("^b\\[", tx$b)
  if (info$rt) {
    set("^lambda\\[", tx$lambda)
    set("^sigma2t\\[", if (qr && eng == "cross") paste(tx$sigma2t, "(ALD scale)") else tx$sigma2t)
  }
  set("^(beta|beta_ability|beta_speed)\\[", tx$beta)
  switch(eng,
    rtirt = { set("^cor_ability_speed$", if (fixed) "U(-1, 1)" else sprintf("c / sqrt(c^2 + v): c ~ %s, v ~ %s", tx$cov, tx$var_speed))
              set("^var_speed$", if (fixed) "fixed 1" else "c^2 + v") },
    latent = { set("^b_ability$", tx$beta); set("^var_speed$", if (qr) paste(tx$var_speed, "(ALD scale)") else tx$var_speed) },
    cross = { set("^rho\\[", tx$rho); set("^var_speed$", if (fixed) "fixed 1" else tx$var_speed) })
  pr
}

#' Summary of a sampled model
#'
#' Item parameters (posterior means), the structural parameters with posterior SD,
#' 95% interval, R-hat and effective sample size, the person covariance, reliability and DIC.
#' @param object A sampled model.
#' @param digits Decimals to print.
#' @param ... Unused.
#' @return An `rtirt_summary`: a list of the tables (see [rtirt_summary]), printed.
#' @export
summary.rtirt <- function(object, digits = 3, ...) {
  if (has_em(object)) return(summary_em(object, digits))
  need_post(object)
  e <- estimates(object)
  pr <- prior_strings(object); e$prior <- unname(pr[e$parameter]); e$prior[is.na(e$prior)] <- ""
  ald <- isTRUE(object$qr) && class(object)[1] == "rtirt_cross"
  p <- rtirt_partable(e, TRUE, ald, "posterior mean, SD and 95% interval; prior of each parameter")
  w <- rtirt_item_wide(p, object$data$items)
  fi <- fit_indices(object)
  S <- object$post$mean$sigma_p
  structure(list(header = rtirt_header(object),
                 fit = rtable(cbind(N = nrow(object$data$Y), max_rhat = suppressWarnings(max(e$rhat, na.rm = TRUE)), fi),
                              paste0("marginal likelihood (persons integrated out); DIC_complete conditions on the persons",
                                     if (isTRUE(object$qr)) "; ALD working likelihood: not comparable with mean models or across quantile levels" else "")),
                 items = rtable(w$est, "item parameters (posterior means)"),
                 se = rtable(w$se, "posterior SDs of the item parameters"),
                 parameters = p,
                 covariance = if (!is.null(S)) rtable(data.frame(variable = rownames(S), S, check.names = FALSE), "person covariance (posterior mean)"),
                 reliability = rtable(as.data.frame(as.list(reliability(object))), "reliability of the EAP scores: var(EAP) / (var(EAP) + mean PSD^2)")),
            class = c("summary.rtirt", "rtirt_summary"))
}

#' Person scores
#'
#' Posterior means (EAP) and posterior SDs of ability and speed, one row per person.
#' @param object A sampled model.
#' @param ... Unused.
#' @return A data frame with `id`, `ability`, `ability_psd` and (RT models) `speed`, `speed_psd`.
#' @export
scores <- function(object, ...) UseMethod("scores")

#' @rdname scores
#' @export
scores.rtirt <- function(object, ...) {
  p <- if (has_em(object)) object$ecm$person else { need_post(object); object$post$person }
  out <- data.frame(id = object$data$id, ability = p$theta$mean, ability_psd = p$theta$sd)
  if (!is.null(p$zeta)) { out$speed <- p$zeta$mean; out$speed_psd <- p$zeta$sd }
  out
}

#' Empirical reliability
#'
#' var(EAP) / (var(EAP) + mean PSD^2) of ability and speed.
#' @param object A sampled model.
#' @param ... Unused.
#' @export
reliability <- function(object, ...) UseMethod("reliability")

#' @rdname reliability
#' @export
reliability.rtirt <- function(object, ...) {
  r <- function(p) { v <- stats::var(p$mean); v / (v + mean(p$sd^2)) }
  p <- if (has_em(object)) object$ecm$person else { need_post(object); object$post$person }
  c(ability = r(p$theta), if (!is.null(p$zeta)) c(speed = r(p$zeta)))
}

#' Fit indices
#'
#' Criteria from the marginal log-likelihood (person parameters integrated out; see
#' [marginal_loglik()]), which compare models with different latent structures:
#' `DIC` (with `pD` = mean deviance - deviance at the posterior means of the item and
#' structural parameters), `WAIC` (`p_waic`) and PSIS-LOO (`LOOIC`, `SE_LOOIC`; needs the
#' 'loo' package). `DIC_complete` is the complete-data DIC of [dic()] (conditional on the
#' person parameters; for comparisons within one model only). Quantile models use the ALD as
#' a working likelihood, so their criteria are not comparable with mean models or across
#' quantile levels. For a model fitted by [ecm()]: the log-likelihood, AIC and BIC.
#' @param object A sampled model, or one fitted by [ecm()].
#' @param ... Unused.
#' @export
fit_indices <- function(object, ...) UseMethod("fit_indices")

#' @rdname fit_indices
#' @export
fit_indices.rtirt <- function(object, ...) {
  if (has_em(object)) { f <- object$ecm
    return(data.frame(method = if (f$map) "ECM (MAP items)" else "ECM", logLik = f$logLik, npar = f$npar, N = f$nobs,
                      AIC = f$AIC, BIC = f$BIC)) }
  need_post(object)
  m <- object$post$fit %||% { message("computing the marginal log-likelihood ..."); marginal_fit(object) }
  cbind(m$table, DIC_complete = dic(object)$dic)
}

# marginal-likelihood criteria (stored by gibbs(fit_indices = TRUE))
marginal_fit <- function(object) {
  L <- marginal_loglik(object)
  M <- post_draws(object, Inf); xn <- colnames(object$post$settings$X) %||% character()
  pm <- colMeans(M)
  d_hat <- -2 * sum(marg_ll_draw(object, draw_pars(pm, object$data$items, xn),
                                 Q = if (isTRUE(object$qr) && .model_info[[class(object)[1]]]$engine == "cross") 61 else 41))
  d_bar <- -2 * mean(rowSums(L))
  mx <- apply(L, 2, max); lppd <- sum(mx + log(colMeans(exp(sweep(L, 2, mx)))))
  p_waic <- sum(apply(L, 2, stats::var))
  tab <- data.frame(Deviance = d_bar, pD = d_bar - d_hat, DIC = 2 * d_bar - d_hat, WAIC = -2 * (lppd - p_waic), p_waic = p_waic,
                    LOOIC = NA_real_, SE_LOOIC = NA_real_)
  lo <- NULL
  if (requireNamespace("loo", quietly = TRUE)) {
    lo <- suppressWarnings(loo::loo(L))
    tab$LOOIC <- lo$estimates["looic", "Estimate"]; tab$SE_LOOIC <- lo$estimates["looic", "SE"]
  }
  list(table = tab, loglik = L, loo = lo)
}

#' Deviance information criterion
#'
#' As `getDic(MCMC)` in Julia: \eqn{\bar D} is the posterior mean of the complete-data deviance
#' (after burn-in; the Julia version averaged over the burn-in as well), \eqn{\hat D} the
#' deviance at the posterior means, \eqn{p_D = \bar D - \hat D}. The deviance conditions on the
#' person parameters (and, for the ALD models, on the mixing variables), so DIC values are
#' not comparable between mean and quantile models or across quantile levels.
#' @param object A sampled model.
#' @param ... Unused.
#' @export
dic <- function(object, ...) UseMethod("dic")

#' @rdname dic
#' @export
dic.rtirt <- function(object, ...) {
  need_post(object)
  ll <- object$post$loglik[(object$post$burnin + 1):object$cond$n_iter, , drop = FALSE]
  d_bar <- -2 * mean(ll); d_hat <- -2 * complete_loglik(object)
  data.frame(d_bar = d_bar, d_hat = d_hat, p_d = d_bar - d_hat, dic = d_bar + (d_bar - d_hat))
}

# complete-data log-likelihood at the posterior means (mirrors the loglik in the C++ samplers)
complete_loglik <- function(object) {
  m <- object$post$mean; D <- object$data; X <- object$post$settings$X; engine <- .model_info[[class(object)[1]]]$engine
  N <- nrow(D$Y); K <- ncol(D$Y)
  ll <- sum(calc_loglik_nodes(D$Y, m$a, m$b, matrix(m$theta)))           # theta as one person-specific node
  ln <- function(x, mu, v) sum(stats::dnorm(x, mu, sqrt(v), log = TRUE))
  q <- object$post$settings$q
  if (engine == "mlirt") return(ll + ln(m$theta, drop(X %*% m$beta), 1))
  if (engine == "cross") {
    mu <- matrix(m$lambda, N, K, byrow = TRUE) - m$zeta - outer(m$theta, m$rho)
    v <- matrix(m$sigma2t, N, K, byrow = TRUE)
    if (object$qr) { k <- ald_k(q); mu <- mu + k[1] * m$nu; v <- k[2] * v * m$nu }
    return(ll + ln(D$log_t, mu, v) + ln(m$theta, 0, 1) + ln(m$zeta, 0, m$sigma_p[2, 2]))
  }
  ll <- ll + ln(D$log_t, matrix(m$lambda, N, K, byrow = TRUE) - m$zeta, matrix(m$sigma2t, N, K, byrow = TRUE))
  if (engine == "latent") {
    mz <- drop(X %*% m$beta) + m$b_theta * m$theta; s <- m$sigma_p[2, 2]
    if (object$qr) { k <- ald_k(q); return(ll + ln(m$theta, 0, 1) + ln(m$zeta, mz + k[1] * m$nu, k[2] * s * m$nu)) }
    return(ll + ln(m$theta, 0, 1) + ln(m$zeta, mz, s))
  }
  S <- m$sigma_p; mu1 <- drop(X %*% m$beta[, 1]); mu2 <- drop(X %*% m$beta[, 2])
  cc <- S[1, 2]; v <- S[2, 2] - cc^2
  ll + ln(m$theta, mu1, 1) + ln(m$zeta, mu2 + cc * (m$theta - mu1), v)
}

#' Convergence check
#'
#' The share of item and structural parameters with effective sample size above `ess` and
#' R-hat below `rhat` (over independent chains), and the worst parameters. With
#' `diagnostics = "rank"` (see [estimates()]) both the bulk and the tail ESS must exceed `ess`,
#' and the default R-hat threshold is 1.01 (Vehtari et al., 2021) instead of 1.05.
#' Julia: `checkConvergence()`.
#' @param object A sampled model.
#' @param ess,rhat Thresholds (`rhat` defaults to 1.05, or 1.01 with `diagnostics = "rank"`).
#' @param diagnostics `"classic"` or `"rank"`, as in [estimates()].
#' @param ... Unused.
#' @export
convergence <- function(object, ...) UseMethod("convergence")

#' @rdname convergence
#' @export
convergence.rtirt <- function(object, ess = 400, rhat = NULL, diagnostics = c("classic", "rank"), ...) {
  if (has_em(object)) {                                                # an ECM fit: the iterations, not chains
    f <- object$ecm
    cat(sprintf("ECM: %s in %d iterations (relative change of the objective below tol)\n", if (f$converged) "converged" else "NOT converged", f$iterations))
    return(invisible(list(converged = f$converged, iterations = f$iterations)))
  }
  diagnostics <- match.arg(diagnostics); rhat <- rhat %||% if (diagnostics == "rank") 1.01 else 1.05
  p <- estimates(object, diagnostics = diagnostics)
  tail_ok <- if (is.null(p$ess_tail)) TRUE else p$ess_tail > ess
  ok <- p$ess > ess & (is.na(p$rhat) | p$rhat < rhat) & tail_ok
  cols <- intersect(c("parameter", "rhat", "ess", "ess_tail"), names(p))
  out <- list(share = mean(ok), worst = p[order(-p$rhat), cols][1:min(5, nrow(p)), ])
  cat(sprintf("%.1f%% of %d parameters with ESS%s > %d and R-hat%s < %.2f\n", 100 * out$share, nrow(p),
              if (diagnostics == "rank") " (bulk and tail)" else "", ess, if (diagnostics == "rank") " (rank-normalized)" else "", rhat))
  print_table(out$worst)
  invisible(out)
}

#' Compare estimates with true values
#'
#' As `comparePara()` in Julia (for simulations with `true_para`): bias, RMSE and correlation
#' of the posterior means (or of the ECM estimates) for one parameter block.
#' @param object A sampled model with `true_para`.
#' @param par Parameter name, e.g. `"a"`, `"b"`, `"ability"`, `"lambda"`, `"rho"`, `"b_ability"`,
#'   `"beta"` (the Julia names `"theta"`, `"b_theta"` and, for `rtirt_latent`, `"rho"` are accepted).
#' @export
compare_para <- function(object, par = "a") {
  if (!has_em(object)) need_post(object)
  if (is.null(object$true_para)) stop("no true_para in the model")
  lab <- par
  if (par == "rho" && inherits(object, "rtirt_latent")) par <- "b_theta"
  par <- switch(par, ability = "theta", b_ability = "b_theta", par)       # internal (Julia) names
  src <- if (has_em(object)) c(object$ecm$par, list(theta = object$ecm$person$theta$mean, zeta = object$ecm$person$zeta$mean)) else object$post$mean
  if (is.null(src[[par]])) stop(sprintf("'%s' is not among the estimates of this fit", lab), call. = FALSE)
  est <- as.numeric(src[[par]]); tru <- as.numeric(object$true_para[[par]])
  if (length(est) != length(tru)) stop(sprintf("'%s': %d estimates vs %d true values", par, length(est), length(tru)))
  data.frame(par = lab, n = length(est), bias = get_bias(est, tru), rmse = get_rmse(est, tru),
             corr = if (length(est) > 2 && stats::sd(est) > 0 && stats::sd(tru) > 0) stats::cor(est, tru) else NA_real_)
}

#' @importFrom coda as.mcmc.list
#' @export
as.mcmc.list.rtirt <- function(x, ...) { need_post(x); post_mcmc(x) }

#' Plots of a sampled model
#'
#' `type = "items"`: posterior means and 95% intervals of the item parameters (ECM fits: estimates
#' and Wald 95% intervals), one panel
#' per parameter type (a, b, lambda, sigma2t, rho). `type = "trace"`: trace plots (all
#' iterations; burn-in shaded) of the structural parameters or of those matching `pars`.
#' @param x A sampled model.
#' @param type `"items"` or `"trace"`.
#' @param pars Regular expression selecting parameters (`"trace"`; at most 12 are drawn).
#' @param ... Passed to [graphics::plot()].
#' @export
plot.rtirt <- function(x, type = c("items", "trace"), pars = NULL, ...) {
  type <- match.arg(type)
  if (type == "trace" || !has_em(x)) need_post(x)                      # items: posterior or Wald 95% intervals (ECM)
  e <- estimates(x)
  if (type == "items") {
    blocks <- intersect(c("a", "b", "lambda", "sigma2t", "rho"), unique(sub("\\[.*", "", e$parameter[grepl("\\[", e$parameter)])))
    blocks <- blocks[vapply(blocks, function(b) sum(startsWith(e$parameter, paste0(b, "["))) == ncol(x$data$Y), TRUE)]
    op <- graphics::par(mfrow = c(1, length(blocks)), mar = c(4, if (length(blocks)) 7 else 4, 2.5, 1)); on.exit(graphics::par(op))
    for (b in blocks) {
      r <- e[startsWith(e$parameter, paste0(b, "[")), ]; K <- nrow(r)
      graphics::plot(r$est, K:1, xlim = range(r$q025, r$q975, if (b == "rho") 0), yaxt = "n", pch = 19,
                     xlab = b, ylab = "", main = b, ...)
      graphics::segments(r$q025, K:1, r$q975, K:1)
      if (b == "rho") graphics::abline(v = 0, lty = 3)
      if (b == blocks[1]) graphics::axis(2, K:1, x$data$items, las = 1, cex.axis = 0.8)
    }
    return(invisible(e))
  }
  if (is.null(pars)) pars <- "^(rho|beta|b_ability|cor_ability_speed|var_speed)"
  sel <- e$parameter[grepl(pars, e$parameter)]
  if (!length(sel)) stop("no parameters match '", pars, "'")
  sel <- utils::head(sel, 12)
  nc <- ceiling(sqrt(length(sel)))
  op <- graphics::par(mfrow = c(ceiling(length(sel) / nc), nc), mar = c(3, 3, 2, 1)); on.exit(graphics::par(op))
  for (p in sel) {
    D <- sapply(x$post$draws, function(m) m[, p]); it <- seq_len(nrow(D)) * x$cond$n_thin
    graphics::matplot(it, D, type = "l", lty = 1, col = seq_len(ncol(D)), xlab = "", ylab = "", main = p, ...)
    graphics::rect(0, graphics::par("usr")[3], x$post$burnin, graphics::par("usr")[4], col = grDevices::adjustcolor("grey", 0.3), border = NA)
  }
  invisible(e[e$parameter %in% sel, ])
}

# ---- several quantile levels ---------------------------------------------------------------

#' Several quantile levels
#'
#' `rtirt_cross(data, quantile = c(.1, .5, .9))` (or `rtirt_latent`) builds one ALD model per
#' level, an `rtirt_qset`; [gibbs()] samples all of them. `coef()` gives a parameters x
#' levels matrix, `estimates()` stacks the tables with a column `q`, and `plot()` draws the
#' quantile coefficients (cross-relations `rho`, or the regression coefficients of
#' `rtirt_latent`) against the level with 95% intervals.
#' @param x,object An `rtirt_qset`.
#' @param pars Regular expression of the parameters to plot (default: `rho` for
#'   `rtirt_cross`, `beta` and `b_ability` for `rtirt_latent`).
#' @param ... Unused, or passed to [graphics::plot()].
#' @name rtirt_qset
NULL

qset_sampled <- function(x) all(vapply(x$models, function(m) !is.null(m$post), TRUE))

#' @rdname rtirt_qset
#' @export
print.rtirt_qset <- function(x, ...) {
  cat(sprintf("%s at quantiles %s\n", x$type, paste(format(x$quantile), collapse = ", ")))
  if (!qset_sampled(x)) { cat("  not sampled yet: fits <- gibbs(models)\n"); return(invisible(x)) }
  cf <- stats::coef(x); key <- grepl(if (x$type == "rtirt_cross") "^rho\\[" else "^(beta|b_ability)", rownames(cf))
  print_table(data.frame(parameter = rownames(cf)[key], cf[key, , drop = FALSE], check.names = FALSE))
  for (m in x$models) rh_warn(m)
  invisible(x)
}

#' @rdname rtirt_qset
#' @export
coef.rtirt_qset <- function(object, ...) {
  out <- sapply(object$models, stats::coef)
  colnames(out) <- sprintf("q=%s", format(object$quantile)); out
}

#' @rdname rtirt_qset
#' @export
estimates.rtirt_qset <- function(object, ...)
  do.call(rbind, Map(function(m, q) cbind(q = q, estimates(m)), object$models, object$quantile))

#' @rdname rtirt_qset
#' @export
summary.rtirt_qset <- function(object, ...) { print(object, ...); invisible(object) }

#' @rdname rtirt_qset
#' @export
plot.rtirt_qset <- function(x, pars = NULL, ...) {
  if (!qset_sampled(x)) stop("run gibbs() first", call. = FALSE)
  if (is.null(pars)) pars <- if (x$type == "rtirt_cross") "^rho\\[" else "^(beta|b_ability)"
  e <- estimates(x); e <- e[grepl(pars, e$parameter), ]
  key <- unique(e$parameter); if (!length(key)) stop("no parameters match '", pars, "'")
  nc <- ceiling(sqrt(length(key)))
  op <- graphics::par(mfrow = c(ceiling(length(key) / nc), nc), mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
  for (k in key) {
    r <- e[e$parameter == k, ]; r <- r[order(r$q), ]
    graphics::plot(r$q, r$est, type = "n", ylim = range(r$q025, r$q975, 0), xlab = "quantile", ylab = "", main = k, ...)
    graphics::polygon(c(r$q, rev(r$q)), c(r$q025, rev(r$q975)), col = "grey85", border = NA)
    graphics::lines(r$q, r$est, type = "b", pch = 19); graphics::abline(h = 0, lty = 3)
  }
  invisible(e)
}

#' @export
gibbs.default <- function(model, ...) stop("model must be a model object, e.g. rtirt_cross(data), or a block_model()", call. = FALSE)

# plain table (as in birt): numbers right-aligned with fixed decimals, text left-aligned
print_table <- function(df, digits = 3) {
  fmt <- function(x) {
    if (!is.numeric(x)) return(ifelse(is.na(x), "", as.character(x)))
    if (is.integer(x)) return(ifelse(is.na(x), "", formatC(x, format = "d", big.mark = "")))
    ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
  }
  cells <- vapply(df, fmt, character(nrow(df)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = nrow(df))
  head <- names(df)
  w <- pmax(nchar(head), apply(cells, 2, function(x) max(nchar(x), 0)))
  left <- !vapply(df, is.numeric, TRUE)
  row <- function(x) sub("\\s+$", "", paste(ifelse(left, sprintf("%-*s", w, x), sprintf("%*s", w, x)), collapse = "  "))
  cat(c(row(head), apply(cells, 1, row)), sep = "\n")
}
