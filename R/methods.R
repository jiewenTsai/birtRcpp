# Methods: summary / coef / estimates / scores / reliability / convergence / fit_indices / plot
# (the accessors of the birt package), plus precis / dic / compare_para of the Julia package.

# the header of print() and summary(): the model as named in ecm() and gibbs(), then its class,
# e.g. model = "cross", quantile = 0.9 (rtirt_cross)
model_label <- function(x) {
  quantile <- if (isTRUE(x$qr)) sprintf(", quantile = %s", format(x$q)) else ""
  sprintf("%s%s (%s)", model_call_text(class(x)[1]), quantile, class(x)[1])
}

#' @export
print.rtirt <- function(x, ...) {
  data <- x$data
  cond <- x$cond
  covariates <- if (.model_info[[class(x)[1]]]$X) sprintf("; covariates: %s", paste(colnames(data$X), collapse = ", ")) else ""
  cat(model_label(x), "\n", sep = "")
  cat(sprintf("  %d persons x %d items%s\n", nrow(data$Y), ncol(data$Y), covariates))
  if (has_em(x)) {
    print_em(x)
    return(invisible(x))
  }
  if (is.null(x$post)) {
    cat("  not sampled yet: fit <- gibbs(model), or ecm(model) for maximum likelihood\n")
    return(invisible(x))
  }
  collapsed <- if (isTRUE(x$post$settings$collapse)) ", collapsed" else ""
  cat(sprintf("  %d chains x %d iterations (burn-in %d, thin %d)%s; %.1f s\n", cond$n_chain, cond$n_iter,
              cond$n_burnin, cond$n_thin, collapsed, x$post$secs))
  rel <- reliability(x)
  cat(sprintf("  Reliability: %s\n", paste(sprintf("%s %.3f", names(rel), rel), collapse = ", ")))
  if (!is.null(x$post$fit)) {
    fit <- x$post$fit$table
    looic <- if (is.na(fit$LOOIC)) "" else sprintf(", LOOIC %s", fmt_fixed(fit$LOOIC))
    cat(sprintf("  Marginal: DIC %s (pD %s), WAIC %s%s\n", fmt_fixed(fit$DIC), fmt_fixed(fit$pD), fmt_fixed(fit$WAIC), looic))
  }
  convergence_warn(x)
  cat("Use summary(), coef(), scores(), estimates(), convergence(), fit_indices(), plot().\n")
  invisible(x)
}

# the warning lines of print(): the criteria of convergence() (ESS > 400 and R-hat < 1.05 for
# every parameter that varies), and a note that one chain gives no R-hat; `level` names the
# quantile of a model in an rtirt_qset
convergence_warn <- function(x, level = "") {
  check <- convergence_check(x)
  where <- if (nzchar(level)) sprintf(" (%s)", level) else ""
  if (check$n_chain < 2) {
    cat(sprintf("  Warning%s: 1 chain, so R-hat cannot be computed (it needs n_chain >= 2); see convergence()\n", where))
  }
  if (check$share < 1) {
    max_rhat <- suppressWarnings(max(check$table$rhat, na.rm = TRUE))
    rhat_text <- if (is.finite(max_rhat)) sprintf("max R-hat %s, ", fmt_fixed(max_rhat)) else ""
    cat(sprintf("  Warning%s: only %.1f%% of %d parameters have ESS > %d and R-hat < %s (%smin ESS %d); run longer (n_iter) or see convergence()\n",
                where, 100 * check$share, check$n, check$ess, format(check$rhat), rhat_text, min(check$table$ess, na.rm = TRUE)))
  }
}

need_post <- function(x) {
  if (is.null(x$post)) {
    if (has_em(x)) {
      stop(sprintf("this is an ecm() fit: it has no chains or draws; sample with gibbs(data, %s, init = fit)",
                   model_call_text(class(x)[1])), call. = FALSE)
    }
    stop("run gibbs() first", call. = FALSE)
  }
}
# draws after burn-in, as an mcmc.list
post_mcmc <- function(x) {
  burnin <- x$post$burnin
  thin <- x$cond$n_thin
  after_burnin <- function(m) {
    kept <- seq_len(nrow(m)) * thin > burnin
    coda::mcmc(m[kept, , drop = FALSE], thin = thin)
  }
  coda::mcmc.list(lapply(x$post$draws, after_burnin))
}

#' Parameter estimates
#'
#' `estimates()` returns the posterior summary of the item and structural parameters as a
#' data frame: posterior mean (`est`), `sd`, 2.5% and 97.5% quantiles (`q025`, `q975`), R-hat and effective sample size (after
#' burn-in, over the independent chains); `sig = "*"` when the 95% interval excludes 0. `sig`
#' is given only for parameters where 0 is a meaningful value: difficulties `b`, regression
#' coefficients `beta`, cross-relations `rho`, `cor_ability_speed` and `b_ability`. It is `NA`
#' for discriminations, variances (`sigma2t`, `var_speed`), which are positive and whose
#' intervals never contain 0, and time intensities `lambda`.
#' `precis()` (the Julia name) returns the same table, printed in a plain layout (assigning it
#' prints nothing); `as.data.frame()` gives it as a plain data frame, e.g. for ggplot2. `coef()` returns the posterior means
#' as a named vector. For a model fitted by [ecm()] the table has the maximum likelihood
#' estimates with standard errors (`est`, `se`, `z`, Wald 95% interval `q025`, `q975`); `z`,
#' like `sig`, is `NA` where 0 is not a meaningful value, and `summary()` then shows no z or p.
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
#' @return `estimates()`: a tibble of class [rtirt_tbl] (printed with 3 decimals; `x$est` and
#'   `as.data.frame(x)` keep the full precision). `coef()`: a named vector, printed with 3 decimals.
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using
#'   multiple sequences. *Statistical Science, 7*, 457-472.
#'
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization,
#'   folding, and localization: An improved R-hat for assessing convergence of MCMC.
#'   *Bayesian Analysis, 16*, 667-718.
#' @export
estimates <- function(object, ...) {
  UseMethod("estimates")
}

#' @rdname estimates
#' @export
estimates.rtirt <- function(object, pars = NULL, method = c("auto", "gibbs", "ecm"), diagnostics = c("classic", "rank"), ...) {
  table <- estimates_rtirt(object, pars, method, diagnostics)
  as_rtirt_tbl(tibble::as_tibble(table))
}

# the parameter table of estimates(), as a data frame
estimates_rtirt <- function(object, pars, method, diagnostics) {
  method <- match.arg(method, c("auto", "gibbs", "ecm"))
  diagnostics <- match.arg(diagnostics, c("classic", "rank"))
  if (method == "ecm" || (method == "auto" && has_em(object))) {
    need_em(object)
    table <- object$ecm$estimates
    if (!is.null(pars)) {
      warn_no_match(pars, table$parameter)
      table <- table[grepl(pars, table$parameter), , drop = FALSE]
    }
    return(table)
  }
  need_post(object)
  if (!is.null(object$post$estimates) && is.null(pars) && diagnostics == "classic") {
    return(object$post$estimates)
  }

  chains <- post_mcmc(object)
  draws <- as.matrix(chains)
  if (!is.null(pars)) {
    warn_no_match(pars, colnames(draws))
    keep <- grepl(pars, colnames(draws))
    if (!any(keep)) {
      return(estimates_rtirt(object, NULL, method, diagnostics)[0, , drop = FALSE])
    }
    draws <- draws[, keep, drop = FALSE]
    chains <- chains[, keep, drop = FALSE]
  }
  post_sd <- apply(draws, 2, stats::sd)
  varies <- post_sd > 0
  interval <- apply(draws, 2, stats::quantile, c(0.025, 0.975))
  table <- data.frame(parameter = colnames(draws), est = colMeans(draws), sd = post_sd,
                      q025 = interval[1, ], q975 = interval[2, ], row.names = NULL)
  if (diagnostics == "rank") {
    rank <- rank_diagnostics(chains, varies)
    table$rhat <- rank[, 1]
    table$ess <- as.integer(round(rank[, 2]))
    table$ess_tail <- as.integer(round(rank[, 3]))
  } else {
    table$rhat <- classic_rhat(chains, varies)
    ess <- tryCatch(coda::effectiveSize(chains), error = function(e) rep(NA_real_, ncol(draws)))
    table$ess <- as.integer(round(ess))
  }
  table$sig <- sig_flag(table$parameter, interval[1, ], interval[2, ])
  table
}

# estimates(pars = ) that selects nothing: say so and list the parameter names
warn_no_match <- function(pars, parameters) {
  if (any(grepl(pars, parameters))) {
    return(invisible())
  }
  blocks <- unique(sub("\\[.*", "", parameters))
  warning(sprintf("no parameters match \"%s\"; this fit has %s", pars, paste(blocks, collapse = ", ")), call. = FALSE)
}

# `sig` of estimates(): "*" when the 95% interval excludes 0, for the parameters where 0 is a
# meaningful value (difficulties, regression coefficients, cross-relations, the correlation);
# NA for the others: discriminations and variances are positive (their intervals, on the log
# scale, never contain 0) and a time intensity of 0 (one second) is no hypothesis
zero_null <- function(parameter) {
  grepl("^(b|beta|beta_ability|beta_speed|rho)\\[", parameter) | parameter %in% c("cor_ability_speed", "b_ability")
}

sig_flag <- function(parameter, lower, upper) {
  flag <- ifelse(!is.na(lower) & (lower > 0 | upper < 0), "*", "")
  flag[!zero_null(parameter)] <- NA_character_
  flag
}

# Gelman-Rubin R-hat of the parameters that vary (NA otherwise, and with one chain)
classic_rhat <- function(chains, varies) {
  rhat <- rep(NA_real_, length(varies))
  if (coda::nchain(chains) > 1 && any(varies)) {
    rhat[varies] <- tryCatch(
      coda::gelman.diag(chains[, varies, drop = FALSE], autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
      error = function(e) NA_real_)
  }
  rhat
}

# rank-normalized R-hat, bulk ESS and tail ESS (Vehtari et al., 2021), one row per parameter
rank_diagnostics <- function(chains, varies) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("diagnostics = \"rank\" needs the 'posterior' package", call. = FALSE)
  }
  one_parameter <- function(j) {
    if (!varies[j]) {
      return(rep(NA_real_, 3))
    }
    draws <- sapply(chains, function(chain) as.numeric(chain[, j]))          # iterations x chains
    c(posterior::rhat(draws), posterior::ess_bulk(draws), posterior::ess_tail(draws))
  }
  t(vapply(seq_along(varies), one_parameter, numeric(3)))
}

#' @rdname estimates
#' @export
precis <- function(object, ...) {
  UseMethod("precis")
}

#' @rdname estimates
#' @export
precis.rtirt <- function(object, pars = NULL, digits = 3, ...) {
  tab <- estimates(object, pars)
  structure(tab, class = c("rtirt_precis", class(tab)), digits = digits)
}

#' @rdname estimates
#' @param x The result of `precis()`.
#' @export
print.rtirt_precis <- function(x, digits = attr(x, "digits"), ...) {
  table <- x
  attr(table, "digits") <- NULL
  class(table) <- "data.frame"
  print_table(table, digits %||% 3)
  invisible(x)
}

#' @rdname estimates
#' @param row.names,optional Unused.
#' @export
as.data.frame.rtirt <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(estimates(x))
}

#' @rdname estimates
#' @export
coef.rtirt <- function(object, ...) {
  table <- estimates(object)
  as_rtirt_values(stats::setNames(table$est, table$parameter))
}

# prior of every parameter of estimates() (as in blavaan's summary), from the samplers in src/
prior_strings <- function(object) {
  info <- .model_info[[class(object)[1]]]
  engine <- info$engine
  settings <- object$settings
  qr <- isTRUE(object$qr)
  fixed <- identical(settings$speed_var, "fixed")
  cov_names <- colnames(object$post$settings$X %||% object$data$X) %||% character()
  names <- par_names(engine, object$data$items, cov_names)
  priors <- stats::setNames(rep("", length(names)), names)
  set <- function(pattern, value) {
    priors[grepl(pattern, names)] <<- value
  }
  text <- prior_text(fit_priors(object), object$data$log_t)

  set("^a\\[", if (identical(settings$itemtype, "1pl")) "fixed 1" else text$a)
  set("^b\\[", text$b)
  if (info$rt) {
    set("^lambda\\[", text$lambda)
    set("^sigma2t\\[", if (qr && engine == "cross") paste(text$sigma2t, "(ALD scale)") else text$sigma2t)
  }
  set("^(beta|beta_ability|beta_speed)\\[", text$beta)
  if (engine == "rtirt") {
    correlation <- sprintf("c / sqrt(c^2 + v): c ~ %s, v ~ %s", text$cov, text$var_speed)
    set("^cor_ability_speed$", if (fixed) "U(-1, 1)" else correlation)
    set("^var_speed$", if (fixed) "fixed 1" else "c^2 + v")
  }
  if (engine == "latent") {
    set("^b_ability$", text$beta)
    set("^var_speed$", if (qr) paste(text$var_speed, "(ALD scale)") else text$var_speed)
  }
  if (engine == "cross") {
    set("^rho\\[", text$rho)
    set("^var_speed$", if (fixed) "fixed 1" else text$var_speed)
  }
  priors
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
  if (has_em(object)) {
    return(summary_em(object, digits))
  }
  need_post(object)
  table <- estimates(object)
  priors <- prior_strings(object)
  table$prior <- unname(priors[table$parameter])
  table$prior[is.na(table$prior)] <- ""
  ald <- isTRUE(object$qr) && class(object)[1] == "rtirt_cross"
  partable <- rtirt_partable(table, TRUE, ald, "posterior mean, SD and 95% interval; prior of each parameter")
  wide <- rtirt_item_wide(partable, object$data$items)
  person_cov <- object$post$mean$sigma_p

  tables <- list(
    header = rtirt_header(object),
    fit = rtable(gibbs_fit_table(object, table), gibbs_fit_title(object)),
    items = rtable(wide$est, "item parameters (posterior means)"),
    se = rtable(wide$se, "posterior SDs of the item parameters"),
    parameters = partable,
    covariance = person_cov_table(person_cov, "person covariance (posterior mean)"),
    reliability = rtable(as.data.frame(as.list(unclass(reliability(object)))),
                         "reliability of the EAP scores: var(EAP) / (var(EAP) + mean PSD^2)"))
  structure(tables, class = c("summary.rtirt", "rtirt_summary"))
}

gibbs_fit_table <- function(object, table) {
  max_rhat <- suppressWarnings(max(table$rhat, na.rm = TRUE))
  cbind(N = nrow(object$data$Y), max_rhat = max_rhat, fit_indices(object))
}

gibbs_fit_title <- function(object) {
  ald_note <- if (isTRUE(object$qr)) "; ALD working likelihood: not comparable with mean models or across quantile levels" else ""
  paste0("marginal likelihood (persons integrated out): compare models by WAIC or LOOIC (smaller is better); ",
         "DIC_complete conditions on the persons (as the Julia package)", ald_note)
}

person_cov_table <- function(person_cov, title) {
  if (is.null(person_cov)) {
    return(NULL)
  }
  rtable(data.frame(variable = rownames(person_cov), person_cov, check.names = FALSE), title)
}

#' Person scores
#'
#' Posterior means (EAP) and posterior SDs of ability and speed, one row per person.
#' @param object A sampled model.
#' @param ... Unused.
#' @return A tibble (class [rtirt_tbl]: printed with 3 decimals, stored in full precision) with
#'   `id`, `ability`, `ability_psd` and (RT models) `speed`, `speed_psd`.
#' @export
scores <- function(object, ...) {
  UseMethod("scores")
}

#' @rdname scores
#' @export
scores.rtirt <- function(object, ...) {
  person <- person_estimates(object)
  table <- data.frame(id = object$data$id, ability = person$theta$mean, ability_psd = person$theta$sd)
  if (!is.null(person$zeta)) {
    table$speed <- person$zeta$mean
    table$speed_psd <- person$zeta$sd
  }
  as_rtirt_tbl(tibble::as_tibble(table))
}

# posterior means and SDs of the persons: ECM (EAP) or Gibbs
person_estimates <- function(object) {
  if (has_em(object)) {
    return(object$ecm$person)
  }
  need_post(object)
  object$post$person
}

#' Empirical reliability
#'
#' var(EAP) / (var(EAP) + mean PSD^2) of ability and speed.
#' @param object A sampled model.
#' @param ... Unused.
#' @export
reliability <- function(object, ...) {
  UseMethod("reliability")
}

#' @rdname reliability
#' @export
reliability.rtirt <- function(object, ...) {
  person <- person_estimates(object)
  eap_reliability <- function(p) {
    eap_var <- stats::var(p$mean)
    eap_var / (eap_var + mean(p$sd^2))
  }
  speed <- if (!is.null(person$zeta)) c(speed = eap_reliability(person$zeta))
  as_rtirt_values(c(ability = eap_reliability(person$theta), speed))
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
fit_indices <- function(object, ...) {
  UseMethod("fit_indices")
}

#' @rdname fit_indices
#' @export
fit_indices.rtirt <- function(object, ...) {
  if (has_em(object)) {
    fit <- object$ecm
    method <- if (fit$map) "ECM (MAP items)" else "ECM"
    return(as_rtirt_tbl(tibble::tibble(method = method, logLik = fit$logLik, npar = fit$npar, N = fit$nobs,
                                       AIC = fit$AIC, BIC = fit$BIC)))
  }
  need_post(object)
  marginal <- object$post$fit
  if (is.null(marginal)) {
    message("computing the marginal log-likelihood ...")
    marginal <- marginal_fit(object)
  }
  as_rtirt_tbl(tibble::as_tibble(cbind(marginal$table, DIC_complete = dic(object)$dic)))
}

# marginal-likelihood criteria (stored by gibbs(fit_indices = TRUE))
marginal_fit <- function(object) {
  loglik <- unclass(marginal_loglik(object))
  draws <- post_draws(object, Inf)
  cov_names <- colnames(object$post$settings$X) %||% character()
  post_mean <- colMeans(draws)
  cross_ald <- isTRUE(object$qr) && .model_info[[class(object)[1]]]$engine == "cross"
  nodes <- if (cross_ald) 61 else 41
  at_mean <- marg_ll_draw(object, draw_pars(post_mean, object$data$items, cov_names), Q = nodes)
  d_hat <- -2 * sum(at_mean)
  d_bar <- -2 * mean(rowSums(loglik))
  person_max <- apply(loglik, 2, max)
  lppd <- sum(person_max + log(colMeans(exp(sweep(loglik, 2, person_max)))))
  p_waic <- sum(apply(loglik, 2, stats::var))
  table <- data.frame(Deviance = d_bar, pD = d_bar - d_hat, DIC = 2 * d_bar - d_hat, WAIC = -2 * (lppd - p_waic),
                      p_waic = p_waic, LOOIC = NA_real_, SE_LOOIC = NA_real_)
  loo_fit <- NULL
  if (requireNamespace("loo", quietly = TRUE)) {
    loo_fit <- suppressWarnings(loo::loo(loglik))
    table$LOOIC <- loo_fit$estimates["looic", "Estimate"]
    table$SE_LOOIC <- loo_fit$estimates["looic", "SE"]
  }
  list(table = table, loglik = loglik, loo = loo_fit)
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
dic <- function(object, ...) {
  UseMethod("dic")
}

#' @rdname dic
#' @export
dic.rtirt <- function(object, ...) {
  need_post(object)
  after_burnin <- (object$post$burnin + 1):object$cond$n_iter
  loglik <- object$post$loglik[after_burnin, , drop = FALSE]
  d_bar <- -2 * mean(loglik)
  d_hat <- -2 * complete_loglik(object)
  as_rtirt_tbl(data.frame(d_bar = d_bar, d_hat = d_hat, p_d = d_bar - d_hat, dic = d_bar + (d_bar - d_hat)))
}

# complete-data log-likelihood at the posterior means (mirrors the loglik in the C++ samplers)
complete_loglik <- function(object) {
  m <- object$post$mean
  data <- object$data
  X <- object$post$settings$X
  engine <- .model_info[[class(object)[1]]]$engine
  n_person <- nrow(data$Y)
  n_item <- ncol(data$Y)
  by_item <- function(v) matrix(v, n_person, n_item, byrow = TRUE)
  normal_ll <- function(x, mu, v) sum(stats::dnorm(x, mu, sqrt(v), log = TRUE))
  q <- object$post$settings$q
  ll <- sum(calc_loglik_nodes(data$Y, m$a, m$b, matrix(m$theta)))           # theta as one person-specific node

  if (engine == "mlirt") {
    return(ll + normal_ll(m$theta, drop(X %*% m$beta), 1))
  }
  if (engine == "cross") {
    mu <- by_item(m$lambda) - m$zeta - outer(m$theta, m$rho)
    v <- by_item(m$sigma2t)
    if (object$qr) {
      k <- ald_k(q)
      mu <- mu + k[1] * m$nu
      v <- k[2] * v * m$nu
    }
    return(ll + normal_ll(data$log_t, mu, v) + normal_ll(m$theta, 0, 1) + normal_ll(m$zeta, 0, m$sigma_p[2, 2]))
  }
  ll <- ll + normal_ll(data$log_t, by_item(m$lambda) - m$zeta, by_item(m$sigma2t))
  if (engine == "latent") {
    speed_mean <- drop(X %*% m$beta) + m$b_theta * m$theta
    speed_var <- m$sigma_p[2, 2]
    if (object$qr) {
      k <- ald_k(q)
      return(ll + normal_ll(m$theta, 0, 1) + normal_ll(m$zeta, speed_mean + k[1] * m$nu, k[2] * speed_var * m$nu))
    }
    return(ll + normal_ll(m$theta, 0, 1) + normal_ll(m$zeta, speed_mean, speed_var))
  }
  person_cov <- m$sigma_p
  ability_mean <- drop(X %*% m$beta[, 1])
  speed_mean <- drop(X %*% m$beta[, 2])
  cov_ability_speed <- person_cov[1, 2]
  speed_resid_var <- person_cov[2, 2] - cov_ability_speed^2
  speed_given_ability <- speed_mean + cov_ability_speed * (m$theta - ability_mean)
  ll + normal_ll(m$theta, ability_mean, 1) + normal_ll(m$zeta, speed_given_ability, speed_resid_var)
}

#' Convergence check
#'
#' The share of item and structural parameters with effective sample size above `ess` and
#' R-hat below `rhat` (over independent chains), and the worst parameters. With
#' `diagnostics = "rank"` (see [estimates()]) both the bulk and the tail ESS must exceed `ess`,
#' and the default R-hat threshold is 1.01 (Vehtari et al., 2021) instead of 1.05.
#' Parameters that do not vary (fixed ones) are left out. With one chain R-hat cannot be
#' computed and only the ESS is checked. `print()` of a [gibbs()] fit warns with the same
#' criteria (default thresholds) when a parameter fails them.
#' Julia: `checkConvergence()`.
#'
#' For several quantile levels (an `rtirt_qset` from `gibbs(..., quantile = c(.25, .75))`)
#' each level is checked; the model of one level is `fits$models[[k]]`.
#' @param object A sampled model, an `rtirt_qset`, or an [ecm()] fit (whether ECM converged).
#' @param ess,rhat Thresholds (`rhat` defaults to 1.05, or 1.01 with `diagnostics = "rank"`).
#' @param diagnostics `"classic"` or `"rank"`, as in [estimates()].
#' @param x An `rtirt_convergence` object.
#' @param ... Unused.
#' @return An object of class `rtirt_convergence`, printed by its `print()` method: `share`
#'   (share of parameters that pass), `worst` (the five parameters with the largest R-hat, or
#'   the smallest ESS with one chain), `n` (parameters checked), `n_chain` and the thresholds;
#'   for an ECM fit `converged` and `iterations`; for an `rtirt_qset` one such object per level
#'   (`levels`) with `quantile`.
#' @examples
#' \donttest{
#' cond <- set_cond(n_subj = 200, n_item = 6)
#' dat <- sim_data(cond, sim_para(cond, "null"), "null")
#' fit <- gibbs(dat, model = "null", n_iter = 1000, n_chain = 2, seed = 1)
#' check <- convergence(fit)   # printed only when not assigned
#' check
#' check$share
#' }
#' @export
convergence <- function(object, ...) {
  UseMethod("convergence")
}

# ESS above `ess` and R-hat below `rhat` for each parameter that varies (R-hat is NA with one chain)
convergence_check <- function(object, ess = 400, rhat = NULL, diagnostics = "classic") {
  rank <- diagnostics == "rank"
  rhat <- rhat %||% if (rank) 1.01 else 1.05
  table <- estimates(object, diagnostics = diagnostics)
  varies <- !is.na(table$ess) & (if ("sd" %in% names(table)) table$sd > 0 else TRUE)
  table <- table[varies, , drop = FALSE]
  tail_ok <- if (!"ess_tail" %in% names(table)) TRUE else table$ess_tail > ess
  ok <- table$ess > ess & (is.na(table$rhat) | table$rhat < rhat) & tail_ok
  list(table = table, ok = ok, share = if (length(ok)) mean(ok) else 1, n = nrow(table), ess = ess, rhat = rhat,
       rank = rank, n_chain = coda::nchain(object$post$draws))
}

#' @rdname convergence
#' @export
convergence.rtirt <- function(object, ess = 400, rhat = NULL, diagnostics = c("classic", "rank"), ...) {
  if (has_em(object)) {                                                # an ECM fit: the iterations, not chains
    fit <- object$ecm
    return(structure(list(converged = fit$converged, iterations = fit$iterations), class = "rtirt_convergence"))
  }
  need_post(object)
  diagnostics <- match.arg(diagnostics)
  check <- convergence_check(object, ess, rhat, diagnostics)
  table <- check$table
  columns <- intersect(c("parameter", "rhat", "ess", "ess_tail"), names(table))
  order_by <- if (check$n_chain > 1) order(-table$rhat, table$ess) else order(table$ess)
  worst <- table[order_by, columns][seq_len(min(5, nrow(table))), ]
  result <- list(share = check$share, worst = worst, n = check$n, n_chain = check$n_chain, ess = check$ess,
                 rhat = check$rhat, rank = check$rank)
  structure(result, class = "rtirt_convergence")
}

#' @rdname convergence
#' @export
convergence.rtirt_qset <- function(object, ...) {
  if (!qset_sampled(object)) {
    stop("run gibbs() first", call. = FALSE)
  }
  levels <- lapply(object$models, convergence, ...)
  structure(list(levels = levels, quantile = object$quantile), class = "rtirt_convergence")
}

#' @rdname convergence
#' @export
print.rtirt_convergence <- function(x, ...) {
  if (!is.null(x$levels)) {                                            # several quantile levels
    for (k in seq_along(x$levels)) {
      cat(sprintf("q = %s: ", format(x$quantile[k])))
      print(x$levels[[k]])
      if (k < length(x$levels)) cat("\n")
    }
    return(invisible(x))
  }
  if (!is.null(x$converged)) {
    status <- if (x$converged) "converged" else "NOT converged"
    cat(sprintf("ECM: %s in %d iterations (relative change of the objective below tol)\n", status, x$iterations))
    return(invisible(x))
  }
  ess_note <- if (x$rank) " (bulk and tail)" else ""
  rhat_note <- if (x$rank) " (rank-normalized)" else ""
  cat(sprintf("%.1f%% of %d parameters with ESS%s > %d and R-hat%s < %s\n", 100 * x$share, x$n,
              ess_note, x$ess, rhat_note, format(x$rhat)))
  if (x$n_chain < 2) {
    cat("  1 chain: R-hat cannot be computed (it needs n_chain >= 2), so only the ESS is checked\n")
  }
  print_table(x$worst)
  invisible(x)
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
  if (!has_em(object)) {
    need_post(object)
  }
  if (is.null(object$true_para)) {
    stop("no true_para in the model")
  }
  label <- par
  if (par == "rho" && inherits(object, "rtirt_latent")) {
    par <- "b_theta"
  }
  par <- switch(par, ability = "theta", b_ability = "b_theta", par)       # internal (Julia) names
  source <- if (has_em(object)) {
    c(object$ecm$par, list(theta = object$ecm$person$theta$mean, zeta = object$ecm$person$zeta$mean))
  } else {
    object$post$mean
  }
  if (is.null(source[[par]])) {
    stop(sprintf("'%s' is not among the estimates of this fit", label), call. = FALSE)
  }
  est <- as.numeric(source[[par]])
  true <- as.numeric(object$true_para[[par]])
  if (length(est) != length(true)) {
    stop(sprintf("'%s': %d estimates vs %d true values", par, length(est), length(true)))
  }
  can_correlate <- length(est) > 2 && stats::sd(est) > 0 && stats::sd(true) > 0
  corr <- if (can_correlate) stats::cor(est, true) else NA_real_
  as_rtirt_tbl(tibble::tibble(par = label, n = length(est), bias = get_bias(est, true), rmse = get_rmse(est, true), corr = corr))
}

#' @importFrom coda as.mcmc.list
#' @export
as.mcmc.list.rtirt <- function(x, ...) {
  need_post(x)
  post_mcmc(x)
}

#' Plots of a sampled model
#'
#' `type = "items"`: posterior means and 95% intervals of the item parameters (ECM fits: estimates
#' and Wald 95% intervals), one panel
#' per parameter type (a, b, lambda, sigma2t, rho). `type = "structural"`: the same for the
#' structural parameters (regression coefficients, `cor_ability_speed`, `b_ability`,
#' `var_speed`), or those matching `pars`. `type = "trace"`: trace plots (all
#' iterations; burn-in shaded) of the structural parameters or of those matching `pars`.
#' For plots of your own (e.g. with ggplot2), `estimates(fit)` (or `as.data.frame(fit)`) has
#' the estimates and intervals as a data frame.
#' @param x A sampled model.
#' @param type `"items"`, `"structural"` or `"trace"`.
#' @param pars Regular expression selecting parameters (`"structural"` and `"trace"`; at most 12
#'   are drawn in `"trace"`).
#' @param ... Passed to [graphics::plot()].
#' @export
plot.rtirt <- function(x, type = c("items", "structural", "trace"), pars = NULL, ...) {
  type <- match.arg(type)
  if (type == "trace" || !has_em(x)) {                 # items: posterior or Wald 95% intervals (ECM)
    need_post(x)
  }
  table <- estimates(x)
  if (type == "items") {
    plot_items(x, table, ...)
    return(invisible(table))
  }
  if (type == "structural") {
    pars <- pars %||% "^(beta|b_ability|cor_ability_speed|var_speed)"
    rows <- table[grepl(pars, table$parameter), ]
    if (!nrow(rows)) {
      stop("no parameters match '", pars, "'", call. = FALSE)
    }
    plot_intervals(rows, ...)
    return(invisible(rows))
  }
  pars <- pars %||% "^(rho|beta|b_ability|cor_ability_speed|var_speed)"
  selected <- table$parameter[grepl(pars, table$parameter)]
  if (!length(selected)) {
    stop("no parameters match '", pars, "'")
  }
  selected <- utils::head(selected, 12)
  plot_traces(x, selected, ...)
  invisible(table[table$parameter %in% selected, ])
}

# one panel per item block with all items: estimates and 95% intervals
plot_items <- function(x, table, ...) {
  indexed <- table$parameter[grepl("\\[", table$parameter)]
  blocks <- intersect(item_blocks, unique(sub("\\[.*", "", indexed)))
  n_items_in <- function(block) sum(startsWith(table$parameter, paste0(block, "[")))
  blocks <- blocks[vapply(blocks, function(block) n_items_in(block) == ncol(x$data$Y), TRUE)]
  left_margin <- if (length(blocks)) 7 else 4
  old_par <- graphics::par(mfrow = c(1, length(blocks)), mar = c(4, left_margin, 2.5, 1))
  on.exit(graphics::par(old_par))
  for (block in blocks) {
    rows <- table[startsWith(table$parameter, paste0(block, "[")), ]
    n_item <- nrow(rows)
    positions <- n_item:1
    graphics::plot(rows$est, positions, xlim = range(rows$q025, rows$q975, if (block == "rho") 0), yaxt = "n", pch = 19,
                   xlab = block, ylab = "", main = block, ...)
    graphics::segments(rows$q025, positions, rows$q975, positions)
    if (block == "rho") {
      graphics::abline(v = 0, lty = 3)
    }
    if (block == blocks[1]) {
      graphics::axis(2, positions, x$data$items, las = 1, cex.axis = 0.8)
    }
  }
}

# estimates and 95% intervals of a few parameters, one row each, 0 marked
plot_intervals <- function(rows, ...) {
  n_par <- nrow(rows)
  positions <- n_par:1
  old_par <- graphics::par(mar = c(4, max(4, 0.5 * max(nchar(rows$parameter))), 2.5, 1))
  on.exit(graphics::par(old_par))
  graphics::plot(rows$est, positions, xlim = range(rows$q025, rows$q975, 0), ylim = c(0.5, n_par + 0.5), yaxt = "n",
                 pch = 19, xlab = "estimate and 95% interval", ylab = "", main = "structural parameters", ...)
  graphics::segments(rows$q025, positions, rows$q975, positions)
  graphics::abline(v = 0, lty = 3)
  graphics::axis(2, positions, rows$parameter, las = 1, cex.axis = 0.8)
}

# trace plots of the chains, burn-in shaded
plot_traces <- function(x, selected, ...) {
  n_col <- ceiling(sqrt(length(selected)))
  old_par <- graphics::par(mfrow = c(ceiling(length(selected) / n_col), n_col), mar = c(4, 3, 2, 1))
  on.exit(graphics::par(old_par))
  for (p in selected) {
    chains <- sapply(x$post$draws, function(m) m[, p])
    iterations <- seq_len(nrow(chains)) * x$cond$n_thin
    graphics::matplot(iterations, chains, type = "l", lty = 1, col = seq_len(ncol(chains)), xlab = "iteration", ylab = "",
                      main = p, ...)
    usr <- graphics::par("usr")
    graphics::rect(0, usr[3], x$post$burnin, usr[4], col = grDevices::adjustcolor("grey", 0.3), border = NA)
  }
}

# ---- several quantile levels ---------------------------------------------------------------

#' Several quantile levels
#'
#' `fits <- gibbs(data, model = "cross", quantile = c(.1, .5, .9))` (or `model = "latent"`)
#' samples one ALD model per level and returns them as an `rtirt_qset`. `coef()` gives a parameters x
#' levels matrix, `estimates()` stacks the tables with a column `q`, and `plot()` draws the
#' quantile coefficients (cross-relations `rho`, or the regression coefficients of
#' `"latent"`) against the level with 95% intervals. [convergence()] checks every level.
#' `summary()` returns an object that prints as `print()`.
#' The fit of one level is `fits$models[[k]]` (an ordinary fit: `summary()`, `estimates()`,
#' `plot(type = "trace")`).
#' @param x,object An `rtirt_qset`.
#' @param pars Regular expression of the parameters to plot (default: `rho` for
#'   `"cross"`, `beta` and `b_ability` for `"latent"`).
#' @param ... Unused, or passed to [graphics::plot()].
#' @name rtirt_qset
NULL

qset_sampled <- function(x) {
  all(vapply(x$models, function(m) !is.null(m$post), TRUE))
}

#' @rdname rtirt_qset
#' @export
print.rtirt_qset <- function(x, ...) {
  cat(sprintf("%s (%s) at quantiles %s\n", model_call_text(x$type), x$type, paste(format(x$quantile), collapse = ", ")))
  if (!qset_sampled(x)) {
    cat("  not sampled yet: fits <- gibbs(models)\n")
    return(invisible(x))
  }
  coefs <- unclass(stats::coef(x))
  key_pattern <- if (x$type == "rtirt_cross") "^rho\\[" else "^(beta|b_ability)"
  key <- grepl(key_pattern, rownames(coefs))
  print_table(data.frame(parameter = rownames(coefs)[key], coefs[key, , drop = FALSE], check.names = FALSE))
  for (k in seq_along(x$models)) {
    convergence_warn(x$models[[k]], sprintf("q = %s", format(x$quantile[k])))
  }
  invisible(x)
}

#' @rdname rtirt_qset
#' @export
coef.rtirt_qset <- function(object, ...) {
  coefs <- sapply(object$models, stats::coef)
  colnames(coefs) <- sprintf("q=%s", format(object$quantile))
  as_rtirt_values(coefs)
}

#' @rdname rtirt_qset
#' @export
estimates.rtirt_qset <- function(object, ...) {
  one_level <- function(model, q) cbind(q = q, estimates_rtirt(model, NULL, "auto", "classic"))
  tables <- Map(one_level, object$models, object$quantile)
  as_rtirt_tbl(tibble::as_tibble(do.call(rbind, tables)))
}

#' @rdname rtirt_qset
#' @export
summary.rtirt_qset <- function(object, ...) {
  structure(list(fits = object), class = "summary.rtirt_qset")
}

#' @rdname rtirt_qset
#' @export
print.summary.rtirt_qset <- function(x, ...) {
  print(x$fits, ...)
  invisible(x)
}

#' @rdname rtirt_qset
#' @export
plot.rtirt_qset <- function(x, pars = NULL, ...) {
  if (!qset_sampled(x)) {
    stop("run gibbs() first", call. = FALSE)
  }
  pars <- pars %||% if (x$type == "rtirt_cross") "^rho\\[" else "^(beta|b_ability)"
  table <- estimates(x)
  table <- table[grepl(pars, table$parameter), ]
  key <- unique(table$parameter)
  if (!length(key)) {
    stop("no parameters match '", pars, "'")
  }
  n_col <- ceiling(sqrt(length(key)))
  old_par <- graphics::par(mfrow = c(ceiling(length(key) / n_col), n_col), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old_par))
  for (k in key) {
    rows <- table[table$parameter == k, ]
    rows <- rows[order(rows$q), ]
    graphics::plot(rows$q, rows$est, type = "n", ylim = range(rows$q025, rows$q975, 0), xlab = "quantile", ylab = "",
                   main = k, ...)
    graphics::polygon(c(rows$q, rev(rows$q)), c(rows$q025, rev(rows$q975)), col = "grey85", border = NA)
    graphics::lines(rows$q, rows$est, type = "b", pch = 19)
    graphics::abline(h = 0, lty = 3)
  }
  invisible(table)
}

#' @export
gibbs.default <- function(data, ...) {
  stop("data must come from input_data(), e.g. gibbs(input_data(...), model = \"cross\")", call. = FALSE)
}
