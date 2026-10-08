# Tests about the response times: RT person fit and conditional independence of times and responses ----

#' Person fit of the response times
#'
#' The person-fit statistic of Marianti, Fox, Avetisyan, Veldkamp & Tijmstra (2014) for the lognormal
#' RT model: \eqn{l^t_j = \sum_i (\log t_{ij} - \lambda_i + \zeta_j [+ \rho_i\theta_j])^2 / \sigma_i^2},
#' chi-square with K degrees of freedom when the model holds. The samplers of [gibbs()] evaluate it in
#' every iteration after burn-in, so the result accounts for the uncertainty of all parameters:
#' * `lt`: posterior mean of \eqn{l^t_j};
#' * `p_flag`: posterior probability that \eqn{l^t_j} exceeds the 95% quantile of chi-square(K), the
#'   probability that the pattern is flagged; Marianti et al. flag a pattern when it is at least .95.
#'
#' Large values point to aberrant timing: rapid guessing, item preknowledge, running out of time,
#' or a speed that changes during the test. Normal log times only (not `rtirt_cross(quantile = )`).
#' @param object A model sampled by [gibbs()].
#' @param flag Posterior probability above which a pattern is flagged.
#' @param ... Unused.
#' @return A data frame with `id`, `lt`, `df`, `p_flag` and `flagged`; attribute `share` (share of
#'   flagged test takers) and `share_mean` (share with `lt` above the 95% quantile).
#' @references Marianti, S., Fox, J.-P., Avetisyan, M., Veldkamp, B. P., & Tijmstra, J. (2014).
#'   Testing for aberrant behavior in response time modeling. *Journal of Educational and Behavioral
#'   Statistics, 39*, 426-451.
#'
#'   Fox, J.-P., & Marianti, S. (2017). Person-fit statistics for joint models for accuracy and
#'   speed. *Journal of Educational Measurement, 54*, 243-262.
#' @examples
#' \donttest{
#' cond <- set_cond(n_subj = 300, n_item = 10)
#' d <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' f <- gibbs(rtirt_null(d), n_iter = 2000, n_chain = 2, seed = 1, verbose = FALSE,
#'            fit_indices = FALSE)
#' pf <- person_fit(f); attr(pf, "share")        # close to 0 under the model
#' }
#' @export
person_fit <- function(object, flag = 0.95, ...) {
  need_post(object)
  if (!.model_info[[class(object)[1]]]$rt) stop("person_fit() needs a model with response times", call. = FALSE)
  if (identical(.model_info[[class(object)[1]]]$engine, "cross") && isTRUE(object$qr))
    stop("the chi-square person fit needs normal log times; rtirt_cross(quantile = ) has ALD times", call. = FALSE)
  pf <- object$post$person$rt_fit
  if (is.null(pf)) stop("refit with this version of birtRcpp: the fit has no person-fit statistics", call. = FALSE)
  K <- ncol(object$data$Y)
  out <- data.frame(id = object$data$id %||% seq_along(pf$lt), lt = pf$lt, df = K, p_flag = pf$exceed, flagged = pf$exceed >= flag)
  structure(out, share = mean(out$flagged), share_mean = mean(out$lt > stats::qchisq(0.95, K)))
}

#' Test of conditional independence between responses and response times
#'
#' The Lagrange multiplier (score) test of van der Linden & Glas (2010): under conditional
#' independence the distribution of the log time on item i does not depend on whether the response is
#' correct, given ability and speed. The alternative shifts the log times of correct responses,
#' \eqn{\log t_{ij} = \lambda_i + \delta_i y_{ij} - \zeta_j [- \rho_i\theta_j] + e_{ij}}, and the
#' test is of \eqn{\delta_i = 0} from the maximum likelihood fit of the model without the shift
#' ([ecm()]). The score of \eqn{\delta_i} for person j is \eqn{y_{ij} E[\log t_{ij} - \lambda_i +
#' \zeta_j (+ \rho_i\theta_j) | \text{data}_j] / \sigma_i^2} (Fisher identity). The information
#' is the observed information of the model extended by the shifts (central differences of the
#' analytic score, as for the standard errors of [ecm()]) with the parameters of the model projected
#' out, so the statistic is chi-square(1) per item, and chi-square(K) for all items together. (The
#' outer product of the casewise scores is simpler but rejects too often: about .065 per item and
#' .10 overall at the 5% level with 500 persons and 10 items.)
#' A positive `delta` means that correct responses took longer than incorrect ones.
#' @param object A model fitted by [ecm()] (normal models with response times).
#' @param ... Unused.
#' @return A data frame of class `rtirt_ci_test` per item: `item`, `delta` (one-step estimate of the
#'   shift), `statistic`, `p.value`, `p.holm`; attribute `overall` (joint test of all items).
#' @references van der Linden, W. J., & Glas, C. A. W. (2010). Statistical tests of conditional
#'   independence between responses and/or response times on test items. *Psychometrika, 75*,
#'   120-139.
#' @examples
#' cond <- set_cond(n_subj = 500, n_item = 8)
#' d <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' ci_test(ecm(rtirt_null(d), se = FALSE))
#' @export
ci_test <- function(object, ...) {
  if (is.null(object$ecm) || !.model_info[[class(object)[1]]]$rt) stop("ci_test() needs an ecm() fit of a model with response times", call. = FALSE)
  if (isTRUE(object$ecm$map)) stop("ci_test() uses the maximum likelihood fit: ecm(model), not ecm(model, prior = \"map\")", call. = FALSE)
  E <- em_setup(object, object$ecm$nodes, FALSE); x0 <- object$ecm$x; np <- length(x0); K <- E$K; N <- E$N
  # score of the model extended by the shifts delta (log t - delta y): nuisance part (Fisher identity, em_score())
  # and the part of delta, sum_i y_ij E[log t - delta_j y - lambda + zeta + rho theta | data] / sigma2_j
  score <- function(z) {
    E2 <- E; E2$L <- E$L - E$Y * matrix(z[np + seq_len(K)], N, K, byrow = TRUE)
    p <- em_unpack(z[seq_len(np)], E2); mo <- em_moments(em_estep(p, E2))
    rho <- if (E$eng == "cross") p$rho else rep(0, K)
    R <- sweep(E2$L, 2, p$lambda) + mo$Ez + outer(mo$Et, rho)
    c(em_score(z[seq_len(np)], E2), colSums(E$Y * R) / p$sigma2t)
  }
  z0 <- c(x0, rep(0, K)); g <- score(z0)[np + seq_len(K)]
  # observed information of the extended model: central differences of the analytic score (as ecm())
  H <- vapply(seq_along(z0), function(k) { h <- 1e-4 * max(1, abs(z0[k])); e <- replace(numeric(length(z0)), k, h)
    (score(z0 + e) - score(z0 - e)) / (2 * h) }, numeric(length(z0)))
  I <- -(H + t(H)) / 2; th <- seq_len(np); dl <- np + seq_len(K)
  Ieff <- I[dl, dl] - I[dl, th] %*% solve(I[th, th], I[th, dl])         # nuisance parameters projected out
  stat <- g^2 / diag(Ieff)
  out <- data.frame(item = E$items, delta = g / diag(Ieff), statistic = stat, p.value = stats::pchisq(stat, 1, lower.tail = FALSE))
  out$p.holm <- stats::p.adjust(out$p.value, "holm")
  ov <- drop(t(g) %*% solve(Ieff, g))
  structure(out, class = c("rtirt_ci_test", "data.frame"), overall = c(statistic = ov, df = K, p.value = stats::pchisq(ov, K, lower.tail = FALSE)))
}

#' @export
print.rtirt_ci_test <- function(x, digits = 3, ...) {
  o <- attr(x, "overall")
  cat("Conditional independence of responses and response times (van der Linden & Glas, 2010)\n")
  cat(sprintf("  all items: LM = %.2f, df = %d, p = %.3g\n\n", o[["statistic"]], as.integer(o[["df"]]), o[["p.value"]]))
  print.data.frame(x, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Tests of item preknowledge
#'
#' Statistics of Sinharay (2020) for test takers who may have known some items in advance (a set of
#' possibly compromised items is given), computed from the item parameters of a fit of a model with
#' response times (`ecm()` estimates, or the posterior means of `gibbs()`), with the lognormal RT model
#' \eqn{\log t_{ij} \sim N(\lambda_i - \tau_j, \sigma_i^2)} and the 2PL for the scores:
#' * `chi_pf`: RT person fit \eqn{\sum_i (\log t_{ij} - \lambda_i + \hat\tau_j)^2/\sigma_i^2} with the
#'   speed estimated from the times, chi-square with K - 1 df (it does not use the compromised set);
#' * `Lambda_s`: signed likelihood ratio (a z statistic) comparing the speed on the compromised items with
#'   the speed on the other items, \eqn{(\hat\tau_c - \hat\tau_{\bar c}) / \sqrt{1/W_c + 1/W_{\bar c}}},
#'   \eqn{W = \sum 1/\sigma_i^2}; large values: faster on the compromised items;
#' * `L_s`: signed likelihood ratio of the item scores comparing the ability (maximum likelihood) on the
#'   two sets; large values: better on the compromised items.
#' The p-values are one-sided for `Lambda_s` and `L_s` (standard normal) and upper-tail for `chi_pf`.
#' Models with the speed entering as \eqn{\lambda_i - \zeta_j} ([rtirt_null()], [rtirt_latreg()],
#' [rtirt_latent()] without quantile RT residuals).
#' @param object A model fitted by [ecm()] or [gibbs()].
#' @param compromised Items in the compromised set (names or a logical / integer index).
#' @param ... Unused.
#' @return A data frame with one row per test taker: the three statistics and their p-values.
#' @references Sinharay, S. (2020). Detection of item preknowledge using response times. *Applied
#'   Psychological Measurement, 44*, 376-392.
#' @examples
#' cond <- set_cond(n_subj = 400, n_item = 20)
#' d <- sim_data(cond, sim_para(cond, "null"), "null")
#' pk <- preknowledge_test(ecm(rtirt_null(d), se = FALSE), compromised = 1:5)
#' colMeans(pk[c("p_chi_pf", "p_Lambda_s", "p_L_s")] < 0.01)    # about .01 under the model
#' @export
preknowledge_test <- function(object, compromised, ...) {
  eng <- .model_info[[class(object)[1]]]$engine
  if (!eng %in% c("rtirt", "latent") || !.model_info[[class(object)[1]]]$rt) stop("preknowledge_test() needs rtirt_null(), rtirt_latreg() or rtirt_latent()", call. = FALSE)
  par <- if (!is.null(object$post)) object$post$mean else { p <- object$ecm$par; list(a = p$a, b = p$b, lambda = p$lambda, sigma2t = p$sigma2t) }
  if (is.null(par$a)) stop("no item parameters: fit the model with ecm() or gibbs()", call. = FALSE)
  D <- object$data; Y <- D$Y; L <- D$log_t; K <- ncol(Y); items <- D$items
  cset <- if (is.logical(compromised)) compromised else if (is.character(compromised)) items %in% compromised else seq_len(K) %in% compromised
  if (!any(cset) || all(cset)) stop("compromised must select some but not all items", call. = FALSE)
  w <- 1 / par$sigma2t; lam <- par$lambda; a <- par$a; b <- par$b
  R <- sweep(-L, 2, lam, "+")                                            # lambda - log t: tau + noise
  tau <- function(s) drop(R[, s, drop = FALSE] %*% w[s]) / sum(w[s])
  th <- tau(rep(TRUE, K))
  chi <- rowSums(sweep((L - matrix(lam, nrow(L), K, byrow = TRUE) + th)^2, 2, w, "*"))
  Lam <- (tau(cset) - tau(!cset)) / sqrt(1 / sum(w[cset]) + 1 / sum(w[!cset]))
  ll <- function(t, y, s) sum(y * stats::plogis(a[s] * (t - b[s]), log.p = TRUE) + (1 - y) * stats::plogis(-a[s] * (t - b[s]), log.p = TRUE))
  mle <- function(y, s) stats::optimize(ll, c(-8, 8), y = y, s = s, maximum = TRUE)
  Ls <- vapply(seq_len(nrow(Y)), function(i) {
    f <- mle(Y[i, ], rep(TRUE, K)); fc <- mle(Y[i, cset], cset); fu <- mle(Y[i, !cset], !cset)
    G <- max(2 * (fc$objective + fu$objective - f$objective), 0)
    sign(fc$maximum - fu$maximum) * sqrt(G)
  }, 0)
  data.frame(id = D$id %||% seq_len(nrow(Y)), chi_pf = chi, p_chi_pf = stats::pchisq(chi, K - 1, lower.tail = FALSE),
             Lambda_s = Lam, p_Lambda_s = stats::pnorm(Lam, lower.tail = FALSE), L_s = Ls, p_L_s = stats::pnorm(Ls, lower.tail = FALSE))
}
