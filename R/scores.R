# Casewise scores of the marginal likelihood and score-based tests of parameter invariance ------
#
# The score of person i is the derivative of log p(y_i, log t_i | parameters), with ability and speed
# integrated out, at the ecm() estimates. By the Fisher identity it is the posterior expectation of the
# complete-data score, so it comes from the same E-step quantities as em_score() (which is its sum over
# persons). Parameters are on the unconstrained scale of ecm(): log a, b, lambda, log sigma2t and the
# structural parameters.
#
# Score-based tests (Merkle & Zeileis, 2013): order the persons by a moderator, cumulate the
# decorrelated scores, and compare the process with a Brownian bridge.

#' Casewise scores of an ecm() fit
#'
#' `estfun()` returns, for a model fitted by [ecm()], the matrix of casewise scores: one row per
#' person, one column per parameter, the derivative of that person's marginal log-likelihood (ability
#' and speed integrated out). The columns sum to zero at the estimates. The parameters are on the
#' unconstrained scale used by `ecm()`: log a, b, lambda, log sigma2t, then the structural parameters
#' (for `rtirt_latreg()`/`rtirt_null()` with free speed variance: the regression coefficients, the
#' covariance `c` and the log residual variance `log_v` of speed given ability). These are the
#' estimating functions used by [score_test()] and by `strucchange::gefp()`.
#' @param x A model fitted by [ecm()] (with `prior = "ml"`).
#' @param ... Unused.
#' @return A persons x parameters matrix.
#' @exportS3Method sandwich::estfun
estfun.rtirt <- function(x, ...) {
  if (is.null(x$ecm)) stop("estfun() needs an ecm() fit: ecm(model)", call. = FALSE)
  if (isTRUE(x$ecm$map)) stop("score-based tests use the maximum likelihood fit: ecm(model), not ecm(model, prior = \"map\")", call. = FALSE)
  E <- em_setup(x, x$ecm$nodes, FALSE); p <- em_unpack(x$ecm$x, E)
  es <- em_estep(p, E); mo <- em_moments(es); N <- E$N; K <- E$K; X <- E$X
  row <- function(v) matrix(v, N, K, byrow = TRUE)
  # accuracy: sum over nodes of posterior weight x (y - P) x d eta / d par
  ab <- row(p$a * p$b); G0 <- G1 <- matrix(0, N, K)
  for (q in seq_len(ncol(es$TH))) {
    th <- es$TH[, q]; g <- es$W[, q] * (E$Y - stats::plogis(outer(th, p$a) - ab))
    G0 <- G0 + g; G1 <- G1 + th * g
  }
  S <- cbind(if (!E$one_pl) (G1 - G0 * row(p$b)) * row(p$a), -G0 * row(p$a))
  nm <- c(if (!E$one_pl) sprintf("a[%s]", E$items), sprintf("b[%s]", E$items))
  if (E$rt) {                                         # log T + zeta = lambda - rho theta + e
    rho <- if (E$eng == "cross") p$rho else rep(0, K); w <- 1 / p$sigma2t
    C <- sweep(E$L, 2, p$lambda)
    SS <- C^2 + 2 * C * mo$Ez + 2 * C * outer(mo$Et, rho) + mo$Ez2 + 2 * outer(mo$Etz, rho) + outer(mo$Et2, rho^2)
    S <- cbind(S, (C + mo$Ez + outer(mo$Et, rho)) * row(w), (-0.5 * row(w) + 0.5 * SS * row(w^2)) * row(p$sigma2t))
    nm <- c(nm, sprintf("lambda[%s]", E$items), sprintf("sigma2t[%s]", E$items))
  }
  xs <- function(prefix) if (E$P) sprintf("%s[%s]", prefix, E$xn) else character()
  str <- switch(E$eng,
    mlirt = list(X * (mo$Et - em_xb(X, p$beta)), xs("beta")),
    rtirt = {
      c <- p$c; v <- p$v; m1 <- em_xb(X, p$beta1); g <- em_xb(X, p$beta2) - c * m1
      Eu <- mo$Ez - c * mo$Et - g
      Eu2 <- mo$Ez2 - 2 * c * mo$Etz + c^2 * mo$Et2 - 2 * g * (mo$Ez - c * mo$Et) + g^2
      Eut <- mo$Etz - c * mo$Et2 - g * mo$Et
      s_c <- (Eut - m1 * Eu) / v; s_v <- -1 / (2 * v) + Eu2 / (2 * v^2)
      list(cbind(X * ((mo$Et - m1) - c / v * Eu), X * Eu / v, if (E$fixed) (1 - c^2) * (s_c - 2 * c * s_v) else cbind(s_c, v * s_v)),
           c(xs("beta_ability"), xs("beta_speed"), if (E$fixed) "atanh_cor" else c("c", "log_v")))
    },
    latent = {
      gm <- p$gamma; s <- p$s; g <- em_xb(X, p$beta)
      Eu <- mo$Ez - gm * mo$Et - g
      Eu2 <- mo$Ez2 - 2 * gm * mo$Etz + gm^2 * mo$Et2 - 2 * g * (mo$Ez - gm * mo$Et) + g^2
      Eut <- mo$Etz - gm * mo$Et2 - g * mo$Et
      list(cbind(X * Eu / s, Eut / s, s * (-1 / (2 * s) + Eu2 / (2 * s^2))), c(xs("beta"), "b_ability", "log_var_speed"))
    },
    cross = list(cbind(-(C * mo$Et + mo$Etz + outer(mo$Et2, p$rho)) * row(w),
                       if (!E$fixed) p$s * (-1 / (2 * p$s) + mo$Ez2 / (2 * p$s^2))),
                 c(sprintf("rho[%s]", E$items), if (!E$fixed) "log_var_speed")))
  S <- cbind(S, str[[1]]); colnames(S) <- c(nm, str[[2]]); rownames(S) <- NULL
  S
}

#' @rdname estfun.rtirt
#' @return `bread()`: the inverse of the average observed information (unconstrained scale), as in
#'   the 'sandwich' package.
#' @exportS3Method sandwich::bread
bread.rtirt <- function(x, ...) {
  if (is.null(x$ecm)) stop("bread() needs an ecm() fit: ecm(model)", call. = FALSE)
  E <- em_setup(x, x$ecm$nodes, FALSE); H <- em_hessian(x$ecm$x, E)
  B <- solve(-(H + t(H)) / 2) * E$N
  nm <- colnames(estfun.rtirt(x)); dimnames(B) <- list(nm, nm); B
}

#' Score-based tests of parameter invariance
#'
#' Tests whether item or structural parameters change along a person-level moderator (DIF,
#' measurement invariance, item position effects) from a single [ecm()] fit: the casewise scores
#' ([estfun.rtirt()]) are ordered by the moderator, cumulated after decorrelation, and compared with a
#' Brownian bridge (Merkle & Zeileis, 2013; for IRT, Schneider, Strobl, Zeileis & Debelak, 2022).
#'
#' **Moderators derived from the response times.** The marginal distribution of the response times
#' does not involve the accuracy parameters (a, b), so their scores have mean zero given any function
#' of the times (for example the total log time): the accuracy parameters can be tested along such a
#' moderator, provided the scores are decorrelated within the tested block (`decorrelate = "block"`,
#' the default). The response time parameters cannot be tested along a moderator built from the
#' times. Speed scores (EAPs) also use the responses and are not such a function.
#'
#' **DIF and impact.** The model identifies ability as N(x'beta, 1) for everyone, so a moderator
#' related to ability (impact) shifts all difficulties and is flagged as instability. To test DIF
#' beyond impact, include the moderator as a covariate of the latent regression
#' (`rtirt_latreg()`, `mlirt()`) and test that fit; `score_test()` gives a message when the
#' moderator correlates with the ability EAPs and is not among the covariates.
#'
#' **Clustered data.** With `cluster`, the scores are summed within clusters, which is exact for a
#' moderator that is constant within clusters (for example a school characteristic); a person-level
#' moderator with clustered persons is refused.
#' @param object A model fitted by [ecm()] (fitted first when it has no ECM fit).
#' @param moderator A vector with one value per person: numeric (continuous or ordered) or a factor.
#' @param parm Parameters to test: `"accuracy"` (a and b of every item, default), `"rt"` (lambda,
#'   sigma2t and rho), `"structural"`, a regular expression on the parameter names, or column indices.
#' @param by_item `TRUE`: test the parameters of each item separately (a and b, or the RT
#'   parameters), with Holm-adjusted p-values, in addition to the joint test.
#' @param test `"auto"` (the double-maximum test `"DM"` for numeric moderators, the Lagrange
#'   multiplier test `"LM"` for factors), `"DM"`, `"CvM"` (Cramér-von Mises), `"maxLM"`
#'   (supremum LM, 10% trimming) or `"LM"`.
#' @param decorrelate `"block"` (within the tested parameters) or `"full"` (all parameters, then
#'   the tested components are selected; not valid for moderators built from the response times).
#' @param cluster Optional cluster identifiers (see Details).
#' @param impact_check Message when the moderator is related to the ability EAPs but is not a
#'   covariate of the model (see Details).
#' @return A data frame of class `rtirt_score_test`: `parameters`, `statistic`, `p.value` and, with
#'   `by_item`, `p.holm`.
#' @examples
#' \donttest{
#' cond <- set_cond(n_subj = 400, n_item = 8)
#' dat <- sim_data(cond, sim_para(cond, "cross"), "cross")
#' fit <- ecm(rtirt_cross(dat), se = FALSE)
#' score_test(fit, rbinom(400, 1, 0.5) == 1, by_item = TRUE)      # a random grouping: no DIF
#' score_test(fit, rowSums(dat$log_t))                            # along the total log time
#' }
#' @export
score_test <- function(object, ...) UseMethod("score_test")

#' @rdname score_test
#' @param ... Arguments of the method (`moderator`, `parm`, ...).
#' @export
score_test.default <- function(object, ...) stop("score_test() needs a model of birtRcpp (fitted by ecm()) or birt", call. = FALSE)

#' @rdname score_test
#' @export
score_test.rtirt <- function(object, moderator, parm = "accuracy", by_item = FALSE, test = c("auto", "DM", "CvM", "maxLM", "LM"),
                       decorrelate = c("block", "full"), cluster = NULL, impact_check = TRUE, ...) {
  test <- match.arg(test); decorrelate <- match.arg(decorrelate)
  if (!requireNamespace("strucchange", quietly = TRUE)) stop("score_test() needs the 'strucchange' package", call. = FALSE)
  if (is.null(object$ecm)) { message("fitting ecm() first"); object <- ecm(object, se = FALSE) }
  S <- estfun.rtirt(object); N <- nrow(S); nm <- colnames(S)
  if (is.logical(moderator) || is.character(moderator)) moderator <- factor(moderator)
  if (length(moderator) != N) stop(sprintf("moderator has %d values for %d persons", length(moderator), N), call. = FALSE)
  if (anyNA(moderator)) stop("moderator has missing values; drop those persons before ecm()", call. = FALSE)
  if (impact_check) {                                  # DIF vs impact: is the moderator related to ability?
    eap <- object$ecm$person$theta$mean; zz <- as.numeric(if (is.factor(moderator)) moderator != levels(moderator)[1] else moderator)
    X <- object$ecm$X; r <- stats::cor(eap, zz)
    explained <- ncol(X) && max(abs(stats::cor(zz, X))) > 0.99
    if (!explained && abs(r) > 0.1 && stats::cor.test(eap, zz)$p.value < 0.001)
      message(sprintf(paste0("the moderator correlates %.2f with the ability EAPs: differences in mean ability (impact) show up as ",
                             "shifts of every b; to test DIF beyond impact, put the moderator in the latent regression ",
                             "(e.g. rtirt_latreg(input_data(..., cov = z))) and test that fit"), r))
  }
  if (!is.null(cluster)) {
    if (length(cluster) != N) stop("cluster needs one value per person", call. = FALSE)
    same <- tapply(seq_len(N), cluster, function(i) length(unique(as.character(moderator[i]))) == 1)
    if (!all(same)) stop("with cluster, the moderator must be constant within clusters (a person-level moderator in clustered data is not covered)", call. = FALSE)
    first <- !duplicated(cluster)
    S <- rowsum(S, cluster, reorder = FALSE); moderator <- moderator[first]; N <- nrow(S)
  }
  pick <- function(p) {
    if (is.numeric(p)) return(p)
    switch(p, accuracy = grep("^[ab]\\[", nm), rt = grep("^(lambda|sigma2t|rho)\\[", nm),
           structural = setdiff(seq_along(nm), grep("^(a|b|lambda|sigma2t|rho)\\[", nm)),
           grep(p, nm))
  }
  idx <- pick(parm); if (!length(idx)) stop("no parameters match parm", call. = FALSE)
  if (test == "auto") test <- if (is.factor(moderator)) "LM" else "DM"
  if (test == "LM" && !is.factor(moderator)) stop("test = \"LM\" is for a factor moderator", call. = FALSE)
  one <- function(j) {
    obj <- structure(list(S = if (decorrelate == "block") S[, j, drop = FALSE] else S), class = "rtirt_scores")
    pp <- if (decorrelate == "block") seq_along(j) else j
    gp <- strucchange::gefp(obj, fit = NULL, scores = function(x, ...) x$S, order.by = moderator, parm = pp, sandwich = FALSE)
    fun <- switch(test, DM = strucchange::maxBB, CvM = strucchange::meanL2BB, maxLM = strucchange::supLM(0.1),
                  LM = strucchange::catL2BB(gp))
    r <- suppressWarnings(strucchange::sctest(gp, functional = fun))
    c(statistic = unname(r$statistic), p.value = unname(r$p.value))
  }
  out <- data.frame(parameters = if (is.character(parm) && length(parm) == 1) parm else "selected",
                    n = length(idx), t(one(idx)), row.names = NULL)
  if (by_item) {
    items <- object$data$items
    per <- lapply(items, function(it) { j <- intersect(idx, grep(sprintf("\\[%s\\]$", it), nm, fixed = FALSE)); if (length(j)) c(it, length(j), one(j)) })
    per <- do.call(rbind, per[!vapply(per, is.null, TRUE)])
    pt <- data.frame(parameters = per[, 1], n = as.integer(per[, 2]), statistic = as.numeric(per[, 3]), p.value = as.numeric(per[, 4]))
    pt$p.holm <- stats::p.adjust(pt$p.value, "holm"); out$p.holm <- NA
    out <- rbind(out, pt)
  }
  structure(out, class = c("rtirt_score_test", "data.frame"), test = test, decorrelate = decorrelate,
            moderator = if (is.factor(moderator)) sprintf("factor with %d levels", nlevels(moderator)) else "numeric",
            clusters = if (!is.null(cluster)) N)
}
#' @exportS3Method stats::coef
coef.rtirt_scores <- function(object, ...) stats::setNames(rep(0, ncol(object$S)), colnames(object$S))

#' @export
print.rtirt_score_test <- function(x, digits = 3, ...) {
  cat(sprintf("Score-based test of parameter invariance (%s test, %s moderator, %s decorrelation%s)\n", attr(x, "test"),
              attr(x, "moderator"), attr(x, "decorrelate"), if (!is.null(attr(x, "clusters"))) sprintf(", %d clusters", attr(x, "clusters")) else ""))
  d <- as.data.frame(unclass(x)); attributes(d)[c("test", "decorrelate", "moderator", "clusters")] <- NULL
  d$statistic <- round(d$statistic, digits); d$p.value <- format.pval(d$p.value, digits = digits, eps = 1e-4)
  if (!is.null(d$p.holm)) d$p.holm <- ifelse(is.na(d$p.holm), "", format.pval(d$p.holm, digits = digits, eps = 1e-4))
  print(d, row.names = FALSE)
  invisible(x)
}
