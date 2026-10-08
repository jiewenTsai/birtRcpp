# Lagrange multiplier (score) tests for declared models: a term added to the linear predictor of an
# observed part, with new parameters that are 0 under the fitted (null) model. Only the null model is
# fitted: the score of the new parameters by Fisher's identity and the information by Louis' formula,
# both at the ECM estimates, with the null parameters projected out of the information.

#' Lagrange multiplier test for a declared model
#'
#' Tests whether a term added to a part of a declared model is needed, from the [ecm()] fit of the
#' model without it (the null model): `score_test(fit, add = Y ~ delta[item] * z)` tests
#' \eqn{\delta_i = 0}{delta_i = 0} in \eqn{\mathrm{logit}\,P(Y_{ji} = 1) = \eta_{ji} + \delta_i z_j}{logit P(Y_ji = 1) = eta_ji + delta_i z_j},
#' that is uniform DIF along z, for all items together and for each item. Only the null model is
#' fitted.
#'
#' **The term.** `add` is a formula `Y ~ term`: `Y` names an observed part of the model, and `term` is
#' added to its linear predictor (`bernoulli_logit()`) or to its mean (`normal()`). New parameters are
#' written with their dimension as in [block_model()] (`delta[item]`, `delta[person]`) or without one
#' for a single value (`Y ~ delta * z`: the same shift for every item). The term must be 0 when the new
#' parameters are 0, e.g. `delta[item] * z`, or `delta[item] * z * theta` for a change of the
#' discrimination (non-uniform DIF). Other names in the term are parameters of the model or data:
#' * data of the model (`block_model(Y = Y, z = z, ...)`: a vector with one value per person, or a
#'   persons x items matrix, e.g. `logT ~ delta[item] * Y` for the dependence of the times on the
#'   responses of van der Linden & Glas, 2010; see [ci_test()]);
#' * or a vector given in `data = list(z = z)`, with the same rules (no missing values).
#'
#' **The statistic** (Rao, 1948; for IRT, Glas, 1998, 1999; Glas & van der Linden, 2010). With
#' \eqn{s}{s} the score of the new parameters \eqn{\delta}{delta} at the null estimates
#' \eqn{\hat\phi}{phi-hat} (with \eqn{\delta = 0}{delta = 0}) and \eqn{I}{I} the observed information of
#' \eqn{(\phi, \delta)}{(phi, delta)},
#' \deqn{LM = s' W^{-1} s, \quad W = I_{\delta\delta} - I_{\delta\phi} I_{\phi\phi}^{-1} I_{\phi\delta},}{LM = s' W^-1 s, W = I_dd - I_dp I_pp^-1 I_pd,}
#' chi-square with as many degrees of freedom as new parameters when the null model holds; it equals the
#' likelihood ratio statistic asymptotically, without fitting the alternative. W is the information
#' about \eqn{\delta}{delta} left after the null parameters are estimated; the score of the null
#' parameters is 0 at their estimates. When a direction of the new parameters is absorbed exactly by the
#' null parameters (e.g. a shift of every item and a free mean of the latent variable), its casewise
#' scores are 0 for every person; the joint test is then made in the other directions (found from the
#' cross-product of the casewise scores) and `df` is their number (the per-element tests are not
#' affected).
#' * The score is the expected complete-data score under the E-step of the fit (Fisher's identity;
#'   Louis, 1982), as for the standard errors of [ecm()].
#' * The information is the observed information by Louis' formula (expected complete-data Hessian
#'   minus the covariance of the complete-data score), on the same quadrature grid. Falk & Monroe
#'   (2018) found that LM tests in IRT with the observed (Hessian) information keep their size, and
#'   those with the cross-product of the casewise scores reject too often.
#' * Per element: when the new parameters have one dimension (e.g. `delta[item]`), each item is tested
#'   alone (the alternative in which only that item's parameters are free): \eqn{s_i' W_{ii}^{-1} s_i}{s_i' W_ii^-1 s_i},
#'   chi-square with one degree of freedom per new parameter, with Holm-adjusted p-values (Holm, 1979).
#'   `estimate` is the one-step estimate \eqn{W_{ii}^{-1} s_i}{W_ii^-1 s_i} of that element (the expected
#'   parameter change; Sörbom, 1989), shown for one new parameter.
#'
#' With `prior = "map"` the score and the information are those of the log posterior of the null
#' parameters; the new parameters get no prior. The usual test is from `prior = "ml"`.
#'
#' **DIF and impact.** As for the score tests of the built-in models ([score_test.rtirt()]): when z is
#' related to the latent variable, a model whose latent mean does not depend on z shows the
#' difference in mean as a shift of every item. Put z in the prior of the latent variable,
#' `theta[person] ~ normal(g * z, 1)`, to test DIF beyond impact; the test of all items together then
#' has one degree of freedom less (the shift that g absorbs).
#'
#' @param object A fit of [ecm()] (or [mml()]) of a [block_model()]: the null model.
#' @param add A formula `Y ~ term` (see Details).
#' @param data Optional named list (or data frame) of person-level or persons x items data used in
#'   the term and not in the model.
#' @param ... Unused.
#' @return A data frame of class `block_lm_test`: one row for the joint test (`all`) and, when the new
#'   parameters have one dimension, one row per element (`statistic`, `df`, `p.value`, `p.holm`,
#'   `estimate`). Attributes: `add`, `score` and `W`.
#' @references Falk, C. F., & Monroe, S. (2018). On Lagrange multiplier tests in multidimensional item response theory: Information matrices and model misspecification. *Educational and Psychological Measurement, 78*(4), 653–678. \doi{10.1177/0013164417714506}
#'
#' Glas, C. A. W. (1998). Detection of differential item functioning using Lagrange multiplier tests. *Statistica Sinica, 8*(3), 647–667.
#'
#' Glas, C. A. W. (1999). Modification indices for the 2-PL and the nominal response model. *Psychometrika, 64*(3), 273–294. \doi{10.1007/BF02294296}
#'
#' Glas, C. A. W., & van der Linden, W. J. (2010). Marginal likelihood inference for a model for item responses and response times. *British Journal of Mathematical and Statistical Psychology, 63*(3), 603–626. \doi{10.1348/000711009X481360}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the EM algorithm. *Journal of the Royal Statistical Society B, 44*(2), 226–233. \doi{10.1111/j.2517-6161.1982.tb01203.x}
#'
#' Rao, C. R. (1948). Large sample tests of statistical hypotheses concerning several parameters with applications to problems of estimation. *Mathematical Proceedings of the Cambridge Philosophical Society, 44*(1), 50–57. \doi{10.1017/S0305004100023987}
#'
#' Sörbom, D. (1989). Model modification. *Psychometrika, 54*(3), 371–384. \doi{10.1007/BF02294623}
#'
#' van der Linden, W. J., & Glas, C. A. W. (2010). Statistical tests of conditional independence between responses and/or response times on test items. *Psychometrika, 75*(1), 120–139. \doi{10.1007/s11336-009-9129-9}
#' @examples
#' set.seed(1); N <- 400; K <- 6
#' z <- rbinom(N, 1, 0.5); th <- rnorm(N); b0 <- seq(-1, 1, length.out = K)
#' eta <- outer(th, rep(1, K)) - matrix(b0, N, K, byrow = TRUE)
#' eta[, 2] <- eta[, 2] + 0.8 * z                                  # DIF in item 2
#' Y <- matrix(rbinom(N * K, 1, plogis(eta)), N)
#' m <- block_model(Y = Y, z = z,
#'   theta[person] ~ normal(g * z, 1), g ~ normal(0, 1),
#'   a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
#'   Y ~ bernoulli_logit(a * (theta - b)))
#' fit <- ecm(m, prior = "ml")
#' score_test(fit, add = Y ~ delta[item] * z)
#' @export
score_test.block_model_ecm <- function(object, add, data = NULL, ...) {
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  if (missing(add)) stop("score_test(): give the term to test, e.g. add = Y ~ delta[item] * z", call. = FALSE)
  if (!isTRUE(object$converged)) stop("score_test(): the fit of the null model did not converge", call. = FALSE)
  ec <- object$ecm; alt <- bm_alternative(object$model, add, data)
  E <- bm_setup(alt$M, NULL, ec$latent, ec$nodes, ec$prior, NULL, isTRUE(ec$analytic), ec$quadrature, from_data = FALSE)
  E$flat <- alt$new
  s <- object$state[setdiff(E$fixed, alt$new)]
  for (v in alt$new) s[[v]] <- numeric(dim_len(alt$M, alt$M$pars[[v]]$dim))
  s <- s[E$fixed]
  C <- ec$C
  if (!identical(sort(E$grid), sort(names(C)))) {                     # the term changed which latent variables are on the grid
    C <- bm_centres(E, s, NULL); for (k in 1:4) C <- bm_centres(E, s, bm_estep(E, s, C))
  }
  obj <- bm_estep(E, s, C)$obj
  if (identical(names(C), names(ec$C)) && abs(obj - object$loglik) > 1e-6 * max(1, abs(object$loglik)))
    stop(sprintf("score_test(): with the new parameters at 0 the %s is %.4f, not %.4f as for the fit: the term must be 0 when its new parameters are 0",
                 object$objective, obj, object$loglik), call. = FALSE)
  g <- bm_score(E, C)(s, NULL)
  H <- bm_information(E, C)(s, NULL)
  I <- -(H + t(H)) / 2
  len <- vapply(E$fixed, function(v) length(s[[v]]), 0L)
  pos <- split(seq_len(sum(len)), rep(E$fixed, len))[E$fixed]
  d <- unlist(pos[alt$new], use.names = FALSE); n <- setdiff(seq_len(sum(len)), d)
  sd_ <- unlist(g[alt$new], use.names = FALSE)
  W <- I[d, d, drop = FALSE] - I[d, n, drop = FALSE] %*% solve(I[n, n, drop = FALSE], I[n, d, drop = FALSE])
  W <- (W + t(W)) / 2
  B <- diag(length(d))                                                # the directions of delta that the data inform
  O <- attr(H, "opg")
  if (!is.null(O)) {                                                  # a direction that the null parameters absorb exactly has casewise scores 0
    Wo <- O[d, d, drop = FALSE] - O[d, n, drop = FALSE] %*% solve(O[n, n, drop = FALSE], O[n, d, drop = FALSE])
    eo <- eigen((Wo + t(Wo)) / 2, symmetric = TRUE); keep <- eo$values > 1e-6 * max(eo$values)
    if (!all(keep)) B <- eo$vectors[, keep, drop = FALSE]
  }
  WB <- crossprod(B, W %*% B); sB <- crossprod(B, sd_)
  ev <- eigen((WB + t(WB)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 0)) stop("score_test(): the observed information of the new parameters given the others is not positive definite; the alternative is not identified", call. = FALSE)
  lm <- sum(sB * solve(WB, sB)); df <- ncol(B)
  out <- data.frame(test = "all", statistic = lm, df = df, p.value = stats::pchisq(lm, df, lower.tail = FALSE), p.holm = NA_real_, estimate = NA_real_)
  dims <- unique(vapply(alt$new, function(v) alt$M$pars[[v]]$dim, ""))
  if (length(dims) == 1 && dims != "scalar") {                         # each element alone: the alternative with only its parameters free
    nk <- dim_len(alt$M, dims); off <- seq_along(alt$new) - 1L
    per <- lapply(seq_len(nk), function(k) {
      j <- (k - 1L) + nk * off + 1L                                    # element k of each new parameter, within d
      Wk <- W[j, j, drop = FALSE]; sk <- sd_[j]
      st <- if (min(diag(Wk)) > 1e-10 * max(abs(diag(W)))) sum(sk * solve(Wk, sk)) else NA_real_   # NA: no information (e.g. z = 0 for this item)
      data.frame(test = if (length(alt$new) == 1) sprintf("%s[%d]", alt$new, k) else sprintf("%s %d", dims, k),
                 statistic = st, df = length(j), p.value = stats::pchisq(st, length(j), lower.tail = FALSE), p.holm = NA_real_,
                 estimate = if (length(j) == 1) sk / Wk[1, 1] else NA_real_)
    })
    per <- do.call(rbind, per); per$p.holm <- stats::p.adjust(per$p.value, "holm")   # NA rows are left out of the adjustment
    out <- rbind(out, per)
  }
  structure(out, class = c("block_lm_test", "data.frame"), add = sprintf("%s ~ %s + %s", alt$y, alt$what, deparse1(alt$term)),
            new = alt$new, score = sd_, W = W, objective = object$objective, df_full = length(d))
}

# the model with the term added to the linear predictor (logistic part) or the mean (normal part) of the
# observed part named on the left of add; new parameters get a placeholder normal(0, 1) prior that is
# left out of the objective (E$flat)
bm_alternative <- function(M, add, data = NULL) {
  err <- function(...) stop(sprintf(...), call. = FALSE)
  if (!inherits(add, "formula") || length(add) != 3) err("score_test(): add must be a formula 'Y ~ term', e.g. add = Y ~ delta[item] * z")
  y <- deparse1(add[[2]]); term <- add[[3]]
  fy <- Filter(function(f) f$observed && f$y == y, M$factors)
  if (!length(fy)) err("score_test(): add = %s ~ ...: %s is not observed data of the model (%s)", y, y, paste(M$used, collapse = ", "))
  fy <- fy[[1]]
  extra <- if (is.null(data)) list() else as.list(data)
  if (length(extra) && (is.null(names(extra)) || any(!nzchar(names(extra))))) err("score_test(): data must be a named list, e.g. data = list(z = z)")
  clash <- intersect(names(extra), c(names(M$data), names(M$pars), names(M$defs)))
  if (length(clash)) err("score_test(): data: %s is already a name in the model", paste(clash, collapse = ", "))
  for (n in names(extra)) if (anyNA(extra[[n]])) err("score_test(): data: %s has missing values", n)
  new <- list()
  strip <- function(e) {                                               # delta[item] -> delta, and note the new parameter
    if (is.call(e) && identical(e[[1]], as.name("["))) {
      n <- as.character(e[[2]]); d <- as.character(e[[3]])
      if (!d %in% c("person", "item")) err("score_test(): add: %s[%s]: the dimension is [person] or [item]", n, d)
      if (n %in% c(names(M$pars), names(M$data), names(extra))) err("score_test(): add: %s is already a name of the model; new parameters need new names", n)
      new[[n]] <<- d
      return(as.name(n))
    }
    if (is.call(e)) for (k in seq_along(e)[-1]) e[[k]] <- strip(e[[k]])
    e
  }
  term <- strip(term)
  if (length(intersect(all.vars(term), names(M$defs)))) term <- do.call(substitute, list(term, M$defs))
  for (n in setdiff(all.vars(term), c(names(M$pars), names(M$data), names(extra), names(new)))) new[[n]] <- "scalar"
  if (!length(new)) err("score_test(): add: the term has no new parameter (write it with its dimension, e.g. delta[item] * z)")
  for (n in names(new)) {                                              # the score of a new parameter must not vanish at 0
    dn <- tryCatch(stats::D(term, n), error = function(e) err("score_test(): add: the term must use + - * / ^ exp log sqrt"))
    both <- intersect(all.vars(dn), setdiff(names(new), n))
    if (length(both)) err("score_test(): add: the term multiplies the new parameters %s and %s (its score is 0 under the null). Data that is not in the model goes in data = list(%s = ...)",
                          n, paste(both, collapse = ", "), paste(names(new)[new == "scalar"], collapse = " = ..., "))
  }
  decls <- lapply(M$factors, bm_decl, M = M)
  k <- fy$id; what <- if (fy$dist == "bernoulli_logit") "eta" else "mean"
  decls[[k]][[3]][[2]] <- call("+", fy[[what]], term)                 # bernoulli_logit(eta + term) or normal(mean + term, sd)
  for (n in names(new)) decls[[length(decls) + 1]] <- call("~", if (new[[n]] == "scalar") as.name(n) else call("[", as.name(n), as.name(new[[n]])),
                                                             quote(normal(0, 1)))
  M1 <- bm_build(c(M$data, extra), decls, baseenv())
  list(M = M1, new = names(new), y = y, term = term, what = deparse1(fy[[what]]))
}

#' @export
print.block_lm_test <- function(x, digits = 3, ...) {
  cat(sprintf("Lagrange multiplier (score) test of %s, %s = 0 (from the %s of the null model)\n", attr(x, "add"),
              paste(attr(x, "new"), collapse = ", "), if (attr(x, "objective") == "log posterior") "posterior mode" else "ML fit"))
  if (x$df[1] < attr(x, "df_full")) cat(sprintf("all: %d of the %d directions are absorbed by the null parameters; tested in the other %d\n",
                                                  attr(x, "df_full") - x$df[1], attr(x, "df_full"), x$df[1]))
  d <- as.data.frame(unclass(x)); attributes(d)[c("add", "new", "score", "W", "objective", "df_full")] <- NULL
  d$statistic <- round(d$statistic, digits); d$p.value <- format.pval(d$p.value, digits = digits, eps = 1e-4)
  d$p.holm <- ifelse(is.na(d$p.holm), "", format.pval(d$p.holm, digits = digits, eps = 1e-4))
  d$estimate <- ifelse(is.na(d$estimate), "", format(round(d$estimate, digits)))
  if (nrow(d) == 1) d$p.holm <- d$estimate <- NULL
  print(d, row.names = FALSE)
  invisible(x)
}
