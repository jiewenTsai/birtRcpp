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
#' (for `model = "latreg"` / `"null"` with free speed variance: the regression coefficients, the
#' covariance `c` and the log residual variance `log_v` of speed given ability). These are the
#' estimating functions used by [score_test()] and by `strucchange::gefp()`.
#' @param x A model fitted by [ecm()] (with `prior = "ml"`).
#' @param ... Unused.
#' @return A persons x parameters matrix.
#' @exportS3Method sandwich::estfun
estfun.rtirt <- function(x, ...) {
  if (is.null(x$ecm)) {
    stop("estfun() needs an ecm() fit: ecm(model)", call. = FALSE)
  }
  if (isTRUE(x$ecm$map)) {
    stop("score-based tests use the maximum likelihood fit: ecm(model), not ecm(model, prior = \"map\")", call. = FALSE)
  }
  E <- em_setup(x, x$ecm$nodes, FALSE)
  p <- em_unpack(x$ecm$x, E)
  es <- em_estep(p, E)
  mo <- em_moments(es)
  by_item <- function(v) matrix(v, E$N, E$K, byrow = TRUE)

  accuracy <- accuracy_scores(p, E, es, by_item)
  scores <- accuracy$scores
  names <- accuracy$names
  if (E$rt) {
    rt <- rt_item_scores(p, E, mo, by_item)
    scores <- cbind(scores, rt$scores)
    names <- c(names, rt$names)
  }
  structural <- structural_scores(p, E, mo, by_item)
  scores <- cbind(scores, structural$scores)
  colnames(scores) <- c(names, structural$names)
  rownames(scores) <- NULL
  scores
}

# accuracy: sum over nodes of posterior weight x (y - P) x d eta / d par (log a and b)
accuracy_scores <- function(p, E, es, by_item) {
  ab <- by_item(p$a * p$b)
  sum_g <- matrix(0, E$N, E$K)
  sum_theta_g <- matrix(0, E$N, E$K)
  for (q in seq_len(ncol(es$TH))) {
    theta <- es$TH[, q]
    g <- es$W[, q] * (E$Y - stats::plogis(outer(theta, p$a) - ab))
    sum_g <- sum_g + g
    sum_theta_g <- sum_theta_g + theta * g
  }
  score_log_a <- if (!E$one_pl) (sum_theta_g - sum_g * by_item(p$b)) * by_item(p$a)
  score_b <- -sum_g * by_item(p$a)
  names_a <- if (!E$one_pl) sprintf("a[%s]", E$items)
  list(scores = cbind(score_log_a, score_b), names = c(names_a, sprintf("b[%s]", E$items)))
}

# log RT items (lambda, log sigma2t): log T + zeta = lambda - rho theta + e
rt_item_scores <- function(p, E, mo, by_item) {
  rho <- if (E$eng == "cross") p$rho else rep(0, E$K)
  w <- 1 / p$sigma2t
  centred <- sweep(E$L, 2, p$lambda)
  sq_resid <- centred^2 + 2 * centred * mo$Ez + 2 * centred * outer(mo$Et, rho) + mo$Ez2 + 2 * outer(mo$Etz, rho) +
    outer(mo$Et2, rho^2)
  score_lambda <- (centred + mo$Ez + outer(mo$Et, rho)) * by_item(w)
  score_log_sigma2t <- (-0.5 * by_item(w) + 0.5 * sq_resid * by_item(w^2)) * by_item(p$sigma2t)
  list(scores = cbind(score_lambda, score_log_sigma2t),
       names = c(sprintf("lambda[%s]", E$items), sprintf("sigma2t[%s]", E$items)))
}

# structural parameters of each engine, on the unconstrained scale of ecm()
structural_scores <- function(p, E, mo, by_item) {
  cov_names <- function(prefix) if (E$P) sprintf("%s[%s]", prefix, E$xn) else character()
  X <- E$X
  switch(E$eng,
    mlirt = list(scores = X * (mo$Et - em_xb(X, p$beta)), names = cov_names("beta")),
    rtirt = structural_scores_rtirt(p, E, mo, cov_names),
    latent = structural_scores_latent(p, E, mo, cov_names),
    cross = structural_scores_cross(p, E, mo, by_item))
}

structural_scores_rtirt <- function(p, E, mo, cov_names) {
  X <- E$X
  c <- p$c
  v <- p$v
  ability_mean <- em_xb(X, p$beta1)
  g <- em_xb(X, p$beta2) - c * ability_mean
  Eu <- mo$Ez - c * mo$Et - g
  Eu2 <- mo$Ez2 - 2 * c * mo$Etz + c^2 * mo$Et2 - 2 * g * (mo$Ez - c * mo$Et) + g^2
  Eut <- mo$Etz - c * mo$Et2 - g * mo$Et
  score_c <- (Eut - ability_mean * Eu) / v
  score_v <- -1 / (2 * v) + Eu2 / (2 * v^2)
  score_cov <- if (E$fixed) (1 - c^2) * (score_c - 2 * c * score_v) else cbind(score_c, v * score_v)
  scores <- cbind(X * ((mo$Et - ability_mean) - c / v * Eu), X * Eu / v, score_cov)
  cov_par_names <- if (E$fixed) "atanh_cor" else c("c", "log_v")
  list(scores = scores, names = c(cov_names("beta_ability"), cov_names("beta_speed"), cov_par_names))
}

structural_scores_latent <- function(p, E, mo, cov_names) {
  X <- E$X
  gm <- p$gamma
  s <- p$s
  g <- em_xb(X, p$beta)
  Eu <- mo$Ez - gm * mo$Et - g
  Eu2 <- mo$Ez2 - 2 * gm * mo$Etz + gm^2 * mo$Et2 - 2 * g * (mo$Ez - gm * mo$Et) + g^2
  Eut <- mo$Etz - gm * mo$Et2 - g * mo$Et
  scores <- cbind(X * Eu / s, Eut / s, s * (-1 / (2 * s) + Eu2 / (2 * s^2)))
  list(scores = scores, names = c(cov_names("beta"), "b_ability", "log_var_speed"))
}

structural_scores_cross <- function(p, E, mo, by_item) {
  w <- 1 / p$sigma2t
  centred <- sweep(E$L, 2, p$lambda)
  score_rho <- -(centred * mo$Et + mo$Etz + outer(mo$Et2, p$rho)) * by_item(w)
  score_speed <- if (!E$fixed) p$s * (-1 / (2 * p$s) + mo$Ez2 / (2 * p$s^2))
  speed_name <- if (!E$fixed) "log_var_speed"
  list(scores = cbind(score_rho, score_speed), names = c(sprintf("rho[%s]", E$items), speed_name))
}

#' @rdname estfun.rtirt
#' @return `bread()`: the inverse of the average observed information (unconstrained scale), as in
#'   the 'sandwich' package.
#' @exportS3Method sandwich::bread
bread.rtirt <- function(x, ...) {
  if (is.null(x$ecm)) {
    stop("bread() needs an ecm() fit: ecm(model)", call. = FALSE)
  }
  E <- em_setup(x, x$ecm$nodes, FALSE)
  hessian <- em_hessian(x$ecm$x, E)
  bread <- solve(-(hessian + t(hessian)) / 2) * E$N
  names <- colnames(estfun.rtirt(x))
  dimnames(bread) <- list(names, names)
  bread
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
#' (`ecm(input_data(..., cov = z), model = "latreg")`, or `model = "irt"` without response times)
#' and test that fit; `score_test()` gives a message when the
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
#' fit <- ecm(dat, model = "cross", se = FALSE)
#' score_test(fit, rbinom(400, 1, 0.5) == 1, by_item = TRUE)      # a random grouping: no DIF
#' score_test(fit, rowSums(dat$log_t))                            # along the total log time
#' }
#' @export
score_test <- function(object, ...) {
  UseMethod("score_test")
}

#' @rdname score_test
#' @param ... Arguments of the method (`moderator`, `parm`, ...).
#' @export
score_test.default <- function(object, ...) {
  stop("score_test() needs a model of birtRcpp (fitted by ecm()) or birt", call. = FALSE)
}

#' @rdname score_test
#' @export
score_test.rtirt <- function(object, moderator, parm = "accuracy", by_item = FALSE, test = c("auto", "DM", "CvM", "maxLM", "LM"),
                       decorrelate = c("block", "full"), cluster = NULL, impact_check = TRUE, ...) {
  test <- match.arg(test)
  decorrelate <- match.arg(decorrelate)
  if (!requireNamespace("strucchange", quietly = TRUE)) {
    stop("score_test() needs the 'strucchange' package", call. = FALSE)
  }
  if (is.null(object$ecm)) {
    message("fitting ecm() first")
    object <- ecm(object, se = FALSE)
  }
  scores <- estfun.rtirt(object)
  n_person <- nrow(scores)
  names <- colnames(scores)
  moderator <- check_moderator(moderator, n_person)
  if (impact_check) {
    impact_message(object, moderator)
  }
  if (!is.null(cluster)) {
    clustered <- cluster_scores(scores, moderator, cluster)
    scores <- clustered$scores
    moderator <- clustered$moderator
    n_person <- nrow(scores)
  }
  columns <- pick_parameters(parm, names)
  if (!length(columns)) {
    stop("no parameters match parm", call. = FALSE)
  }
  if (test == "auto") {
    test <- if (is.factor(moderator)) "LM" else "DM"
  }
  if (test == "LM" && !is.factor(moderator)) {
    stop("test = \"LM\" is for a factor moderator", call. = FALSE)
  }
  run_test <- function(j) score_process_test(scores, j, moderator, test, decorrelate)

  label <- if (is.character(parm) && length(parm) == 1) parm else "selected"
  table <- data.frame(parameters = label, n = length(columns), t(run_test(columns)), row.names = NULL)
  if (by_item) {
    table <- rbind_item_tests(table, object$data$items, columns, names, run_test)
  }
  moderator_text <- if (is.factor(moderator)) sprintf("factor with %d levels", nlevels(moderator)) else "numeric"
  structure(table, class = c("rtirt_score_test", "data.frame"), test = test, decorrelate = decorrelate,
            moderator = moderator_text, clusters = if (!is.null(cluster)) n_person)
}

check_moderator <- function(moderator, n_person) {
  if (is.logical(moderator) || is.character(moderator)) {
    moderator <- factor(moderator)
  }
  if (length(moderator) != n_person) {
    stop(sprintf("moderator has %d values for %d persons", length(moderator), n_person), call. = FALSE)
  }
  if (anyNA(moderator)) {
    stop("moderator has missing values; drop those persons before ecm()", call. = FALSE)
  }
  moderator
}

# DIF vs impact: is the moderator related to ability (and not already a covariate)?
impact_message <- function(object, moderator) {
  eap <- object$ecm$person$theta$mean
  z <- as.numeric(if (is.factor(moderator)) moderator != levels(moderator)[1] else moderator)
  X <- object$ecm$X
  r <- stats::cor(eap, z)
  explained <- ncol(X) && max(abs(stats::cor(z, X))) > 0.99
  if (!explained && abs(r) > 0.1 && stats::cor.test(eap, z)$p.value < 0.001) {
    message(sprintf(paste0("the moderator correlates %.2f with the ability EAPs: differences in mean ability (impact) show up as ",
                           "shifts of every b; to test DIF beyond impact, put the moderator in the latent regression ",
                           "(e.g. ecm(input_data(..., cov = z), model = \"latreg\"), z centred) and test that fit"), r))
  }
}

# clustered data: sum the scores within clusters (the moderator must be constant within clusters)
cluster_scores <- function(scores, moderator, cluster) {
  n_person <- nrow(scores)
  if (length(cluster) != n_person) {
    stop("cluster needs one value per person", call. = FALSE)
  }
  constant <- tapply(seq_len(n_person), cluster, function(i) length(unique(as.character(moderator[i]))) == 1)
  if (!all(constant)) {
    stop("with cluster, the moderator must be constant within clusters (a person-level moderator in clustered data is not covered)", call. = FALSE)
  }
  first <- !duplicated(cluster)
  list(scores = rowsum(scores, cluster, reorder = FALSE), moderator = moderator[first])
}

# columns of the scores tested together: a parameter group, a regular expression or indices
pick_parameters <- function(parm, names) {
  if (is.numeric(parm)) {
    return(parm)
  }
  item_pattern <- "^(a|b|lambda|sigma2t|rho)\\["
  switch(parm,
    accuracy = grep("^[ab]\\[", names),
    rt = grep("^(lambda|sigma2t|rho)\\[", names),
    structural = setdiff(seq_along(names), grep(item_pattern, names)),
    grep(parm, names))
}

# one score-based test of the columns j (Merkle & Zeileis, 2013, via strucchange)
score_process_test <- function(scores, j, moderator, test, decorrelate) {
  block <- decorrelate == "block"
  process_scores <- if (block) scores[, j, drop = FALSE] else scores
  process <- structure(list(S = process_scores), class = "rtirt_scores")
  tested <- if (block) seq_along(j) else j
  gefp <- strucchange::gefp(process, fit = NULL, scores = function(x, ...) x$S, order.by = moderator, parm = tested,
                            sandwich = FALSE)
  functional <- switch(test,
    DM = strucchange::maxBB,
    CvM = strucchange::meanL2BB,
    maxLM = strucchange::supLM(0.1),
    LM = strucchange::catL2BB(gefp))
  result <- suppressWarnings(strucchange::sctest(gefp, functional = functional))
  c(statistic = unname(result$statistic), p.value = unname(result$p.value))
}

# one test per item (its tested parameters), with Holm-adjusted p-values
rbind_item_tests <- function(table, items, columns, names, run_test) {
  one_item <- function(item) {
    j <- intersect(columns, grep(sprintf("\\[%s\\]$", item), names, fixed = FALSE))
    if (length(j)) c(item, length(j), run_test(j))
  }
  per_item <- lapply(items, one_item)
  per_item <- do.call(rbind, per_item[!vapply(per_item, is.null, TRUE)])
  item_table <- data.frame(parameters = per_item[, 1], n = as.integer(per_item[, 2]),
                           statistic = as.numeric(per_item[, 3]), p.value = as.numeric(per_item[, 4]))
  item_table$p.holm <- stats::p.adjust(item_table$p.value, "holm")
  table$p.holm <- NA
  rbind(table, item_table)
}
#' @exportS3Method stats::coef
coef.rtirt_scores <- function(object, ...) {
  stats::setNames(rep(0, ncol(object$S)), colnames(object$S))
}

#' @export
print.rtirt_score_test <- function(x, digits = 3, ...) {
  clusters <- if (!is.null(attr(x, "clusters"))) sprintf(", %d clusters", attr(x, "clusters")) else ""
  cat(sprintf("Score-based test of parameter invariance (%s test, %s moderator, %s decorrelation%s)\n", attr(x, "test"),
              attr(x, "moderator"), attr(x, "decorrelate"), clusters))
  table <- as.data.frame(unclass(x))
  attributes(table)[c("test", "decorrelate", "moderator", "clusters")] <- NULL
  table$statistic <- round(table$statistic, digits)
  table$p.value <- format.pval(table$p.value, digits = digits, eps = 1e-4)
  if (!is.null(table$p.holm)) {
    table$p.holm <- ifelse(is.na(table$p.holm), "", format.pval(table$p.holm, digits = digits, eps = 1e-4))
  }
  print(table, row.names = FALSE)
  invisible(x)
}
