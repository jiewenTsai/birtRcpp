# Components of the EM algorithm: ecm() is built from them, and they compose EMs for new models.
# Persons j are rows, items i columns; Q quadrature nodes.

#' Components of the EM algorithm
#'
#' The steps of [ecm()] as functions that can be combined into an EM for another model (the
#' counterpart of the sampling blocks, [blocks] and [compose]: `draw_*()` draws a block from its
#' full conditional, `update_*()` maximizes, or takes a monotone MM step, in the same block).
#'
#' **E-step**
#' * `calc_nodes(Q)`: Gauss–Hermite nodes `z` and log weights `logw` of N(0, 1).
#' * `calc_loglik_nodes(Y, a, b, theta)`: the 2PL log-likelihood of every person at every node,
#'   \eqn{\sum_i y_{ij}\eta - \log(1 + e^\eta)}, \eqn{\eta = a_i(\theta - b_i)}, an N x Q matrix.
#'   `theta` is the vector of nodes when they are common to all persons (two matrix products,
#'   Bock–Aitkin), or an N x Q matrix of person-specific nodes (latent regressions; C++).
#' * `calc_normal_given_nodes(o, w, m, v)`: a normal latent variable u (e.g. speed) with prior
#'   N(m, v) at each node (`m`: N x Q, or a vector of Q), observed through normal
#'   pseudo-observations \eqn{o_{ij} = u_j + e_{ij}}, \eqn{e_{ij} \sim N(0, 1/w_i)}: integrated in
#'   closed form. Returns `logI` (N x Q, add to the log-likelihood), the conditional `mean` (N x Q) and
#'   `var` of u.
#' * `calc_node_weights(logL, logw)`: posterior weights `W` of the nodes (N x Q) and the
#'   log-likelihood `ll` of every person.
#' * `calc_artificial_data(W, Y)`: Bock–Aitkin's expected counts: `n` (Q x K, persons at each node
#'   answering each item) and `r` (Q x K, correct answers).
#' * `calc_item_sums(Y, a, b, theta, W, score)`: weighted item sums over persons and nodes, either the
#'   Pólya–Gamma minorizer sums \eqn{S_r = \sum W E[\omega|\eta]\theta^r} (`score = FALSE`) or the
#'   score sums \eqn{\sum W (y - P)\theta^r} (`score = TRUE`), r = 0, 1, 2; common nodes go through
#'   the artificial data.
#'
#' **M-step**
#' * `update_items_logistic(S, T0, T1, a, b, prior_a, prior_b)`: the Pólya–Gamma (Jaakkola–Jordan)
#'   minorizer step for 2PL items: a 2 x 2 weighted least-squares problem per item from the sums `S`
#'   of `calc_item_sums()` and \eqn{T_0 = \sum_j (y_{ij} - 1/2)}, \eqn{T_1 = \sum_j (y_{ij} - 1/2)
#'   E[\theta_j]}; the marginal log-likelihood never decreases. With `prior_a` = c(mean, sd) and
#'   `prior_b` = c(mean, sd, lower, upper) (or a list whose mean has one value per item) the conditional modes (MAP).
#'
#' **Driver and information**
#' * `em_squarem(x, Fmap, ...)`: SQUAREM (Varadhan & Roland, 2008; SqS3) around a fixed-point map
#'   `Fmap(x)` returning `list(x = , obj = )` (one EM step from x and the log-likelihood at x), with a
#'   fall-back to plain EM steps when the extrapolation lowers the objective.
#' * `calc_hessian(score, x)`: the symmetrized numerical Jacobian of an analytic score function
#'   (central differences); \eqn{(-H)^{-1}} is the covariance of the estimates.
#' @param Q Number of nodes.
#' @param Y Responses (N x K, 0/1).
#' @param a,b Item parameters.
#' @param theta Common nodes (a vector) or person-specific nodes (N x Q).
#' @param o,w,m,v Pseudo-observations (N x K), their precisions (per item), prior means (N x Q or a
#'   vector of Q) and prior variance of the normal latent variable.
#' @param logL,logw Log-likelihood at the nodes (N x Q) and log weights of the nodes.
#' @param W Posterior weights of the nodes (N x Q).
#' @param score Score sums instead of minorizer sums.
#' @param S,T0,T1 Item sums (a list with `S0`, `S1`, `S2`) and the totals of `y - 1/2`.
#' @param prior_a,prior_b Normal priors of a (on (0, Inf)) and b (MAP), or `NULL` (maximum likelihood).
#' @param x Starting value / point.
#' @param Fmap Fixed-point map.
#' @param max_iter,tol,accelerate Iteration limit, convergence tolerance (relative change of the
#'   objective and `sqrt(tol)` for the largest change of x), use SQUAREM.
#' @param progress Optional function called as `progress(iteration, objective)`.
#' @param h Relative step of the central differences.
#' @return See Details.
#' @references Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item
#'   parameters: Application of an EM algorithm. *Psychometrika, 46*, 443-459.
#'
#'   Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods for accelerating the
#'   convergence of any EM algorithm. *Scandinavian Journal of Statistics, 35*, 335-353.
#' @examples
#' set.seed(1)
#' N <- 500; K <- 10; a <- runif(K, 0.8, 2); b <- rnorm(K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(rnorm(N), a), 2, a * b))), N)
#' # an EM for the 2PL from the components
#' nd <- calc_nodes(41); kap <- Y - 0.5
#' Fmap <- function(x) {
#'   a <- exp(x[1:K]); b <- x[K + 1:K]
#'   nw <- calc_node_weights(calc_loglik_nodes(Y, a, b, nd$z), nd$logw)
#'   S <- calc_item_sums(Y, a, b, nd$z, nw$W)
#'   new <- update_items_logistic(S, colSums(kap), colSums(kap * drop(nw$W %*% nd$z)), a, b)
#'   list(x = c(log(new$a), new$b), obj = sum(nw$ll))
#' }
#' fit <- em_squarem(c(rep(0, K), rep(0, K)), Fmap)
#' round(cbind(a, a_hat = exp(fit$x[1:K]), b, b_hat = fit$x[K + 1:K]), 2)
#' @name em_components
NULL

#' @rdname em_components
#' @export
calc_nodes <- function(Q = 41) { gh <- statmod::gauss.quad.prob(Q, "normal"); list(z = gh$nodes, logw = log(gh$weights)) }

#' @rdname em_components
#' @export
calc_loglik_nodes <- function(Y, a, b, theta) {
  if (is.matrix(theta)) return(.em_loglik_nodes(as_mat(Y), as.numeric(a), as.numeric(b), as_mat(theta)))
  Y <- as.matrix(Y); M <- !is.na(Y); Y0 <- ifelse(M, Y, 0)
  eta <- outer(theta, a) - matrix(a * b, length(theta), length(a), byrow = TRUE)              # Q x K, common to all persons
  base <- outer(drop(Y0 %*% a), theta) - drop(Y0 %*% (a * b))                                 # sum_i y eta
  if (all(M)) base - matrix(rowSums(softplus(eta)), nrow(Y), length(theta), byrow = TRUE)
  else base - (M * 1) %*% t(softplus(eta))
}

#' @rdname em_components
#' @export
calc_normal_given_nodes <- function(o, w, m, v) {
  if (is.null(dim(m))) m <- matrix(m, nrow(o), length(m), byrow = TRUE)
  normal_given_sums(drop(o %*% w), drop(o^2 %*% w), w, m, v)
}

# The same integral from the sums s1 = sum_i w_i o_i and s2 = sum_i w_i o_i^2 (one value per person, or
# per person and node when o depends on the node; m conforms to them). Completing the square: u | o ~
# N(mean, var), var = 1 / (1/v + sum w); the log integral is what ecm(), marginal_loglik() and
# calc_normal_given_nodes() add for speed.
normal_given_sums <- function(s1, s2, w, m, v) {
  P <- 1 / v + sum(w); mu <- (m / v + s1) / P
  list(logI = -0.5 * length(w) * log(2 * pi) + 0.5 * sum(log(w)) - 0.5 * log(v * P) - 0.5 * (s2 + m^2 / v - P * mu^2), mean = mu, var = 1 / P)
}

#' @rdname em_components
#' @export
calc_node_weights <- function(logL, logw = 0) {
  L <- if (length(logw) == 1 && logw == 0) logL else logL + matrix(logw, nrow(logL), ncol(logL), byrow = TRUE)
  mx <- apply(L, 1, max); ll <- mx + log(rowSums(exp(L - mx)))           # log-sum-exp per person
  list(W = exp(L - ll), ll = ll)
}

#' @rdname em_components
#' @export
calc_artificial_data <- function(W, Y) {
  Y <- as.matrix(Y); M <- !is.na(Y)
  list(n = crossprod(W, M * 1), r = crossprod(W, ifelse(M, Y, 0)))
}

#' @rdname em_components
#' @export
calc_item_sums <- function(Y, a, b, theta, W, score = FALSE) {
  if (is.matrix(theta)) { r <- .em_item_sums(as_mat(Y), as.numeric(a), as.numeric(b), as_mat(theta), W, score); return(lapply(r, as.vector)) }
  ad <- calc_artificial_data(W, Y); Q <- length(theta); K <- length(a)
  eta <- outer(theta, a) - matrix(a * b, Q, K, byrow = TRUE)
  G <- if (score) ad$r - ad$n * stats::plogis(eta) else ad$n * pg_mean(eta)
  list(S0 = colSums(G), S1 = colSums(G * theta), S2 = if (score) numeric(K) else colSums(G * theta^2))
}

#' @rdname em_components
#' @export
update_items_logistic <- function(S, T0, T1, a, b, prior_a = NULL, prior_b = NULL) {
  S0 <- S$S0; S1 <- S$S1; S2 <- S$S2
  if (is.null(prior_a) && is.null(prior_b)) {             # eta = a theta - d: 2 x 2 weighted least squares, keep a > 0
    dt <- S2 * S0 - S1^2; an <- (S0 * T1 - S1 * T0) / dt; d <- (S1 * T1 - S2 * T0) / dt
    a0 <- a; d0 <- a * b
    for (k in which(!(an > 1e-3))) {                      # the minorizer rises along the segment
      t <- 1; repeat { t <- t / 2; at <- a0[k] + t * (an[k] - a0[k]); if (at > 1e-3 || t < 1e-8) break }
      an[k] <- max(at, 1e-3); d[k] <- d0[k] + t * (d[k] - d0[k])
    }
    return(list(a = an, b = d / an))
  }
  vb <- prior_b[[2]]^2; va <- prior_a[[2]]^2              # conditional modes of b, then a
  bn <- pmin(pmax((a * (a * S1 - T0) + prior_b[[1]] / vb) / (1 / vb + a^2 * S0), prior_b[[3]]), prior_b[[4]])
  nl <- T1 - bn * T0; pl <- S2 - 2 * bn * S1 + bn^2 * S0
  list(a = pmax((prior_a[[1]] / va + nl) / (1 / va + pl), 1e-3), b = bn)
}

#' @rdname em_components
#' @export
em_squarem <- function(x, Fmap, max_iter = 1000, tol = 1e-10, accelerate = TRUE, progress = NULL) {
  trace <- numeric(); conv <- FALSE; it <- 0
  r <- Fmap(x)
  while (it < max_iter) {
    it <- it + 1; trace[it] <- r$obj
    if (accelerate) {                                    # SqS3 with a monotone fall-back
      r2 <- Fmap(r$x); d1 <- r$x - x; d2 <- r2$x - r$x - d1
      al <- -sqrt(sum(d1^2) / max(sum(d2^2), 1e-300)); al <- min(al, -1)
      xa <- x - 2 * al * d1 + al^2 * d2
      ra <- tryCatch(Fmap(xa), error = function(e) NULL)
      if (!is.null(ra) && is.finite(ra$obj) && ra$obj >= r$obj) { xn <- xa; rn <- ra } else { xn <- r$x; rn <- r2 }
    } else { xn <- r$x; rn <- Fmap(xn) }
    if (!is.null(progress)) progress(it, rn$obj)
    if (!all(is.finite(xn)) || !is.finite(rn$obj)) {
      warning(sprintf("EM stopped at iteration %d: the step gave non-finite values (a parameter went to 0 or infinity)", it), call. = FALSE)
      break
    }
    dx <- max(abs(xn - x)); x <- xn
    done <- abs(rn$obj - r$obj) <= tol * abs(rn$obj) && dx <= sqrt(tol)
    r <- rn
    if (done) { conv <- TRUE; break }
  }
  trace[it + 1] <- r$obj
  list(x = x, obj = r$obj, iterations = it, converged = conv, trace = trace)
}

#' @rdname em_components
#' @export
calc_hessian <- function(score, x, h = 1e-4) {
  n <- length(x); H <- matrix(0, n, n)
  for (j in seq_len(n)) {
    e <- replace(numeric(n), j, h * max(1, abs(x[j])))
    H[, j] <- (score(x + e) - score(x - e)) / (2 * e[j])
  }
  (H + t(H)) / 2
}
