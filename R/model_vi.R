# CAVI (coordinate ascent variational inference) for declared models. The third use of the same full
# conditionals: a Gibbs step draws from p(x | rest), an ECM step moves x to the mode of the averaged
# conditional, and a CAVI step sets q(x) proportional to exp E_q[log p(x | rest)], which has the form
# of the full conditional with its natural parameters replaced by their expectations (Blei,
# Kucukelbir & McAuliffe, 2017; Durante & Rigon, 2019 for the Polya-Gamma logistic part). The person
# parameters get the free-form optimum q(theta_j) on the adaptive grid of ecm() (a normal-linear
# latent such as speed is integrated exactly given the others), so vi() and ecm() share
# the E-step; with every other q a point mass vi() is ecm().

#' Variational Bayes
#'
#' `vi()` fits a model by mean-field variational Bayes (CAVI); see [vi.block_model()].
#' @param model A [block_model()].
#' @param ... Passed to the method.
#' @export
vi <- function(model, ...) UseMethod("vi")

#' @export
vi.default <- function(model, ...) stop("vi() fits declared models: model must come from block_model()", call. = FALSE)

#' Variational Bayes (CAVI) for a declared model
#'
#' Fits a [block_model()] by coordinate ascent variational inference with the same steps as
#' [gibbs()] and [ecm()]. The three engines use the same full conditionals:
#'
#' | step | Gibbs (`gibbs()`) | ECM (`ecm()`) | CAVI (`vi()`) |
#' |---|---|---|---|
#' | logistic part | draw \eqn{\omega \sim PG(1, \eta)}{omega ~ PG(1, eta)} | \eqn{E[\omega]}{E[omega]} at \eqn{\eta}{eta} | \eqn{q(\omega) = PG(1, c)}{q(omega) = PG(1, c)}, \eqn{c^2 = E_q[\eta^2]}{c^2 = E_q[eta^2]} |
#' | person parameters | draw | integrate (adaptive quadrature) | free-form \eqn{q(\theta_j)}{q(theta_j)} on the same grid |
#' | `step_gibbs(x)`, normal prior | draw \eqn{N(m/p, 1/p)}{N(m/p, 1/p)} | move to \eqn{m/p}{m/p} | \eqn{q(x) = N(E m / E p, 1 / E p)}{q(x) = N(E m / E p, 1 / E p)} (truncated if bounded) |
#' | `step_gibbs(v)`, variance | draw IG | mode | \eqn{q(v)}{q(v)} = IG with expected squared residuals |
#'
#' Expectations of products such as \eqn{E_q[(a(\theta - b))^2]}{E_q[(a (theta - b))^2]} are exact: the
#' expressions are multilinear, and under the factorized q every product is a product of first and
#' second moments. The objective is the evidence lower bound (ELBO); every update raises it. A latent
#' variable integrated exactly (normal and linear, e.g. speed) uses the steps of [ecm()] when its parts are
#' affine in all latent variables: q(u | grid) from a quadratic per person, and in the updates each person's
#' mean and covariance of the latent variables through 2L points (exact for these quadratic parts; Julier &
#' Uhlmann, 2004).
#'
#' Known behavior of mean-field VB, also here: posterior means are close to the posterior, but posterior
#' SDs of the item parameters are too small (the correlations between parameters are ignored); use
#' [gibbs()] for intervals. `step_mala()` has no CAVI update (the parameter needs a
#' conditionally conjugate full conditional) and `step_shift()` is skipped.
#'
#' @inheritParams ecm.block_model
#' @param verbose Print the ELBO every 25 sweeps.
#' @param quadrature `"adaptive"` (default): q(person) on nodes centred and scaled per person. `"ba81"`: the same
#'   nodes for every person (Bock & Aitkin, 1981; default 41, 21 or 11 per latent variable), with the logistic
#'   parts summed over persons at each node, as in [ecm()]; faster for many persons (the ELBO is that of the
#'   common nodes; use the adaptive nodes when a latent variable has a narrow posterior).
#' @return A list of class `block_model_vi`: `summary` (posterior means and SDs of q), `latent` (posterior
#'   means and SDs of the person parameters), `q`, `elbo` (trace, one value per sweep), `elbo_drop` (the
#'   largest relative drop between sweeps of one round: 0 when every update is right), `iterations`, `converged`, `time`,
#'   `model`, `steps`; `algorithm(fit)` lists the updates.
#' @references Blei, D. M., Kucukelbir, A., & McAuliffe, J. D. (2017). Variational inference: A review for
#'   statisticians. *Journal of the American Statistical Association, 112*, 859–877.
#'
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application
#' of an EM algorithm. *Psychometrika, 46*(4), 443–459. \doi{10.1007/BF02293801}
#'
#' Durante, D., & Rigon, T. (2019). Conditionally conjugate mean-field variational Bayes for logistic
#' models. *Statistical Science, 34*, 472–485.
#'
#' Julier, S. J., & Uhlmann, J. K. (2004). Unscented filtering and nonlinear estimation. *Proceedings of the
#' IEEE, 92*(3), 401–422. \doi{10.1109/JPROC.2003.823141}
#' @examples
#' set.seed(1); N <- 300; K <- 6
#' th <- rnorm(N); a0 <- runif(K, 0.8, 2); b0 <- seq(-1.5, 1.5, length.out = K)
#' Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(th, a0), 2, a0 * b0))), N)
#' m <- block_model(Y = Y,
#'   theta[person] ~ normal(0, 1),
#'   a[item] ~ normal(1, 1, lower = 0),
#'   b[item] ~ normal(0, 3),
#'   Y ~ bernoulli_logit(a * (theta - b)))
#' fit <- vi(m, list(step_gibbs(theta, a, b)))
#' fit
#' algorithm(fit)
#' @export
vi.block_model <- function(model, steps = NULL, latent = NULL, nodes = NULL, init = NULL, cycles = 6, max_iter = 1000, tol = 1e-9,
                      verbose = FALSE, quadrature = c("adaptive", "ba81"), ...) {
  if (length(list(...))) warning("unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  t0 <- Sys.time(); quadrature <- match.arg(quadrature)
  E <- vb_setup(bm_setup(model, steps, latent, nodes, "map", init, TRUE, quadrature, from_data = FALSE))   # CAVI from the prior means: the ELBO rises from the first sweep
  Q <- vb_init(E)
  C <- bm_centres(E, E$s, NULL)
  elbo <- numeric(); worst <- 0; where <- ""; last <- -Inf; conv <- FALSE
  for (cy in seq_len(cycles)) {
    start <- length(elbo)
    for (it in seq_len(max_iter)) {
      e <- vb_estep(E, Q, C); val <- sum(e$ll) + vb_rest(E, Q); elbo <- c(elbo, val)
      k <- length(elbo)
      if (k - start > 1) { d <- (elbo[k - 1] - val) / max(1, abs(val)); if (d > worst) { worst <- d; where <- sprintf("round %d, sweep %d", cy, it) } }
      if (verbose && it %% 25 == 0) message(sprintf("CAVI round %d, sweep %d: ELBO %.6f", cy, it, val))
      if (k - start > 1 && abs(val - elbo[k - 1]) < tol * abs(val)) break
      for (st in E$steps) for (v in st$vars) Q <- vb_update(E, Q, e, v)      # with q(person) of this sweep held fixed
    }
    e <- vb_estep(E, Q, C); C <- bm_centres(E, NULL, e)
    if (abs(elbo[length(elbo)] - last) < 1e-8 * abs(elbo[length(elbo)])) { conv <- TRUE; break }
    last <- elbo[length(elbo)]
  }
  e <- vb_estep(E, Q, C)
  latent_sum <- lapply(stats::setNames(E$latent, E$latent), function(v) {
    if (E$sigma && v %in% E$ana) { X <- matrix(e$mu[[v]], model$N); X2 <- X^2 + matrix(e$vv[[v]], model$N) }
    else { X <- matrix(e$X[[v]], model$N); X2 <- X^2 }
    m <- rowSums(e$W * X); data.frame(mean = m, sd = sqrt(pmax(rowSums(e$W * X2) - m^2, 0))) })
  structure(list(summary = vb_table(E, Q), latent = latent_sum, q = Q, elbo = elbo, iterations = length(elbo), converged = conv, cycles = cy,
                 elbo_drop = c(worst = worst, where = where), model = model, steps = E$steps,
                 ecm = list(latent = E$latent, nodes = E$Q, C = C, quadrature = quadrature),
                 time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), class = "block_model_vi")
}

# the moment-packed parts of every factor (expectations of squares and cross products)
vb_setup <- function(E) {
  M <- E$M
  for (st in E$steps) if (st$type == "mala")
    stop(sprintf("vi(): step_mala(%s) has no CAVI update (the full conditional is not conjugate); use step_gibbs() or ecm()/gibbs()",
                 paste(st$vars, collapse = ", ")), call. = FALSE)
  for (f in M$factors) {
    if (f$dist %in% c("half_t", "inv_gamma")) next
    if (is.null(f$ir) || is.null(f$ir[[if (f$dist == "normal") "r" else "eta"]]))
      stop(sprintf("vi(): '%s' is not multilinear (no exp(), log(), division by a parameter); use ecm() or gibbs()", deparse1(f$decl)), call. = FALSE)
    if (f$dist == "normal" && !(is.numeric(f$sd) || (is.call(f$sd) && identical(f$sd[[1]], as.name("sqrt")) && is.name(f$sd[[2]]) &&
                                                     as.character(f$sd[[2]]) %in% names(M$pars) && M$pars[[as.character(f$sd[[2]])]]$dist != "normal")))
      stop(sprintf("vi(): the sd in '%s' must be a number or sqrt(a variance parameter)", deparse1(f$decl)), call. = FALSE)
    ir <- if (f$dist == "normal") f$ir$r$ir else f$ir$eta$ir
    p <- list(sq = ml_pack_mom(ml_mul(ir, ir), M), lin = ml_pack_mom(ir, M), c1 = list(), cc = list(), cr = list())
    for (v in intersect(f$vars, names(M$pars))) {
      c1 <- ml_deriv(ir, v); r0 <- ml_drop(ir, v)
      p$c1[[v]] <- ml_pack_mom(c1, M); p$cc[[v]] <- ml_pack_mom(ml_mul(c1, c1), M); p$cr[[v]] <- ml_pack_mom(ml_mul(c1, r0), M)
    }
    if (f$dist == "normal") p$var <- if (is.numeric(f$sd)) NULL else as.character(f$sd[[2]])
    E$vb[[f$id]] <- p
  }
  E
}

vb_init <- function(E) {
  M <- E$M; s0 <- E$s; Q <- list(norm = list(), var = list())
  for (v in E$fixed) {
    P <- M$pars[[v]]; x <- s0[[v]]
    if (P$dist == "normal") Q$norm[[v]] <- list(m = x, V = 0 * x + 0.01, H = 0 * x)
    else Q$var[[v]] <- list(alpha = 0 * x + 3, beta = 2 * x, au = 0 * x + 1, bu = 0 * x + 1)
  }
  vb_moments(Q)
}
vb_moments <- function(Q) {
  Q$mom <- list(m = lapply(Q$norm, `[[`, "m"), m2 = lapply(Q$norm, function(q) q$m^2 + q$V))
  Q
}

# moments of every parameter on the stacked state: q moments for the others, node values for the latent
vb_mom <- function(E, Q, X, reps) {
  m <- Q$mom$m; m2 <- Q$mom$m2
  for (v in names(m)) if (E$M$pars[[v]]$dim == "person") { m[[v]] <- rep(m[[v]], reps); m2[[v]] <- rep(m2[[v]], reps) }
  for (v in names(X)) { m[[v]] <- X[[v]]; m2[[v]] <- X[[v]]^2 }
  list(m = m, m2 = m2)
}

# E[1/sd^2] and E[log sd^2] of a normal factor, in its shape
vb_wl <- function(E, Q, f, M) {
  p <- E$vb[[f$id]]
  if (is.null(p$var)) { w <- 1 / f$sd^2; return(list(w = full_of(w, f, M), l = full_of(log(1 / w), f, M))) }
  q <- Q$var[[p$var]]; pk <- f$ir$sd2
  list(w = ml_eval(pk, stats::setNames(list(q$alpha / q$beta), p$var), M, f$shape),
       l = ml_eval(pk, stats::setNames(list(log(q$beta) - digamma(q$alpha)), p$var), M, f$shape))
}

# E_q[log density] of a factor, one value per (stacked) person; om: c of q(omega) for a logistic part
vb_rows <- function(E, Q, f, mom, M) {
  fx <- M$factors[[f$id]]; p <- E$vb[[f$id]]; m <- fmask(fx, M)
  v <- if (fx$dist == "normal") { wl <- vb_wl(E, Q, fx, M); m * (-0.5 * wl$w * ml_eval_mom(p$sq, mom, M, fx$shape) - 0.5 * wl$l - 0.5 * log(2 * pi)) }
       else { cc <- sqrt(pmax(ml_eval_mom(p$sq, mom, M, fx$shape), 0)); m * (fx$kappa * ml_eval_mom(p$lin, mom, M, fx$shape) - cc / 2 - log1p(exp(-cc))) }
  if (fx$shape == "NK") rowSums(v) else v
}

# the optimal q(person parameters) on the frozen grid C: weights, node values, sum_j log normalizer
vb_estep <- function(E, Q, C) {
  N <- E$M$N; ng <- E$ng; ns <- E$ns; Xg <- list(); corr <- 0
  for (l in seq_along(E$grid)) {
    v <- E$grid[l]; z <- rep(E$Z[, l], each = N); sdv <- rep(C[[v]]$sd, ng)
    Xg[[v]] <- rep(C[[v]]$mu, ng) + sdv * z; corr <- corr + log(sdv) + z^2 / 2 + 0.5 * log(2 * pi)
  }
  mg <- vb_mom(E, Q, Xg, ng); ll <- corr
  for (f in E$fgrid) if (!(f$id %in% E$collapse)) ll <- ll + vb_rows(E, Q, f, mg, E$Mg)
  if (length(E$collapse)) { sn <- bm_nodes(E, E$s, C)[E$grid]; ll <- ll + as.vector(vb_collapsed_ll(E, Q, sn)) }
  mu <- list(); vv <- list()
  for (u in E$ana) {                                                    # q(u | grid) is normal: E_q log p is quadratic in u
    if (E$affine) { a <- vb_exact_affine(E, Q, u, Xg); mu[[u]] <- a$mu; vv[[u]] <- a$vv; ll <- ll + a$ll; next }
    mu0 <- mg; mu0$m[[u]] <- rep(0, N * ng); mu0$m2[[u]] <- rep(0, N * ng); p <- 0; h <- 0
    for (f in E$fana[[u]]) {
      fx <- E$Mg$factors[[f$id]]; pk <- E$vb[[f$id]]; w <- fmask(fx, E$Mg) * vb_wl(E, Q, fx, E$Mg)$w
      p <- p + reduce_to(w * ml_eval_mom(pk$cc[[u]], mu0, E$Mg, fx$shape), "person", fx)
      h <- h - reduce_to(w * ml_eval_mom(pk$cr[[u]], mu0, E$Mg, fx$shape), "person", fx)
    }
    mu[[u]] <- h / p; vv[[u]] <- 1 / p
    mu0$m[[u]] <- mu[[u]]; mu0$m2[[u]] <- mu[[u]]^2
    for (f in E$fana[[u]]) ll <- ll + vb_rows(E, Q, f, mu0, E$Mg)
    ll <- ll + 0.5 * log(2 * pi / p)
  }
  LL <- matrix(ll, N, ng) + matrix(E$lw, N, ng, byrow = TRUE)
  mx <- apply(LL, 1, max); lj <- mx + log(rowSums(exp(LL - mx))); Wg <- exp(LL - lj)
  if (E$sigma) out <- list(W = Wg, ll = lj, X = Xg, mu = mu, vv = vv, sig = bm_sigma(E, Wg, Xg, mu, vv))   # exact variables through their moments
  else {
    X <- lapply(Xg, function(x) rep(x, ns))
    for (k in seq_along(E$ana)) { u <- E$ana[k]; X[[u]] <- as.vector(outer(mu[[u]], rep(1, ns)) + outer(sqrt(vv[[u]]), E$S[, k])) }
    out <- list(W = matrix(rep(Wg / ns, ns), N), ll = lj, X = X)
  }
  if (length(E$collapse)) { out$Mc <- vb_collapse(E, Wg); out$Xc <- sn }
  out
}

# The exact step for affine parts (as bm_exact_affine() for ECM): every part of u is linear in the latent
# variables jointly with a variance free of them, so g(x) = E_q[log p(parts of u) | latent values x] is a
# quadratic polynomial in x = (u, the grid variables of these parts), one per person. Its coefficients come from
# g at 0, at -/+ each unit vector and at the sums of two unit vectors (sums over items once per person); at each
# node g in u is alpha + beta u - p u^2 / 2, so q(u | node) = N(beta / p, 1 / p) and the integral over u is
# alpha + beta^2 / (2 p) + log(2 pi / p) / 2.
vb_exact_affine <- function(E, Q, u, Xg) {
  M <- E$M; N <- M$N; n <- N * E$ng
  gl <- intersect(E$grid, unique(unlist(lapply(E$fana[[u]], `[[`, "vars")))); lv <- c(u, gl); L <- length(lv)
  g <- function(x) { mom <- vb_mom(E, Q, stats::setNames(lapply(x, function(v) rep(v, N)), lv), 1)
    out <- 0; for (f in E$fana[[u]]) out <- out + vb_rows(E, Q, f, mom, M); out }
  unit <- function(k, h = 1) replace(numeric(L), k, h)
  c0 <- g(numeric(L)); gp <- lapply(seq_len(L), function(k) g(unit(k))); gm <- lapply(seq_len(L), function(k) g(unit(k, -1)))
  b <- lapply(seq_len(L), function(k) (gp[[k]] - gm[[k]]) / 2)
  A <- matrix(list(), L, L)
  for (k in seq_len(L)) A[[k, k]] <- gp[[k]] + gm[[k]] - 2 * c0
  for (k in seq_len(L)) for (l in seq_len(L)[-seq_len(k)]) A[[k, l]] <- A[[l, k]] <- g(unit(k) + unit(l)) - c0 - b[[k]] - b[[l]] - (A[[k, k]] + A[[l, l]]) / 2
  rp <- function(x) rep_len(x, n)
  p <- rp(-A[[1, 1]]); be <- rp(b[[1]]); al <- rp(c0)
  for (k in seq_along(gl) + 1L) { x <- Xg[[lv[k]]]; be <- be + rp(A[[1, k]]) * x; al <- al + rp(b[[k]]) * x
    for (l in seq_along(gl) + 1L) al <- al + 0.5 * rp(A[[k, l]]) * x * Xg[[lv[l]]] }
  list(mu = be / p, vv = 1 / p, ll = al + be^2 / (2 * p) + 0.5 * log(2 * pi / p))
}

# Common nodes (quadrature = "ba81"): E_q[log p] of a collapsed part at every person and node (N x nodes), by
# matrix products. Logistic: sum_i (y - 1/2) E[eta] - c/2 - log(1 + exp(-c)), c^2 = E[eta^2]; normal:
# -1/2 E[1/sd^2] (y^2 + 2 y E[r | y = 0] + E[r^2 | y = 0]) - 1/2 E[log sd^2] - 1/2 log(2 pi).
vb_collapsed_ll <- function(E, Q, sn) {
  LL <- 0; top <- seq_len(E$ng); Mc <- E$Mc0; mom <- vb_mom(E, Q, sn, 1)
  for (id in E$collapse) { fx <- Mc$factors[[id]]; cs <- E$cs[[as.character(id)]]; p <- E$vb[[id]]
    at <- function(x) full_of(x, fx, Mc)[top, , drop = FALSE]
    if (fx$dist == "bernoulli_logit") {
      lin <- at(ml_eval_mom(p$lin, mom, Mc, fx$shape)); cc <- sqrt(pmax(at(ml_eval_mom(p$sq, mom, Mc, fx$shape)), 0))
      LL <- LL + cs$Y %*% t(lin) - cs$M %*% t(lin / 2 + cc / 2 + log1p(exp(-cc)))
    } else {
      wl <- vb_wl(E, Q, fx, Mc); w <- at(wl$w); l <- at(wl$l)
      r0 <- at(ml_eval_mom(p$lin, mom, Mc, fx$shape)); r2 <- at(ml_eval_mom(p$sq, mom, Mc, fx$shape))     # the data of Mc0 are 0
      LL <- LL - 0.5 * (cs$Y2 %*% t(w) + 2 * cs$Y %*% t(w * r0) + cs$M %*% t(w * r2 + l + log(2 * pi)))
    }
  }
  LL
}
# the artificial data of the collapsed parts (bm_collapse()); a CAVI update multiplies kappa by the mask, so kappa = y - 1/2
vb_collapse <- function(E, Wg) {
  Mc <- bm_collapse(E, Wg)
  for (id in E$collapse) if (Mc$factors[[id]]$dist == "bernoulli_logit") Mc$factors[[id]]$kappa <- Mc$data[[Mc$factors[[id]]$y]] - 0.5
  Mc
}

# E_q[log prior] of the non-latent parameters plus the entropies of their q
vb_rest <- function(E, Q) {
  M <- E$M; tot <- 0
  for (v in E$fixed) {
    P <- M$pars[[v]]; f <- E$pf[[v]]
    if (P$dist == "normal") {
      tot <- tot + sum(vb_rows(E, Q, f, Q$mom, M)) + sum(Q$norm[[v]]$H)
      if (is.finite(P$lower) || is.finite(P$upper)) {                   # normalizing constant of a truncated prior (constant mean and sd)
        mu <- fval(f, Q$mom$m, M, "mean"); sdv <- fval(f, Q$mom$m, M, "sd")
        tot <- tot - sum(full_of(log(stats::pnorm((P$upper - mu) / sdv) - stats::pnorm((P$lower - mu) / sdv)), f, M))
      }
    } else {
      q <- Q$var[[v]]; El <- log(q$beta) - digamma(q$alpha); Ei <- q$alpha / q$beta
      if (P$dist == "inv_gamma") tot <- tot + sum(P$hyper[1] * log(P$hyper[2]) - lgamma(P$hyper[1]) - (P$hyper[1] + 1) * El - P$hyper[2] * Ei)
      else {
        nu <- P$hyper[1]; A <- P$hyper[2]; Elu <- log(q$bu) - digamma(q$au); Eiu <- q$au / q$bu
        tot <- tot + sum((nu / 2) * (log(nu) - Elu) - lgamma(nu / 2) - (nu / 2 + 1) * El - nu * Eiu * Ei) +
          sum(-log(A^2) / 2 - lgamma(0.5) - 1.5 * Elu - Eiu / A^2) + sum(vb_ig_entropy(q$au, q$bu))
      }
      tot <- tot + sum(vb_ig_entropy(q$alpha, q$beta))
    }
  }
  tot
}
vb_ig_entropy <- function(a, b) a + log(b) + lgamma(a) - (1 + a) * digamma(a)

# q(x) for a normal prior: truncated normal with location h/p and precision p
vb_qnormal <- function(h, p, lo, hi) {
  mu <- h / p; sg <- 1 / sqrt(p)
  if (!is.finite(lo) && !is.finite(hi)) return(list(m = mu, V = 1 / p, H = 0.5 * log(2 * pi * exp(1)) + log(sg)))
  al <- (lo - mu) / sg; be <- (hi - mu) / sg
  lZ <- if (!is.finite(hi)) stats::pnorm(al, lower.tail = FALSE, log.p = TRUE) else if (!is.finite(lo)) stats::pnorm(be, log.p = TRUE)
        else log(pmax(stats::pnorm(be) - stats::pnorm(al), 1e-300))
  ra <- if (is.finite(lo)) exp(stats::dnorm(al, log = TRUE) - lZ) else 0
  rb <- if (is.finite(hi)) exp(stats::dnorm(be, log = TRUE) - lZ) else 0
  ta <- if (is.finite(lo)) al * ra else 0; tb <- if (is.finite(hi)) be * rb else 0
  list(m = mu + sg * (ra - rb), V = pmax(sg^2 * (1 + ta - tb - (ra - rb)^2), 1e-12 * sg^2),
       H = 0.5 * log(2 * pi * exp(1)) + log(sg) + lZ + (ta - tb) / 2)
}

# one CAVI update of q(v), with the person weights of the last E-step
vb_update <- function(E, Q, e, v) {
  M <- E$M; P <- M$pars[[v]]
  where <- function(f) {
    if (f$id %in% E$collapse) return(list(M = e$Mc, f = e$Mc$factors[[f$id]], mom = vb_mom(E, Q, e$Xc, 1), w = NULL, lat = FALSE))
    if (E$sigma && any(f$vars %in% E$ana))                              # a part of an exact variable: at the 2L points of each person
      return(list(M = E$Msig, f = E$Msig$factors[[f$id]], mom = vb_mom(E, Q, e$sig$X, e$sig$n), w = e$sig$w, lat = TRUE))
    if (any(f$vars %in% E$latent)) list(M = E$Ms, f = E$Ms$factors[[f$id]], mom = vb_mom(E, Q, e$X, E$ng * E$ns), w = as.vector(e$W), lat = TRUE)
    else list(M = M, f = f, mom = Q$mom, w = NULL, lat = FALSE)
  }
  if (P$dist == "normal") {
    prec <- 0; num <- 0
    for (f in factors_of(M, v)) {
      g <- where(f); fx <- g$f; pk <- E$vb[[f$id]]; m <- fmask(fx, g$M); if (g$lat) m <- m * g$w
      cc <- ml_eval_mom(pk$cc[[v]], g$mom, g$M, fx$shape); cr <- ml_eval_mom(pk$cr[[v]], g$mom, g$M, fx$shape)
      if (fx$dist == "normal") { w <- m * vb_wl(E, Q, fx, g$M)$w; Pp <- w * cc; H <- -w * cr }
      else { om <- pg_mean(sqrt(pmax(ml_eval_mom(pk$sq, g$mom, g$M, fx$shape), 0)))     # q(omega) updated first: c^2 = E[eta^2]
             Pp <- m * om * cc; H <- m * (fx$kappa * ml_eval_mom(pk$c1[[v]], g$mom, g$M, fx$shape) - om * cr) }
      prec <- prec + bm_fold(E, reduce_to(Pp, P$dim, fx), P$dim, g$lat); num <- num + bm_fold(E, reduce_to(H, P$dim, fx), P$dim, g$lat)
    }
    Q$norm[[v]] <- vb_qnormal(num, prec, P$lower, P$upper)
    return(vb_moments(Q))
  }
  n <- dim_len(M, P$dim); al <- numeric(n); be <- numeric(n)
  for (f in factors_of(M, v)) {
    if (!f$observed && f$y == v) next
    g <- where(f); fx <- g$f; m <- fmask(fx, g$M); if (g$lat) m <- m * g$w
    al <- al + bm_fold(E, reduce_to(m * full_of(1, fx, g$M), P$dim, fx), P$dim, g$lat) / 2
    be <- be + bm_fold(E, reduce_to(m * ml_eval_mom(E$vb[[f$id]]$sq, g$mom, g$M, fx$shape), P$dim, fx), P$dim, g$lat) / 2
  }
  q <- Q$var[[v]]
  Q$var[[v]] <- if (P$dist == "inv_gamma") list(alpha = P$hyper[1] + al, beta = P$hyper[2] + be, au = q$au, bu = q$bu) else {
    nu <- P$hyper[1]; A <- P$hyper[2]
    au <- rep((nu + 1) / 2, n); bu <- nu * q$alpha / q$beta + 1 / A^2                    # q(u) given q(v) (Huang & Wand auxiliary)
    list(alpha = nu / 2 + al, beta = nu * au / bu + be, au = au, bu = bu) }
  Q
}

vb_table <- function(E, Q) {
  nm <- function(v, n) if (n == 1) v else sprintf("%s[%d]", v, seq_len(n))
  do.call(rbind, c(
    lapply(names(Q$norm), function(v) data.frame(parameter = nm(v, length(Q$norm[[v]]$m)), mean = Q$norm[[v]]$m, sd = sqrt(Q$norm[[v]]$V))),
    lapply(names(Q$var), function(v) { q <- Q$var[[v]]
      data.frame(parameter = nm(v, length(q$alpha)), mean = ifelse(q$alpha > 1, q$beta / (q$alpha - 1), NA),
                 sd = ifelse(q$alpha > 2, q$beta / ((q$alpha - 1) * sqrt(q$alpha - 2)), NA)) })))
}

#' @rdname vi.block_model
#' @param x,object A fit of `vi()`.
#' @param ... Unused.
#' @export
print.block_model_vi <- function(x, ...) {
  cat(sprintf("CAVI for the declared model: %s after %d sweeps (%d rounds), ELBO %.4f, %.1f s\n",
              if (x$converged) "converged" else "NOT converged", x$iterations, x$cycles, x$elbo[length(x$elbo)], x$time))
  print(format(utils::head(x$summary, 30), digits = 3), row.names = FALSE)
  if (nrow(x$summary) > 30) cat(sprintf("  ... %d more (fit$summary)\n", nrow(x$summary) - 30))
  cat("Means and SDs of q. Mean-field SDs of item parameters are too small; use gibbs() for intervals.\n")
  invisible(x)
}

#' Check a CAVI fit
#'
#' The CAVI counterpart of [check_ecm()] and the Geweke check: coordinate ascent is right when
#'
#' 1. **the ELBO never decreases** from one sweep to the next (from `n_start` starting values drawn
#'    from the declared priors, on a fixed grid); and
#' 2. **the fit is a maximum of the ELBO**: moving the mean of any q (normal factors) a little up or
#'    down lowers it.
#'
#' @inheritParams vi.block_model
#' @param n_start Number of starting values.
#' @param sweeps Sweeps per start in check 1.
#' @param seed Random seed.
#' @return A data frame with one row per check (class `em_check`, printed like [check_ecm()]).
#' @export
check_vi <- function(model, steps = NULL, latent = NULL, nodes = NULL, n_start = 5, sweeps = 20, seed = NULL, quadrature = c("adaptive", "ba81")) {
  if (!is.null(seed)) set.seed(seed)
  quadrature <- match.arg(quadrature)
  E <- vb_setup(bm_setup(model, steps, latent, nodes, "map", NULL, TRUE, quadrature))
  elbo <- function(Q, C) sum(vb_estep(E, Q, C)$ll) + vb_rest(E, Q)
  drop <- 0
  for (r in seq_len(n_start)) {
    s <- simulate_model(E$M)[E$fixed]
    for (v in E$fixed) { P <- E$M$pars[[v]]
      if (P$dist != "normal") s[[v]] <- pmin(pmax(s[[v]], 0.05), 20) else if (is.finite(P$lower)) s[[v]] <- pmax(s[[v]], P$lower + 0.2) }
    E$s <- s; Q <- vb_init(E); C <- bm_centres(E, s, NULL); C <- bm_centres(E, NULL, vb_estep(E, Q, C))
    prev <- elbo(Q, C)
    for (k in seq_len(sweeps)) {
      e <- vb_estep(E, Q, C)
      for (st in E$steps) for (v in st$vars) Q <- vb_update(E, Q, e, v)
      cur <- elbo(Q, C); drop <- max(drop, (prev - cur) / max(1, abs(prev))); prev <- cur
    }
  }
  fit <- vi(model, steps, latent, nodes, quadrature = quadrature); Q <- fit$q; C <- fit$ecm$C; f0 <- elbo(Q, C); gain <- 0
  for (v in names(Q$norm)) for (k in seq_along(Q$norm[[v]]$m)) for (sgn in c(-1, 1)) {
    Q2 <- Q; h <- 0.05 * sqrt(Q$norm[[v]]$V[k]); Q2$norm[[v]]$m[k] <- Q$norm[[v]]$m[k] + sgn * h
    gain <- max(gain, (elbo(vb_moments(Q2), C) - f0) / max(1, abs(f0)))
  }
  for (v in names(Q$var)) for (slot in c("beta", if (E$M$pars[[v]]$dist == "half_t") "bu")) for (k in seq_along(Q$var[[v]][[slot]])) for (sgn in c(-1, 1)) {
    Q2 <- Q; Q2$var[[v]][[slot]][k] <- Q$var[[v]][[slot]][k] * (1 + sgn * 0.01)      # the scale of q of a variance (and of its auxiliary)
    gain <- max(gain, (elbo(Q2, C) - f0) / max(1, abs(f0)))
  }
  out <- data.frame(check = c("the ELBO never decreases from sweep to sweep (relative drop)",
                              "the fit is a maximum: moving a q mean or the scale of a variance q raises the ELBO by (relative)"),
                    worst = c(drop, gain), limit = c(1e-8, 1e-8), detail = "")
  out$passed <- out$worst <= out$limit
  structure(out[, c("check", "worst", "limit", "passed", "detail")], class = c("em_check", "data.frame"), passed = all(out$passed), n_start = n_start,
            what = "CAVI check")
}
