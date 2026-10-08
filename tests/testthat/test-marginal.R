# R reference of the closed-form cross-ALD integral in src/marginal.cpp:
# log int prod_k ALD(b_ik + zeta; q, s_k) N(zeta; 0, v) dzeta, for every row i of B (N x K)
log_ald_rt_ref <- function(B, s, v, q) {
  N <- nrow(B); K <- ncol(B); sd <- sqrt(v)
  kap <- -B; ord <- t(apply(kap, 1, order)); idx <- cbind(rep(seq_len(N), K), as.vector(ord))
  ks <- matrix(kap[idx], N, K)                                  # sorted kinks
  is <- matrix((1 / s)[as.vector(ord)], N, K); bs <- matrix(B[idx], N, K) * is
  above <- function(M) cbind(rowSums(M), rowSums(M) - t(apply(M, 1, cumsum)))   # N x (K + 1)
  al <- -q * sum(1 / s) + above(is)                             # slope in zeta
  be <- -q * rowSums(B %*% diag(1 / s, K)) + above(bs)          # intercept
  lo <- cbind(-Inf, ks); hi <- cbind(ks, Inf)
  m <- al * v
  za <- (lo - m) / sd; zb <- (hi - m) / sd
  lp <- ifelse(za > 0, stats::pnorm(za, lower.tail = FALSE, log.p = TRUE) + log1p(-exp(stats::pnorm(zb, lower.tail = FALSE, log.p = TRUE) - stats::pnorm(za, lower.tail = FALSE, log.p = TRUE))),
               stats::pnorm(zb, log.p = TRUE) + log1p(-exp(stats::pnorm(za, log.p = TRUE) - stats::pnorm(zb, log.p = TRUE))))
  L <- be + al^2 * v / 2 + lp
  sum(log(q * (1 - q) / s)) + calc_node_weights(L)$ll
}

# brute-force 2-D integration over (theta, zeta) of the complete-data likelihood
brute <- function(fit, p, persons, nz = 321) {
  D <- fit$data; X <- fit$post$settings$X; eng <- birtRcpp:::.model_info[[class(fit)[1]]]$engine
  qr <- isTRUE(fit$qr); q <- fit$post$settings$q; K <- ncol(D$Y)
  xb <- function(beta, i) if (length(beta) && ncol(X)) sum(X[i, ] * beta) else 0
  ald <- function(u, s) log(q * (1 - q) / s) - u * (q - (u < 0)) / s
  th <- seq(-7, 7, length.out = 281); ht <- diff(th[1:2])
  sapply(persons, function(i) {
    mt <- switch(eng, mlirt = xb(p$beta, i), rtirt = xb(p$beta_theta, i), 0)
    irt <- sapply(th + mt, function(t) { e <- p$a * (t - p$b); sum(D$Y[i, ] * e - log1p(exp(e))) })
    lt <- irt + dnorm(th, 0, 1, log = TRUE)
    if (eng == "mlirt") return(log(sum(exp(lt)) * ht))
    vs <- p$var_speed
    zz <- seq(-8, 8, length.out = nz) * sqrt(vs); hz <- diff(zz[1:2])
    tot <- sapply(seq_along(th), function(k) {
      t <- th[k] + mt
      m <- switch(eng, rtirt = xb(p$beta_speed, i) + p$cor * sqrt(vs) * (t - mt), latent = xb(p$beta, i) + p$b_theta * t, cross = 0)
      v <- switch(eng, rtirt = vs - p$cor^2 * vs, vs)
      z <- m + zz * sqrt(v / vs)
      lz <- if (eng == "latent" && qr) ald(z - m, vs) else dnorm(z, m, sqrt(v), log = TRUE)
      rt <- sapply(z, function(zeta) {
        u <- D$log_t[i, ] - p$lambda + zeta + if (eng == "cross") p$rho * t else 0
        if (eng == "cross" && qr) sum(ald(u, p$sigma2t)) else sum(dnorm(u, 0, sqrt(p$sigma2t), log = TRUE))
      })
      a <- rt + lz; mx <- max(a); lt[k] + mx + log(sum(exp(a - mx)) * hz * sqrt(v / vs))
    })
    mx <- max(tot); mx + log(sum(exp(tot - mx)) * ht)
  })
}

test_that("marginal log-likelihood equals brute-force integration for every model", {
  set.seed(5)
  for (cfg in list(list("mlirt", NULL), list("null", NULL), list("latreg", NULL), list("latent", NULL),
                   list("latent", 0.3), list("cross", NULL), list("cross", 0.7))) {
    type <- cfg[[1]]; qv <- cfg[[2]]
    cond <- set_cond(n_subj = 60, n_item = 5, n_feat = 2)
    dat <- sim_data(cond, sim_para(cond, type), type)
    ctor <- switch(type, mlirt = mlirt, null = rtirt_null, latreg = rtirt_latreg, latent = rtirt_latent, cross = rtirt_cross)
    m <- if (is.null(qv)) ctor(dat) else ctor(dat, quantile = qv)
    f <- gibbs(m, n_iter = 200, n_chain = 1, seed = 1, verbose = FALSE, priors = rtirt_priors_julia())   # a formula check: short chains, stable priors
    M <- birtRcpp:::post_draws(f, 1); p <- birtRcpp:::draw_pars(M[1, ], f$data$items, colnames(f$post$settings$X) %||% character())
    ours <- birtRcpp:::marg_ll_draw(f, p, Q = if (type == "cross" && !is.null(qv)) 201 else 61)[1:4]
    br <- brute(f, p, 1:4, nz = if (is.null(qv)) 321 else 3201); cat(sprintf("%s %s: max |diff| %.2e\n", type, format(qv), max(abs(ours - br))))
    expect_lt(max(abs(ours - br)), if (is.null(qv)) 1e-6 else 2e-4)
  }
})

test_that("C++ and R versions of the cross-ALD integral agree; fit indices are marginal", {
  B <- matrix(stats::rnorm(40), 8, 5); s <- stats::runif(5, 0.2, 1)
  expect_equal(birtRcpp:::.log_ald_rt(B, s, 0.7, 0.3), log_ald_rt_ref(B, s, 0.7, 0.3), tolerance = 1e-10)
  cond <- set_cond(n_subj = 200, n_item = 6, n_feat = 2)
  dat <- sim_data(cond, sim_para(cond, "latent"), "latent")
  f <- gibbs(rtirt_latent(dat), n_iter = 600, n_chain = 2, seed = 1, verbose = FALSE)
  fi <- fit_indices(f)
  expect_true(all(c("DIC", "pD", "WAIC", "LOOIC", "DIC_complete") %in% names(fi)))
  expect_true(fi$pD > 0 && fi$pD < 60)             # about the number of item and structural parameters (27)
  expect_lt(fi$DIC, fi$DIC_complete)              # complete-data deviance includes the person parameters
  g <- gibbs(rtirt_latent(dat), n_iter = 300, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
  expect_null(g$post$fit); expect_message(fit_indices(g), "marginal")
})
