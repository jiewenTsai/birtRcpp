# Example 12. Speed that changes during the test (Fox & Marianti, 2016; LNIRTQ of the LNIRT package),
# reproduced on the Amsterdam Chess Test data (van der Maas & Wagenmakers, 2005; 259 players, 40 items).
# Persons j, items i:
#   P(y_ij = 1) = Phi(a_i theta_j - b_i)
#   log t_ij = lambda_i - phi_i (zeta_j0 + zeta_j1 x_i + zeta_j2 x_i^2) + e_ij,  x_i = (i - 1) / K
#   speed components independent with variances s_0, s_1, s_2; ability regressed on them:
#   theta_j = mu + sum_q c_q zeta_jq + u_j, u ~ N(0, v)   (Cov(theta, zeta_q) = c_q s_q)
#   items (a, b, phi, lambda) ~ N(mu_I, Sigma_I) with E(b) = 0; identification as LNIRTQ:
#   prod(a) = 1, mean(b) = 0, prod(phi) = 1, mean(zeta_q) = 0 (all as exact group moves)
# The constant-speed model (zeta_j1 = zeta_j2 = 0) is the sampler of ex11_lnirt.R.
# Blocks: draw_probit_z() (missing responses are left untruncated); the speed vector of every person and
# the four parameters of every item as one draw_mvn() call each (batched multivariate normal draws);
# draw_invgamma(), draw_invwishart(), draw_mean_mvn(); missing log times are masked.
library(birtRcpp)

speed_sampler <- function(Y, LT, n_iter = 10000, n_burn = 1000) {
  N <- nrow(Y); K <- ncol(Y); x <- (seq_len(K) - 1) / K; Xn <- cbind(1, x, x^2)
  M <- !is.na(LT); LT0 <- ifelse(M, LT, 0)
  th <- rnorm(N); Z3 <- matrix(rnorm(3 * N, 0, 0.1), N); a <- rep(1, K); b <- rep(0, K); phi <- rep(1, K)
  lam <- colSums(LT0) / colSums(M); s2 <- rep(0.3, K); s <- rep(0.1, 3); cc <- rep(0, 3); v <- 1; mu <- 0
  V0 <- diag(c(0.05, 0.3, 0.05, 0.3)); mu0 <- c(1, 0, 1, 0); muI <- c(1, 0, 1, mean(lam)); SI <- diag(4)
  out <- matrix(NA, n_iter, 11, dimnames = list(NULL, c("var_theta", "var_zeta0", "var_zeta1", "var_zeta2", "cov_theta_zeta0",
           "cov_theta_zeta1", "cov_theta_zeta2", "cor_a_b", "cor_b_lambda", "cor_phi_lambda", "mu_theta")))
  items <- array(NA, c(n_iter, K, 4))
  row <- function(v) matrix(v, N, K, byrow = TRUE)
  for (it in seq_len(n_iter)) {
    # latent responses (missing y untruncated), ability given speed
    Z <- draw_probit_z(Y, outer(th, a) - row(b))
    g <- calc_normal_part(row(a), Z, 1, -row(b), by = "person")
    th <- draw_gauss(g$prec + 1 / v, g$num + (mu + drop(Z3 %*% cc)) / v)
    # speed vectors: prior N(0, diag(s)) and theta | zeta; log t = lambda - (phi_i Xn_i) zeta_j + e
    A <- Xn * phi; W <- M * row(1 / s2); P0 <- diag(1 / s) + outer(cc, cc) / v
    Q <- array(apply(W, 1, function(w) P0 + crossprod(A, A * w)), c(3, 3, N))
    H <- -crossprod(A, t(W * (LT0 - row(lam)))) + outer(cc, th - mu) / v
    Z3 <- t(draw_mvn(Q, H))
    # mean(zeta_q) = 0 as a group move: lambda_i - phi_i x_i^q m_q keeps log t unchanged
    m <- colMeans(Z3); Z3 <- sweep(Z3, 2, m); lam <- lam - phi * drop(Xn %*% m); mu <- mu + sum(cc * m)
    sp <- Z3 %*% t(Xn)                                                       # speed on each item, N x K
    # items: 4-variate normal per item (probit regression on (theta, -1), RT regression on (-speed, 1))
    w0 <- colSums(W); w1 <- colSums(W * sp); w2 <- colSums(W * sp^2); l0 <- colSums(W * LT0); l1 <- colSums(W * sp * LT0)
    Qi <- array(solve(SI), c(4, 4, K)); Qi[1:2, 1:2, ] <- Qi[1:2, 1:2, ] + c(sum(th^2), -sum(th), -sum(th), N)
    Qi[3, 3, ] <- Qi[3, 3, ] + w2; Qi[3, 4, ] <- Qi[3, 4, ] - w1; Qi[4, 3, ] <- Qi[4, 3, ] - w1; Qi[4, 4, ] <- Qi[4, 4, ] + w0
    X4 <- draw_mvn(Qi, drop(solve(SI, muI)) + rbind(colSums(th * Z), -colSums(Z), -l1, l0), lower = c(0, -Inf, 0, -Inf))
    ok <- !is.na(X4[1, ]); a[ok] <- X4[1, ok]; b[ok] <- X4[2, ok]; phi[ok] <- X4[3, ok]; lam[ok] <- X4[4, ok]
    # identification as group moves: prod(a) = 1 (theta * c), mean(b) = 0 (theta + d, b + a d), prod(phi) = 1 (zeta * c)
    ca <- exp(mean(log(a))); a <- a / ca; th <- th * ca; mu <- mu * ca; cc <- cc * ca
    d <- -mean(b) / mean(a); b <- b + a * d; th <- th + d; mu <- mu + d
    cp <- exp(mean(log(phi))); phi <- phi / cp; Z3 <- Z3 * cp; s <- s * cp^2; cc <- cc / cp
    sp <- Z3 %*% t(Xn)
    # residual variances: (SS + 10) / chi2, as LNIRTQ
    R <- (LT0 - row(lam) + sp * row(phi)) * M
    s2 <- draw_invgamma(colSums(M) / 2, (colSums(R^2) + 10) / 2)
    # person model: speed variances (SS + 1) / chi2, regression of theta on the speed components
    s <- draw_invgamma(rep(N / 2, 3), (colSums(Z3^2) + 1) / 2)
    D <- cbind(1, Z3); bb <- draw_mvn(crossprod(D) / v + diag(c(10, rep(1e-4, 3))), drop(crossprod(D, th)) / v)   # mu ~ N(0, .1) as LNIRTQ
    mu <- bb[1]; cc <- bb[-1]
    v <- draw_invgamma(N / 2, (sum((th - D %*% bb)^2) + 10) / 2)
    # item level: mu_I (E(b) fixed at 0) and Sigma_I
    ab <- cbind(a, b, phi, lam)
    muI <- draw_mean_mvn(ab, SI, mu0, V0); muI[2] <- 0
    SI <- draw_invwishart(4 + K, crossprod(sweep(ab, 2, muI)) + V0)
    cr <- stats::cov2cor(SI); cts <- cc * s
    out[it, ] <- c(v + sum(cc^2 * s), s, cts, cr[1, 2], cr[2, 4], cr[3, 4], mu)
    items[it, , ] <- ab
  }
  list(draws = out, items = items, n_burn = n_burn)
}

if (!requireNamespace("LNIRT", quietly = TRUE)) stop("this example uses the data of the LNIRT package")
data("AmsterdamChess", package = "LNIRT")
ac <- get("AmsterdamChess")
Y <- as.matrix(ac[paste0("Y", 1:40)]); Y[Y == 9] <- NA
RT <- as.matrix(ac[paste0("RT", 1:40)]); RT[RT == 10000] <- NA; LT <- log(RT)

set.seed(1); t0 <- Sys.time()
f <- speed_sampler(Y, LT, n_iter = 10000, n_burn = 1000)
secs <- as.numeric(Sys.time() - t0, units = "secs")
d <- f$draws[-(1:f$n_burn), ]
e <- colMeans(d)
cor_q <- e[5:7] / sqrt(e["var_theta"] * e[2:4])
cat(sprintf("%.0f s for 10000 iterations\n", secs))
cat("Variable speed (Fox & Marianti, 2016; LNIRTQ): person covariance and item correlations\n")
print(round(rbind(LNIRTQ = c(var_theta = .423, var_zeta0 = .061, var_zeta1 = .112, var_zeta2 = .062, cov0 = .115, cov1 = -.004, cov2 = -.014,
                             cor0 = .716, cor1 = -.018, cor2 = -.086, cor_a_b = .366, cor_b_lambda = .787, cor_phi_lambda = -.535),
                  blocks = c(e[1:7], cor_q, e[8:10])), 3))
cat("posterior SD (blocks):", paste(sprintf("%s %.3f", names(e)[1:7], apply(d[, 1:7], 2, sd)), collapse = ", "), "\n")
cat("published: cor(theta, zeta0) = .72 [HPD .637, .785], cor(theta, zeta1) = -.02, cor(theta, zeta2) = -.09\n")
ep <- coda::effectiveSize(f$items[-(1:f$n_burn), , 3])
cat(sprintf("smallest ESS of the time discriminations: %.0f of %d (LNIRTQ, Fox et al. 2023: 466 for phi_1)\n", min(ep), nrow(d)))

# constant speed: the LNIRT model, with the sampler of ex11_lnirt.R
src <- readLines(system.file("examples", "ex11_lnirt.R", package = "birtRcpp"))
eval(parse(text = src[1:(grep("^if \\(!requireNamespace", src) - 1)]))   # defines lnirt_sampler()
set.seed(2); fc <- lnirt_sampler(Y, LT, n_iter = 5000, n_burn = 1000)
dc <- fc$draws[-(1:1000), ]
cat("Constant speed (Fox & Marianti, 2016: var_theta about .32, var_zeta .085, rho .65; LNIRT rerun .330, .085, .639)\n")
print(round(c(var_theta = mean(dc[, "var_theta"]), var_zeta = mean(dc[, "var_zeta"]),
              rho = mean(dc[, "cov"] / sqrt(dc[, "var_theta"] * dc[, "var_zeta"]))), 3))
ess <- coda::effectiveSize(d)
cat(sprintf("smallest ESS of the person and item-level parameters: %.0f (%s) of %d draws\n", min(ess), names(which.min(ess)), nrow(d)))
