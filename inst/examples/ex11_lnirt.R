# Example 11. The joint model of LNIRT (van der Linden, 2007; Klein Entink, Fox & van der Linden, 2009;
# Fox, Klotzke & Simsek, 2023), reproducing Fox & Marianti (2017, Table 2) on the credential data.
# Persons j, items i:
#   P(y_ij = 1) = Phi(a_i theta_j - b_i)                                  (probit)
#   log t_ij = lambda_i - phi_i zeta_j + e_ij,  e_ij ~ N(0, sigma2_i)      (phi_i: time discrimination)
#   (theta_j, zeta_j) ~ N(0, Sigma_P);  (a_i, b_i, phi_i, lambda_i) ~ N(mu_I, Sigma_I)
#   identification as LNIRT: mean(theta) = mean(zeta) = 0, prod(a) = prod(phi) = 1
#   priors as LNIRT: Sigma_P ~ IW(2, I), mu_I ~ N((1, 0, 1, 0), V0), Sigma_I ~ IW(4, V0),
#   V0 = diag(.05, .3, .05, .3), sigma2_i: (SS_i + 10) / chi2
# Blocks: draw_probit_z() (Albert & Chib, 1993) makes the probit model linear, so every step is normal:
#   * (a, b) first with the latent responses integrated out, draw_items_mh(): easy items mix slowly
#     under data augmentation alone (the latent responses and (a, b) hold each other);
#   * persons: theta | zeta and zeta | theta, calc_normal_part() and draw_gauss();
#   * items: (a, b, phi, lambda) of every item jointly as a 4-variate normal, draw_mvn() with a, phi > 0;
#   * the identification as exact group moves: theta - m with b - a m, zeta - m with lambda - phi m,
#     a / c with theta * c, phi / d with zeta * d leave the likelihood unchanged;
#   * sigma2 with draw_invgamma(), Sigma_P and Sigma_I with draw_invwishart(), mu_I with draw_mean_mvn();
#   * the RT person fit l^t_j = sum_i ((log t_ij - lambda_i + phi_i zeta_j) / sigma_i)^2 of Marianti et al.
#     (2014) with calc_person_fit(); its posterior mean is compared with the chi-square 95% quantile.
library(birtRcpp)

lnirt_sampler <- function(Y, LT, n_iter = 5000, n_burn = 1000, collapse = TRUE) {
  N <- nrow(Y); K <- ncol(Y); M <- !is.na(LT); LT0 <- ifelse(M, LT, 0); nM <- colSums(M)
  th <- rnorm(N); ze <- rnorm(N); a <- rep(1, K); b <- rep(0, K); phi <- rep(1, K)
  lam <- colSums(LT0) / nM; s2 <- rep(0.3, K)
  SP <- diag(2); V0 <- diag(c(0.05, 0.3, 0.05, 0.3)); mu0 <- c(1, 0, 1, 0); muI <- c(1, 0, 1, mean(lam)); SI <- diag(4)
  keep <- (n_burn + 1):n_iter; nk <- length(keep)
  out <- matrix(NA, n_iter, 3 + 10 + 4, dimnames = list(NULL, c("var_theta", "cov", "var_zeta",
           "SI_a", "SI_b", "SI_phi", "SI_lambda", "cor_a_b", "cor_a_phi", "cor_a_lambda", "cor_b_phi", "cor_b_lambda", "cor_phi_lambda",
           "mu_a", "mu_b", "mu_phi", "mu_lambda")))
  item_draws <- array(NA, c(n_iter, K, 4)); s_ab <- matrix(0, K, 4); s_s2 <- numeric(K); s_lt <- numeric(N); s_th <- matrix(0, N, 2)
  row <- function(v) matrix(v, N, K, byrow = TRUE); n_acc <- 0
  for (it in seq_len(n_iter)) {
    # (a, b) given (phi, lambda) with the latent responses integrated out: prior N(m_i, C) from N(mu_I, Sigma_I)
    if (collapse) {
      B <- SI[1:2, 3:4] %*% solve(SI[3:4, 3:4])
      m <- sweep(t(B %*% t(sweep(cbind(phi, lam), 2, muI[3:4]))), 2, muI[1:2], "+")
      mh <- draw_items_mh(Y, th, a, b, m, solve(SI[1:2, 1:2] - B %*% SI[3:4, 1:2]), link = "probit", form = "intercept")
      a <- mh$a; b <- mh$b; n_acc <- n_acc + sum(attr(mh, "accepted"))
    }
    # latent responses, then ability given speed: Z = a theta - b + e
    Z <- draw_probit_z(Y, outer(th, a) - row(b))
    cv <- SP[1, 1] - SP[1, 2]^2 / SP[2, 2]; g <- calc_normal_part(row(a), Z, 1, -row(b), by = "person")
    th <- draw_gauss(g$prec + 1 / cv, g$num + SP[1, 2] / SP[2, 2] * ze / cv)
    m <- mean(th); th <- th - m; b <- b - a * m                         # mean(theta) = 0 as a group move
    # speed given ability: log t = lambda - phi zeta + e
    cv <- SP[2, 2] - SP[1, 2]^2 / SP[1, 1]; g <- calc_normal_part(row(-phi), LT, row(1 / s2), row(lam), by = "person")
    ze <- draw_gauss(g$prec + 1 / cv, g$num + SP[1, 2] / SP[1, 1] * th / cv)
    m <- mean(ze); ze <- ze - m; lam <- lam - phi * m                   # mean(zeta) = 0 as a group move
    # items: one 4-variate normal per item (prior N(mu_I, Sigma_I), probit regression on (theta, -1),
    # RT regression on (-zeta, 1)), all items in one draw_mvn() call
    Wk <- M / row(s2)
    w0 <- colSums(Wk); w1 <- colSums(Wk * ze); w2 <- colSums(Wk * ze^2); l0 <- colSums(Wk * LT0); l1 <- colSums(Wk * ze * LT0)
    Q <- array(solve(SI), c(4, 4, K)); Q[1:2, 1:2, ] <- Q[1:2, 1:2, ] + c(sum(th^2), -sum(th), -sum(th), N)
    Q[3, 3, ] <- Q[3, 3, ] + w2; Q[3, 4, ] <- Q[3, 4, ] - w1; Q[4, 3, ] <- Q[4, 3, ] - w1; Q[4, 4, ] <- Q[4, 4, ] + w0
    H <- drop(solve(SI, muI)) + rbind(colSums(th * Z), -colSums(Z), -l1, l0)
    x <- draw_mvn(Q, H, lower = c(0, -Inf, 0, -Inf)); ok <- !is.na(x[1, ])
    a[ok] <- x[1, ok]; b[ok] <- x[2, ok]; phi[ok] <- x[3, ok]; lam[ok] <- x[4, ok]
    # identification prod(a) = prod(phi) = 1 as exact group moves
    ca <- exp(mean(log(a))); a <- a / ca; th <- th * ca
    cp <- exp(mean(log(phi))); phi <- phi / cp; ze <- ze * cp
    # residual variances, person and item covariances, item means
    fitted <- row(lam) - outer(ze, phi)
    s2 <- draw_invgamma(nM / 2, (colSums((LT0 - fitted)^2 * M) + 10) / 2)
    SP <- draw_invwishart(2 + N, crossprod(cbind(th, ze)) + diag(2))
    ab <- cbind(a, b, phi, lam)
    muI <- draw_mean_mvn(ab, SI, mu0, V0)
    SI <- draw_invwishart(4 + K, crossprod(sweep(ab, 2, muI)) + V0)
    cr <- stats::cov2cor(SI)
    out[it, ] <- c(SP[1, 1], SP[1, 2], SP[2, 2], diag(SI), cr[1, 2], cr[1, 3], cr[1, 4], cr[2, 3], cr[2, 4], cr[3, 4], muI)
    item_draws[it, , ] <- ab
    if (it > n_burn) { s_ab <- s_ab + ab; s_s2 <- s_s2 + s2; s_lt <- s_lt + calc_person_fit(LT, fitted, s2)$lt; s_th <- s_th + cbind(th, ze) }
  }
  list(draws = out, items = item_draws, item_mean = s_ab / nk, sigma2 = s_s2 / nk, lt = s_lt / nk, person = s_th / nk, n_obs = rowSums(M),
       n_burn = n_burn, accept = n_acc / (n_iter * K))
}

if (!requireNamespace("LNIRT", quietly = TRUE)) stop("this example uses the data of the LNIRT package")
data("CredentialForm1", package = "LNIRT")
cf <- get("CredentialForm1")
s <- cf$Pretest == 6                                                         # Fox & Marianti (2017): N = 723
Y <- as.matrix(cf[s, paste0("iraw.", 1:170)])
RT <- as.matrix(cf[s, paste0("idur.", 1:170)]); RT[RT == 0] <- NA; LT <- log(RT)

set.seed(1); t0 <- Sys.time()
f <- lnirt_sampler(Y, LT, n_iter = 5000, n_burn = 1000)
secs <- as.numeric(Sys.time() - t0, units = "secs")
keep <- (f$n_burn + 1):nrow(f$draws); d <- f$draws[keep, ]
pub <- c(var_theta = .093, cov = .022, var_zeta = .022, SI_a = .287, SI_b = .24, SI_phi = .115, SI_lambda = .103,
         cor_a_b = -.236, cor_a_phi = .501, cor_a_lambda = .076, cor_b_phi = -.421, cor_b_lambda = .464, cor_phi_lambda = -.322)
est <- colMeans(d)[names(pub)]
cat(sprintf("%.0f s for 5000 iterations\n", secs))
cat("Fox & Marianti (2017, Table 2) and this sampler (posterior means; rho = cov / sqrt(var_theta var_zeta))\n")
print(round(rbind(published = c(pub, rho = .486), blocks = c(est, rho = mean(d[, "cov"] / sqrt(d[, "var_theta"] * d[, "var_zeta"])))), 3))
cat("ranges of b and lambda (published: b -1.71 to .79, lambda 2.89 to 4.85)\n")
print(round(rbind(b = range(f$item_mean[, 2]), lambda = range(f$item_mean[, 4])), 2))
lt_crit <- qchisq(0.95, f$n_obs)                                           # chi-square with the number of observed times
cat(sprintf("RT person fit: l^t above its 95%% chi-square quantile for %d of %d test takers (%.1f%%; published 124, 17.15%%)\n",
            sum(f$lt > lt_crit), length(f$lt), 100 * mean(f$lt > lt_crit)))
ess <- coda::effectiveSize(f$items[keep, , 1]); essb <- coda::effectiveSize(f$items[keep, , 2])
cat(sprintf("smallest ESS of a and b over the items: %.0f and %.0f (of %d draws)\n", min(ess), min(essb), length(keep)))
