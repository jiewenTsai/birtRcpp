# Example 9. Compensatory multidimensional IRT with response times (Man, Harring, Jiao & Zhan, 2019)
#   logit P(y_ij = 1) = a_i1 theta_j1 + a_i2 theta_j2 + d_i        (a_im = 0 where item i does not load)
#   log t_ij = lambda_i - zeta_j + e_ij,  e_ij ~ N(0, sigma2_i)        (time discrimination 1 / sigma_i; person j, item i)
#   persons: (theta_1, theta_2, zeta) trivariate normal, Var(theta_m) = 1, written as
#            theta ~ N(0, R), R = [1 r; r 1],  zeta | theta ~ N(c1 theta_1 + c2 theta_2, v)
#   items:   (d_i, lambda_i) bivariate normal with mean mu_I and covariance S_I ~ inverse Wishart
# Steps (all Gibbs except the correlation r):
#   * draw_pg(): Polya-Gamma variables for the matrix of logits.
#   * Given omega the logit is linear in every parameter: normal pseudo-observations z = (y - 1/2) / omega
#     with precisions omega. calc_normal_part() gives the precision and linear term of any parameter
#     entering eta linearly; draw_tnorm() draws the loadings (a >= 0), draw_gauss() the intercepts and
#     the abilities (one dimension at a time).
#   * Speed, time intensities: calc_normal_part() and draw_gauss(); residual variances: draw_sigma2();
#     speed regression (c1, c2) with draw_mvn(), v with draw_var(); r with draw_cor().
#   * Item level (d, lambda): conditional normal priors given the other one; (mu_I, Sigma_I) with
#     draw_mean_mvn() and draw_invwishart().
# rt = FALSE drops the response times (the multidimensional IRT model alone).
library(birtRcpp)

mirt_rt_sampler <- function(Y, L, Q, n_iter = 3000, rt = TRUE, priors = rtirt_priors()) {
  N <- nrow(Y); K <- ncol(Y)
  th <- matrix(rnorm(2 * N), N); ze <- rep(0, N); A <- Q * 1; d <- rep(0, K); r <- 0; cc <- c(0, 0); v <- 1
  lam <- if (rt) colMeans(L) else rep(0, K); s2 <- rep(0.3, K)
  muI <- c(0, mean(lam)); SI <- diag(c(1, 1)); pa <- priors$a
  out <- matrix(NA, n_iter, 3 * K + 6, dimnames = list(NULL, c(sprintf("a1[%d]", 1:K), sprintf("a2[%d]", 1:K), sprintf("d[%d]", 1:K),
                                                             "r", "cor_th1_zeta", "cor_th2_zeta", "var_zeta", "cor_d_lambda", "sd_d")))
  th_sum <- th_ss <- matrix(0, N, 2)
  row <- function(v) matrix(v, N, K, byrow = TRUE)
  for (it in seq_len(n_iter)) {
    eta <- th %*% t(A) + row(d)
    om <- draw_pg(eta); z <- (Y - 0.5) / om                                # Polya-Gamma pseudo-observations
    # loadings a_im, one dimension at a time: x = theta_jm (truncated normal, prior N+)
    for (m in 1:2) {
      i <- Q[, m] == 1; x <- matrix(th[, m], N, K)
      g <- calc_normal_part(x, z, om, eta - x * row(A[, m]))
      p <- 1 / pa[2]^2 + g$prec[i]
      new <- draw_tnorm((pa[1] / pa[2]^2 + g$num[i]) / p, 1 / sqrt(p), 0, Inf)
      eta[, i] <- eta[, i] + outer(th[, m], new - A[i, m]); A[i, m] <- new
    }
    # intercepts d_i: x = 1, prior from the item level d_i | lambda_i
    md <- if (rt) muI[1] + SI[1, 2] / SI[2, 2] * (lam - muI[2]) else rep(muI[1], K)
    vd <- if (rt) SI[1, 1] - SI[1, 2]^2 / SI[2, 2] else SI[1, 1]
    g <- calc_normal_part(1, z, om, eta - row(d))
    new <- draw_gauss(1 / vd + g$prec, md / vd + g$num); eta <- eta + row(new - d); d <- new
    # abilities, one dimension at a time: x = a_im, prior theta_jm | theta_jo, plus zeta | theta
    for (m in 1:2) {
      o <- 3 - m; g <- calc_normal_part(row(A[, m]), z, om, eta - outer(th[, m], A[, m]), by = "person")
      p <- 1 / (1 - r^2) + g$prec + if (rt) cc[m]^2 / v else 0
      h <- r * th[, o] / (1 - r^2) + g$num + if (rt) cc[m] * (ze - cc[o] * th[, o]) / v else 0
      new <- draw_gauss(p, h); eta <- eta + outer(new - th[, m], A[, m]); th[, m] <- new
    }
    if (rt) {
      # speed: log t = lambda - zeta + e, so x = -1 with w = 1 / sigma2; prior zeta | theta
      mz <- drop(th %*% cc); g <- calc_normal_part(-1, L, row(1 / s2), row(lam), by = "person")
      ze <- draw_gauss(1 / v + g$prec, mz / v + g$num)
      # time intensities: x = 1, prior lambda_i | d_i from the item level
      ml <- muI[2] + SI[1, 2] / SI[1, 1] * (d - muI[1]); vl <- SI[2, 2] - SI[1, 2]^2 / SI[1, 1]
      g <- calc_normal_part(1, L, row(1 / s2), -ze %o% rep(1, K))
      lam <- draw_gauss(1 / vl + g$prec, ml / vl + g$num)
      s2 <- draw_sigma2(L, -ze %o% rep(1, K), lam, s2, priors)
      # location move (lambda + delta, zeta + delta): exact, with the item-level prior
      dl <- draw_gauss(K / vl + N / v, sum(ml - lam) / vl + sum(mz - ze) / v); lam <- lam + dl; ze <- ze + dl
      # regression of speed on the abilities, then its residual variance
      cc <- draw_mvn(diag(2) / priors$cov[["sd"]]^2 + crossprod(th) / v, drop(crossprod(th, ze)) / v)
      v <- draw_var(N / 2, sum((ze - th %*% cc)^2) / 2, v, priors$var_speed)
    }
    r <- draw_cor(th[, 1], th[, 2], r)                                     # abilities with variance 1
    # item level: mu_I ~ N(0, 100 I), Sigma_I ~ inverse Wishart(4, I)
    if (rt) {
      X <- cbind(d, lam); muI <- draw_mean_mvn(X, SI, c(0, 0), diag(2) * 100)
      SI <- draw_invwishart(4 + K, diag(2) + crossprod(sweep(X, 2, muI)))
    } else {
      muI[1] <- draw_mean_mvn(cbind(d), SI[1, 1, drop = FALSE], 0, matrix(100))
      SI[1, 1] <- draw_invgamma(1 + K / 2, 1 + sum((d - muI[1])^2) / 2)
    }
    Rm <- matrix(c(1, r, r, 1), 2); s_tz <- drop(Rm %*% cc); vz <- drop(t(cc) %*% Rm %*% cc) + v
    out[it, ] <- c(A[, 1], A[, 2], d, r, s_tz / sqrt(vz), vz, SI[1, 2] / sqrt(SI[1, 1] * SI[2, 2]), sqrt(SI[1, 1]))
    if (it > n_iter / 2) { th_sum <- th_sum + th; th_ss <- th_ss + th^2 }
  }
  n <- n_iter / 2
  list(draws = out, theta_mean = th_sum / n, theta_sd = sqrt(th_ss / n - (th_sum / n)^2))
}

# the simulation design of Man et al. (2019): two dimensions, items 1-2 on dimension 1, 3-4 on
# dimension 2, the rest on both; correlation .3 between the abilities, -.3 between ability and speed
set.seed(9)
N <- 1500; K <- 15
Q <- cbind(c(1, 1, 0, 0, rep(1, 11)), c(0, 0, 1, 1, rep(1, 11)))
Sp <- matrix(c(1, .3, -.3 * .5, .3, 1, -.3 * .5, -.3 * .5, -.3 * .5, .25), 3)   # (theta1, theta2, zeta), SD(zeta) = .5
P <- matrix(rnorm(3 * N), N) %*% chol(Sp); theta <- P[, 1:2]; zeta <- P[, 3]
A_true <- Q * matrix(runif(2 * K, 0.7, 1.8), K)
SI_true <- matrix(c(1, .3 * .4, .3 * .4, .16), 2)                                # (d, lambda): SDs 1 and .4, correlation .3
IL <- matrix(rnorm(2 * K), K) %*% chol(SI_true); d_true <- IL[, 1]; lam_true <- 4 + IL[, 2]
sig_true <- 1 / runif(K, 1, 2)                                                   # time discrimination U(1, 2)
Y <- matrix(rbinom(N * K, 1, plogis(theta %*% t(A_true) + matrix(d_true, N, K, byrow = TRUE))), N)
L <- matrix(lam_true, N, K, byrow = TRUE) - zeta + matrix(rnorm(N * K), N) * rep(sig_true, each = N)

set.seed(1); fj <- mirt_rt_sampler(Y, L, Q)
set.seed(1); fi <- mirt_rt_sampler(Y, L, Q, rt = FALSE)
keep <- 1501:3000; e <- colMeans(fj$draws[keep, ])
cat("person and item-level structure (sample: the values in the simulated data; only 15 items inform the item level)\n")
print(round(rbind(true = c(r = .3, cor_th1_zeta = -.3, cor_th2_zeta = -.3, var_zeta = .25, cor_d_lambda = .3, sd_d = 1),
                  sample = c(cor(theta)[1, 2], cor(theta, zeta), var(zeta), cor(d_true, lam_true), sd(d_true)),
                  est = e[c("r", "cor_th1_zeta", "cor_th2_zeta", "var_zeta", "cor_d_lambda", "sd_d")],
                  ess = coda::effectiveSize(fj$draws[keep, c("r", "cor_th1_zeta", "cor_th2_zeta", "var_zeta", "cor_d_lambda", "sd_d")])), 3))
cat("item parameters: RMSE (loadings where the item loads)\n")
rmse <- function(f) { e <- colMeans(f$draws[keep, ]); a <- c(e[1:K], e[K + 1:K])[c(Q) == 1]
  c(a = sqrt(mean((a - A_true[Q == 1])^2)), d = sqrt(mean((e[2 * K + 1:K] - d_true)^2))) }
print(round(rbind(joint = rmse(fj), `IRT only` = rmse(fi)), 3))
cat("abilities: RMSE of the posterior means and mean posterior SD\n")
th_err <- function(f) c(rmse1 = sqrt(mean((f$theta_mean[, 1] - theta[, 1])^2)), rmse2 = sqrt(mean((f$theta_mean[, 2] - theta[, 2])^2)),
                        post_sd = mean(f$theta_sd))
print(round(rbind(joint = th_err(fj), `IRT only` = th_err(fi)), 3))
