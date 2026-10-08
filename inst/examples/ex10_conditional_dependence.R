# Example 10. Conditional dependence between response time and accuracy (Bolsinova, De Boeck & Tijmstra, 2017)
# Persons j, items i. The residual log time of a response enters the logit of the same response:
#   logit P(y_ij = 1) = a_i (theta_j - b_i) + phi_i r_ij,   r_ij = log t_ij - lambda_i + zeta_j
#   log t_ij = lambda_i - zeta_j + e_ij,  e_ij ~ N(0, sigma2_i),  theta ~ N(0, 1),  zeta | theta ~ N(c theta, v)
# r_ij > 0: slower than expected for this person and item. phi_i > 0: slow responses are more often
# correct. This is the intercept part of their model (they also let the slope depend on the residual,
# and standardize it by sigma_i; the unstandardized residual keeps every full conditional conjugate).
# Given the Polya-Gamma variables the logit is linear in a, b, phi, theta, zeta and lambda, so every
# step is a normal (or truncated normal) draw, now with an accuracy part also for zeta and lambda:
#   * draw_pg(): omega for the matrix of logits, pseudo-observations z = (y - 1/2) / omega;
#   * calc_normal_part(): precision and linear term of each parameter, from the pseudo-observations
#     (weights omega) and from the log times (weights 1 / sigma2);
#   * draw_tnorm() for a, draw_gauss() for b, phi, theta, zeta, lambda; draw_sigma2(); draw_var();
#   * draw_shift() still applies: (lambda + delta, zeta + delta) changes neither log t nor r.
# phi_model = FALSE fixes phi = 0 (the usual hierarchical model with conditional independence).
library(birtRcpp)

cd_sampler <- function(Y, L, n_iter = 3000, phi_model = TRUE, priors = rtirt_priors()) {
  N <- nrow(Y); K <- ncol(Y)
  th <- rnorm(N); ze <- rep(0, N); a <- rep(1, K); b <- rep(0, K); phi <- rep(0, K)
  lam <- colMeans(L); s2 <- apply(L, 2, var); cc <- 0; v <- 1
  pa <- priors$a; pb <- priors$b; pl <- priors$lambda; sphi <- 1                 # phi_i ~ N(0, 1)
  out <- matrix(NA, n_iter, 3 * K + 2, dimnames = list(NULL, c(sprintf("a[%d]", 1:K), sprintf("b[%d]", 1:K), sprintf("phi[%d]", 1:K), "c", "v")))
  row <- function(v) matrix(v, N, K, byrow = TRUE)
  for (it in seq_len(n_iter)) {
    R <- L - row(lam) + ze                                                      # residual log times
    lin <- function() sweep(outer(th, a), 2, a * b) + R * row(phi)
    eta <- lin(); om <- draw_pg(eta); z <- (Y - 0.5) / om
    # item parameters: b (x = -a), a (x = theta - b, truncated), phi (x = r)
    g <- calc_normal_part(row(-a), z, om, eta + row(a * b))
    b <- draw_gauss(1 / pb[["sd"]]^2 + g$prec, pb[["mean"]] / pb[["sd"]]^2 + g$num); eta <- lin()
    x <- outer(th, b, "-"); g <- calc_normal_part(x, z, om, eta - x * row(a))
    p <- 1 / pa[2]^2 + g$prec; a <- draw_tnorm((pa[1] / pa[2]^2 + g$num) / p, 1 / sqrt(p), 0, Inf); eta <- lin()
    if (phi_model) { g <- calc_normal_part(R, z, om, eta - R * row(phi)); phi <- draw_gauss(1 / sphi^2 + g$prec, g$num); eta <- lin() }
    # theta: accuracy part (x = a_i) + speed regression zeta | theta ~ N(c theta, v)
    g <- calc_normal_part(row(a), z, om, eta - outer(th, a), by = "person")
    th <- draw_gauss(1 + g$prec + cc^2 / v, g$num + cc * ze / v); eta <- lin()
    # zeta: RT part (x = -1, log t = lambda - zeta) + accuracy part (x = phi_i) + prior N(c theta, v)
    g1 <- calc_normal_part(-1, L, row(1 / s2), row(lam), by = "person")
    g2 <- calc_normal_part(row(phi), z, om, eta - outer(ze, phi), by = "person")
    ze <- draw_gauss(1 / v + g1$prec + g2$prec, cc * th / v + g1$num + g2$num)
    # lambda: RT part (x = 1) + accuracy part (x = -phi_i) + prior
    R <- L - row(lam) + ze; eta <- lin()
    g1 <- calc_normal_part(1, L, row(1 / s2), -ze %o% rep(1, K))
    g2 <- calc_normal_part(row(-phi), z, om, eta + row(lam * phi))
    lam <- draw_gauss(1 / pl[["sd"]]^2 + g1$prec + g2$prec, pl[["mean"]] / pl[["sd"]]^2 + g1$num + g2$num)
    s2 <- draw_sigma2(L, -ze %o% rep(1, K), lam, s2, priors)
    sh <- draw_shift(list(lambda = lam, zeta = ze), prior = list(lambda = c(pl[["mean"]], pl[["sd"]]), zeta = list(mean = cc * th, sd = sqrt(v))),
                     lower = list(lambda = pl[["lower"]])); lam <- sh$lambda; ze <- sh$zeta
    # speed regression on ability
    cc <- draw_gauss(1 / priors$cov[["sd"]]^2 + sum(th^2) / v, sum(th * ze) / v)
    v <- draw_var(N / 2, sum((ze - cc * th)^2) / 2, v, priors$var_speed)
    out[it, ] <- c(a, b, phi, cc, v)
  }
  out
}

# simulation in the direction Bolsinova et al. found: slow responses help on hard items (phi > 0) and
# hurt on easy items (phi < 0)
set.seed(10)
N <- 1000; K <- 12
b_true <- seq(-1.5, 1.5, length.out = K); a_true <- runif(K, 0.8, 1.8); phi_true <- 0.5 * b_true
lam_true <- rnorm(K, 4, 0.3); sig_true <- runif(K, 0.4, 0.6)
theta <- rnorm(N); zeta <- -0.4 * theta + rnorm(N, 0, 0.5)
E <- matrix(rnorm(N * K), N) * rep(sig_true, each = N)
L <- matrix(lam_true, N, K, byrow = TRUE) - zeta + E
Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(theta, a_true), 2, a_true * b_true) + E * rep(phi_true, each = N))), N)

set.seed(1); fd <- cd_sampler(Y, L)
keep <- 1501:3000; ph <- fd[keep, 2 * K + 1:K]
cat("phi: true value, posterior mean and 95% interval\n")
print(round(cbind(b = b_true, phi_true, est = colMeans(ph), lower = apply(ph, 2, quantile, 0.025), upper = apply(ph, 2, quantile, 0.975)), 2))
cat("the pattern of Bolsinova et al.: phi against the estimated difficulty\n")
b_hat <- colMeans(fd[keep, K + 1:K])
print(round(coef(lm(colMeans(ph) ~ b_hat)), 2))                     # true slope 0.5
