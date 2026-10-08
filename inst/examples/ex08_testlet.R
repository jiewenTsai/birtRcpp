# Example 8. The testlet model from the building blocks (?blocks)
# Bradlow, Wainer & Wang (1999); logistic version: Wang, Bradlow & Wainer (2002).
#   logit P(y_ij = 1) = a_i (theta_j - b_i - gamma_j,d(i)),  theta_j ~ N(0, 1),  gamma_jd ~ N(0, s_d)
#   (person j, item i)
# Items of one testlet d share the person-by-testlet effect gamma_jd. Within testlet d the item
# response is a 2PL in e_jd = theta_j - gamma_jd, so the 2PL blocks apply testlet by testlet:
#   * draw_omega(), draw_items_pg() with "ability" e_d;
#   * irt_part() gives the accuracy part of e_d as exp(num e - prec e^2 / 2); written in theta it has
#     precision prec and linear term num + prec gamma_d, written in gamma_d precision prec and
#     linear term prec theta - num;
#   * draw_gauss() for theta and gamma, draw_var() for the testlet variances (half-t prior).
# testlet = FALSE fixes gamma = 0: the usual 2PL, which ignores the local dependence.
library(birtRcpp)

testlet_sampler <- function(Y, testlet, n_iter = 2000, testlet_model = TRUE, priors = rtirt_priors()) {
  N <- nrow(Y); K <- ncol(Y); D <- max(testlet)
  th <- rnorm(N); ga <- matrix(0, N, D); a <- rep(1, K); b <- rep(0, K); s <- rep(0.5, D)
  out <- matrix(NA, n_iter, 2 * K + D, dimnames = list(NULL, c(sprintf("a[%d]", 1:K), sprintf("b[%d]", 1:K), sprintf("s_testlet[%d]", 1:D))))
  th_sum <- th_ss <- numeric(N)
  for (it in seq_len(n_iter)) {
    P <- H <- matrix(0, N, D)                          # accuracy part of e_d: precision and linear term
    for (d in 1:D) {
      j <- testlet == d; e <- th - ga[, d]
      om <- draw_omega(a[j], b[j], e)
      ab <- draw_items_pg(Y[, j, drop = FALSE], om, e, a[j], priors); a[j] <- ab$a; b[j] <- ab$b
      ip <- irt_part(Y[, j, drop = FALSE], om, a[j], b[j]); P[, d] <- ip$prec; H[, d] <- ip$num
    }
    th <- draw_gauss(1 + rowSums(P), rowSums(H + P * ga))                 # prior N(0, 1)
    if (testlet_model) {
      for (d in 1:D) ga[, d] <- draw_gauss(1 / s[d] + P[, d], P[, d] * th - H[, d])
      for (d in 1:D) s[d] <- draw_var(N / 2, sum(ga[, d]^2) / 2, s[d], priors$var_speed)
    }
    out[it, ] <- c(a, b, if (testlet_model) s else rep(0, D))
    if (it > n_iter / 2) { th_sum <- th_sum + th; th_ss <- th_ss + th^2 }
  }
  n <- n_iter / 2
  list(draws = out, theta_mean = th_sum / n, theta_sd = sqrt(th_ss / n - (th_sum / n)^2))
}

set.seed(8)
N <- 1000; testlet <- rep(1:6, each = 4); K <- length(testlet)
s_true <- c(0.2, 0.5, 1, 0.2, 0.5, 1)                  # testlet variances
a_true <- runif(K, 0.7, 2); b_true <- rnorm(K)
theta <- rnorm(N); gamma <- sapply(s_true, function(v) rnorm(N, 0, sqrt(v)))
Y <- matrix(rbinom(N * K, 1, plogis(sweep(theta - gamma[, testlet], 2, b_true) * rep(a_true, each = N))), N)

set.seed(1); ft <- testlet_sampler(Y, testlet)
set.seed(1); f2 <- testlet_sampler(Y, testlet, testlet_model = FALSE)
keep <- 1001:2000
cat("testlet variances\n")
print(round(rbind(true = s_true, est = colMeans(ft$draws[keep, 2 * K + 1:6])), 2))
cat("item parameters: correlation with the true values and largest error\n")
err <- function(f, k) c(cor = cor(colMeans(f$draws[keep, k]), c(a_true, b_true)[k]), max_abs_err = max(abs(colMeans(f$draws[keep, k]) - c(a_true, b_true)[k])))
print(round(rbind(`testlet a` = err(ft, 1:K), `2PL a` = err(f2, 1:K), `testlet b` = err(ft, K + 1:K), `2PL b` = err(f2, K + 1:K)), 3))
cat("ability: mean posterior SD, and RMSE of the posterior mean\n")
print(round(rbind(testlet = c(post_sd = mean(ft$theta_sd), rmse = sqrt(mean((ft$theta_mean - theta)^2))),
                  `2PL` = c(post_sd = mean(f2$theta_sd), rmse = sqrt(mean((f2$theta_mean - theta)^2)))), 3))
