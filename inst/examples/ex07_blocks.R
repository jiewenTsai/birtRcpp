# Example 7. Writing a sampler in R from the building blocks (?blocks)
# The rtirt_cross() sampler (normal log RT) of gibbs(), step by step in R. Every draw runs in C++;
# the R code only computes the conditional means and precisions, so new schemes (other blockings,
# cut or semi-modular inference, nested chains) can be written the same way.
# Persons j, items i:  logit P(y_ij = 1) = a_i (theta_j - b_i),
#   log t_ij = lambda_i - zeta_j - rho_i theta_j + e_ij,  e_ij ~ N(0, sigma2_i),  zeta_j ~ N(0, s).
# The same model by ECM needs no code: ecm(m) (maximum marginal likelihood), or block_model() + ecm()
# (tutorial birtRcpp_blocks, section 6). Here ECM comes first, then the hand-written Gibbs sampler.
library(birtRcpp)

cross_sampler <- function(model, n_iter = 2000, seed = 1, priors = rtirt_priors()) {
  D <- model$data; Y <- D$Y; L <- D$log_t; N <- nrow(Y); K <- ncol(Y)
  # the start of gibbs(init = "random"): the same random numbers in the same order
  set.seed(seed)
  th <- rnorm(N); a <- rep(1, K); b <- rnorm(K, 0, 0.1); ze <- 0.3 * rnorm(N)
  lam <- unname(colMeans(L)); s2 <- unname(apply(L, 2, var)); s2 <- pmax(s2, 0.1 * median(s2), 1e-4)  # variance floor
  rho <- rnorm(K, 0, 0.1); s <- 1
  # priors: lambda ~ N(mean, sd^2) above lower; rho ~ N(0, vr)
  pl <- priors$lambda; mu_l <- pl[["mean"]]; sd_l <- pl[["sd"]]; lo_l <- pl[["lower"]]; vr <- priors$rho[["sd"]]^2
  W <- matrix(1 / s2, N, K, byrow = TRUE); Z0 <- matrix(0, N, K)
  out <- matrix(NA, n_iter, 5 * K + 1)
  for (it in seq_len(n_iter)) {
    # cross-relations rho_k: y = lambda - zeta - log T = rho_k theta + e
    yr <- matrix(lam, N, K, byrow = TRUE) - ze - L
    rho <- draw_gauss(1 / vr + colSums(W * th^2), colSums(W * th * yr))
    # speed variance
    s <- draw_var(N / 2, sum(ze^2) / 2, s, priors$var_speed)
    # accuracy: Polya-Gamma variables, then b and a
    om <- draw_omega(a, b, th); ab <- draw_items_pg(Y, om, th, a, priors); a <- ab$a; b <- ab$b
    # ability: prior N(0, 1) + accuracy part + RT part
    ip <- irt_part(Y, om, a, b)
    yt <- matrix(lam, N, K, byrow = TRUE) - L - ze
    th <- draw_gauss(ip$prec + 1 + drop(W %*% rho^2), ip$num + drop((W * yt) %*% rho))
    # time intensities and residual variances (m = -zeta - rho theta)
    m <- -outer(th, rho) - ze
    lam <- draw_lambda(L, m, Z0, W, mu_l, sd_l, lo_l)
    s2 <- draw_sigma2(L, m, lam, s2, priors); W <- matrix(1 / s2, N, K, byrow = TRUE)
    # speed: prior N(0, s) + RT part
    yz <- matrix(lam, N, K, byrow = TRUE) - L - outer(th, rho)
    ze <- draw_gauss(1 / s + rowSums(W), rowSums(W * yz))
    # exact group moves: location (lambda + c, zeta + c), cross-relation (rho - d, zeta + d theta)
    sh <- draw_shift(list(lambda = lam, zeta = ze), prior = list(lambda = c(mu_l, sd_l), zeta = c(0, sqrt(s))), lower = list(lambda = lo_l))
    lam <- sh$lambda; ze <- sh$zeta
    d <- draw_gauss(K / vr + sum(th^2) / s, sum(rho) / vr - sum(ze * th) / s); rho <- rho - d; ze <- ze + d * th
    out[it, ] <- c(a, b, lam, s2, rho, s)
  }
  colnames(out) <- c(sprintf("a[%s]", D$items), sprintf("b[%s]", D$items), sprintf("lambda[%s]", D$items),
                     sprintf("sigma2t[%s]", D$items), sprintf("rho[%s]", D$items), "var_speed")
  out
}

set.seed(7)
cond <- set_cond(n_subj = 400, n_item = 8)
m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
e <- ecm(m)                                                        # ECM: estimates and standard errors
e
g <- gibbs(m, n_iter = 2000, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
cpp <- as.matrix(g$post$draws)
r <- cross_sampler(m, n_iter = nrow(cpp), seed = 1)                # the same chain, written in R
max(abs(r - cpp))                                                  # equal to rounding at every iteration
keep <- seq(nrow(cpp) / 2 + 1, nrow(cpp))
k <- c(1:3, 33:35, 41)                                             # a[1:3], rho[1:3], var_speed
round(cbind(R = colMeans(r[keep, ]), Cpp = colMeans(cpp[keep, ]), ECM = estimates(e)$est[match(colnames(cpp), estimates(e)$parameter)])[k, ], 3)
