al_tex <- function(x, ...) paste(as.character(algorithm(x, ..., format = "latex")), collapse = "\n")   # algorithm() as LaTeX (ASCII)

sim_2pl_model <- function(N = 200, K = 5, seed = 1) {
  set.seed(seed); th <- rnorm(N); a <- runif(K, 0.8, 2); b <- seq(-1.5, 1.5, length.out = K)
  Y <- matrix(rbinom(N * K, 1, plogis(sweep(outer(th, a), 2, a * b))), N)
  list(Y = Y, m = block_model(Y = Y, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
                              Y ~ bernoulli_logit(a * (theta - b))))
}

test_that("block_model() declares parameters, data and dimensions, and says what is wrong", {
  d <- sim_2pl_model()
  expect_s3_class(d$m, "block_model"); expect_named(d$m$pars, c("theta", "a", "b"))
  expect_match(paste(as.character(equations(d$m, "latex")), collapse = "\n"), "operatorname\\{logit\\} P\\(Y_\\{ji\\} = 1\\) &= a_\\{i\\}", fixed = FALSE)
  Y <- d$Y
  expect_error(block_model(Y = Y, a[item] ~ normal(0, 1), Y ~ bernoulli_logit(a * (theta - 1))), "unknown name.*theta")
  expect_error(block_model(Y = Y, theta[person] ~ normal(0, 1), a[item] ~ normal(theta, 1), Y ~ bernoulli_logit(a * theta)),
               "one value per person.*one value per item")
  expect_error(block_model(Y = Y, theta[person] ~ gamma(1, 1)), "distribution must be")
  expect_error(block_model(Y = Y, theta[person] ~ normal(0, 1), theta[person] ~ normal(0, 2)), "declared twice")
  expect_error(block_model(Y = Y, Z = matrix(0, 3, 3), theta[person] ~ normal(0, 1)), "persons x items")
  expect_error(block_model(Y = Y + 1, theta[person] ~ normal(0, 1), Y ~ bernoulli_logit(theta)), "0/1")
})

test_that("steps are checked against the model", {
  m <- sim_2pl_model()$m
  expect_error(gibbs(m, list(step_gibbs(theta, a))), "no step updates b")
  expect_error(gibbs(m, list(step_gibbs(theta, a, b), step_mala(theta, a))), "one dimension")
  expect_error(gibbs(m, list(step_gibbs(theta, a, c))), "not a parameter")
  expect_error(gibbs(m, list(step_gibbs(theta, a, b), step_shift(theta = 1, a = 1))), "changes the logistic part|not linear")
  m2 <- block_model(Y = sim_2pl_model()$Y, theta[person] ~ normal(0, 1), la[item] ~ normal(0, 1), b[item] ~ normal(0, 3),
                    Y ~ bernoulli_logit(exp(la) * (theta - b)))
  expect_error(gibbs(m2, list(step_gibbs(theta, la, b))), "not linear in la.*step_mala")
  fit <- gibbs(m2, n_iter = 50, n_chain = 1)                                   # default: Gibbs where possible, MALA otherwise
  expect_match(al_tex(fit), "\\mathrm{la}_{i} &\\leftarrow \\text{MALA on}", fixed = TRUE)
})

test_that("the Polya-Gamma variables are drawn when a step needs them, and stacked steps pass the Geweke check", {
  m <- sim_2pl_model(N = 20, K = 4, seed = 2)$m
  al <- al_tex(m, steps = list(step_mala(theta), step_gibbs(a, b)))
  expect_match(al, "\\omega_{ji} \\mid \\cdot &\\sim \\mathrm{PG}(1, a_{i}\\left(\\theta_{j} - b_{i}\\right))", fixed = TRUE)
  expect_match(al, "redrawn because MALA moved theta", fixed = TRUE)
  expect_match(al, "m_{i} &= -\\sum_{j} a_{i}\\left(\\kappa_{ji} - \\omega_{ji}\\,a_{i}\\,\\theta_{j}\\right)", fixed = TRUE)   # the linear term of b
  fit <- gibbs(m, list(step_mala(theta), step_gibbs(a, b)), n_iter = 100, n_chain = 1)
  expect_equal(unname(fit$n_omega), 100)                                                 # once per iteration (after MALA)
  skip_on_cran()
  for (st in list(list(step_mala(theta), step_gibbs(a, b)), list(step_mala(a, b), step_gibbs(theta), step_shift(theta = 1, b = 1)))) {
    r <- check_sampler(m, steps = st, n_iter = 2e4, seed = 3)
    expect_true(attr(r, "passed"))
  }
})

test_that("a declared 2PL agrees with the same sampler written by hand", {
  skip_on_cran()
  d <- sim_2pl_model(N = 400, K = 5, seed = 4)
  fit <- gibbs(d$m, list(step_mala(theta), step_gibbs(a, b), step_shift(theta = 1, b = 1)), n_iter = 4000, n_chain = 2, seed = 1)
  hand <- function(s, x) {
    s$omega <- draw_omega(s$a, s$b, s$theta)
    s[c("a", "b")] <- draw_items_pg(x$Y, s$omega, s$theta, s$a)
    s$theta <- draw(irt_part(x$Y, s$omega, s$a, s$b) + normal_prior(0, 1)); s }
  ref <- run_sampler(list(Y = d$Y), function(x) list(theta = rnorm(400, 0, 0.1), a = rep(1, 5), b = rep(0, 5)), hand, c("a", "b"),
                     n_iter = 4000, n_chain = 2, seed = 2)
  s1 <- summary(fit); s2 <- summary(ref)
  z <- (s1$mean - s2$mean) / sqrt(s1$sd^2 / s1$ess + s2$sd^2 / s2$ess)
  expect_lt(max(abs(z)), 4)
  expect_true("theta" %in% names(fit$mean))
})

test_that("normal parts, variances and shifts of a declared RT model pass the Geweke check", {
  skip_on_cran()
  set.seed(5); N <- 20; K <- 4
  L <- matrix(rnorm(N * K, 3), N); Y <- matrix(rbinom(N * K, 1, 0.5), N); L[1, 2] <- NA
  m <- block_model(Y = Y, logT = L,
    theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, sqrt(s)), s ~ half_t(3, 1),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(3, 1), rho[item] ~ normal(0, 1),
    sigma2[item] ~ inv_gamma(3, 1),
    Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  expect_match(al_tex(m), "\\sigma_{i}^2 \\mid \\cdot &\\sim \\mathrm{IG}(3 + \\frac{n_{i}}{2}", fixed = TRUE)
  expect_match(al_tex(m), "u \\mid s &\\sim \\mathrm{IG}(2, \\frac{3}{s} + 1)", fixed = TRUE)                 # the half-t auxiliary
  expect_output(print(algorithm(m)), "Gibbs: one iteration")
  st <- list(step_gibbs(rho, s, a, b, theta, lambda, sigma2, zeta), step_shift(lambda = 1, zeta = 1), step_shift(rho = -1, zeta = theta))
  r <- check_sampler(m, steps = st, n_iter = 2e4, seed = 4)
  expect_true(attr(r, "passed"))
})

test_that("ecm() runs the same steps as EM and agrees with a hand-written EM", {
  d <- sim_2pl_model(N = 300, K = 5, seed = 6)
  m <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
                   eta <- a * (theta - b), Y ~ bernoulli_logit(eta))
  expect_match(paste(as.character(equations(m, "latex")), collapse = "\n"), "a_\\{i\\}\\\\left\\(\\\\theta_\\{j\\} - b_\\{i\\}\\\\right\\)")   # the definition is substituted
  fit <- ecm(m, list(step_gibbs(theta, a, b), step_shift(theta = 1, b = 1)), prior = "ml")   # theta integrated, shift skipped
  st <- function(s, x) { e <- estep_irt(x$Y, s$a, s$b); s[c("a", "b")] <- mstep_items(e, s$a, s$b); s$.loglik <- e$loglik; s }
  ml <- run_em(list(Y = d$Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), st, loglik = function(s, x) estep_irt(x$Y, s$a, s$b)$loglik, positive = "a")
  expect_lt(max(abs(c(fit$state$a, fit$state$b) - c(ml$state$a, ml$state$b))), 1e-4)
  expect_lt(max(abs(fit$estimates$se / ml$estimates$se - 1)), 1e-3)
  expect_named(fit$latent, "theta"); expect_length(fit$latent$theta$mean, 300)
  al <- al_tex(fit)
  expect_match(al, "\\bar{\\omega}_{ji} &= \\frac{\\operatorname{tanh}", fixed = TRUE)
  expect_match(al, "p_{i} &= \\sum_{j} \\mathbb{E}\\left[\\bar{\\omega}_{ji}\\,\\left(\\theta_{j} - b_{i}\\right)^2\\right]", fixed = TRUE)
  expect_match(al, "conditional mode", fixed = TRUE)
  fg <- ecm(m, list(step_mala(a, b)), prior = "ml")                                 # gradient steps reach the same maximum
  expect_lt(max(abs(c(fg$state$a, fg$state$b) - c(ml$state$a, ml$state$b))), 1e-3)
  expect_error(ecm(m, list(step_gibbs(a))), "no step updates b")
})

test_that("check_ecm() passes a right EM and fails a wrong one", {
  d <- sim_2pl_model(N = 200, K = 4, seed = 7)
  r <- check_ecm(d$m, list(step_gibbs(a, b)), n_start = 2, seed = 1)
  expect_s3_class(r, "em_check"); expect_true(attr(r, "passed"))
  set.seed(7); N <- 400; K <- 6; Y <- matrix(rbinom(N * K, 1, plogis(outer(rnorm(N, 0, 1.4), rep(1, K)) - rep(seq(-1, 1, length.out = K), each = N))), N)
  ll <- function(s, d) estep_irt(d$Y, rep(1, K), s$b, sd = s$sigma)$loglik
  mk <- function(sig) function(s, d) { e <- estep_irt(d$Y, rep(1, K), s$b, sd = s$sigma); s$.loglik <- e$loglik
    s$b <- mstep_items(e, rep(1, K), s$b, one_pl = TRUE)$b; s$sigma <- sig(e); s }
  good <- check_ecm(data = list(Y = Y), init = list(b = rep(0, K), sigma = 1), step = mk(function(e) sqrt(mean(e$eap^2 + e$psd^2))),
                   loglik = ll, positive = "sigma", n_start = 3, seed = 2)
  expect_true(attr(good, "passed"))
  bad <- check_ecm(data = list(Y = Y), init = list(b = rep(0, K), sigma = 1), step = mk(function(e) sd(e$eap)),
                  loglik = ll, positive = "sigma", n_start = 3, seed = 2)
  expect_false(attr(bad, "passed"))
})

test_that("ecm() integrates a normal-linear latent exactly and reproduces ecm() for the cross model", {
  skip_on_cran()
  set.seed(3); cond <- set_cond(n_subj = 200, n_item = 6); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  ref <- ecm(rtirt_cross(sim)); er <- estimates(ref)
  m <- block_model(Y = unname(sim$Y), logT = unname(sim$log_t),
    theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, sqrt(v)), v ~ half_t(3, 1),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
    sigma2[item] ~ half_t(3, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  f <- ecm(m, prior = "ml", init = list(lambda = colMeans(sim$log_t)))
  expect_match(al_tex(f), "integrated exactly at each node", fixed = TRUE)
  key <- sub("sigma2t", "sigma2", sub("var_speed", "v", sub("\\[item0?([0-9]+)\\]", "[\\1]", er$parameter)))
  i <- match(key, f$estimates$parameter)
  expect_lt(abs(f$loglik - as.numeric(logLik(ref))), 1e-3)
  expect_lt(max(abs(f$estimates$est[i] - er$est)), 1e-3)
  expect_lt(max(abs(f$estimates$se[i] / er$se - 1)), 0.01)
  fc <- ecm(m, prior = "ml", init = list(lambda = colMeans(sim$log_t)), quadrature = "ba81")    # same nodes for everyone: Y collapsed
  expect_match(al_tex(fc), "summed over persons at each node", fixed = TRUE)
  expect_lt(max(abs(fc$estimates$est[i] - er$est)), 1e-3)
  expect_lt(max(abs(fc$estimates$se[i] / er$se - 1)), 0.01)
})

test_that("the M-step of an exact latent variable through 2L points equals the m -/+ sd points per node", {
  set.seed(4); cond <- set_cond(n_subj = 150, n_item = 5); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  LT <- unname(sim$log_t); LT[sample(length(LT), 30)] <- NA
  m <- block_model(Y = unname(sim$Y), logT = LT, theta[person] ~ normal(0, 1), zeta[person] ~ normal(r * theta, 1), r ~ normal(0, 1),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), phi[item] ~ normal(1, 1, lower = 0),
    sigma2[item] ~ half_t(3, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - phi * zeta, sqrt(sigma2)))
  E <- bm_setup(m, NULL, NULL, NULL, "map", list(lambda = colMeans(LT, na.rm = TRUE)))
  expect_true(E$sigma)
  E2 <- E; E2$sigma <- FALSE; E2$ns <- 2L; E2$Ms <- bm_stackM(m, 2L * E$ng)              # h/p -/+ 1/sqrt(p) at every node
  E3 <- E2; E3$affine <- FALSE                                                            # and the general exact step, node by node
  s <- E$s; C <- bm_centres(E, s, NULL); e <- bm_estep(E, s, C); C <- bm_centres(E, s, e)
  e1 <- bm_estep(E, s, C)
  for (Ek in list(E2, E3)) {
    e2 <- bm_estep(Ek, s, C)
    expect_equal(e1$obj, e2$obj, tolerance = 1e-10)
    expect_equal(bm_qvalue(E, s, e1), bm_qvalue(Ek, s, e2), tolerance = 1e-10)
    s1 <- bm_step(E, C)(s, NULL); s2 <- bm_step(Ek, C)(s, NULL)
    for (v in E$fixed) expect_equal(s1[[v]], s2[[v]], tolerance = 1e-8)
    expect_equal(unlist(bm_score(E, C)(s, NULL)), unlist(bm_score(Ek, C)(s, NULL)), tolerance = 1e-8)
  }
})

test_that("Newton and quasi-Newton (nlminb) reach the maximum of ECM", {
  d <- sim_2pl_model(N = 300, K = 5, seed = 9)
  m <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), la[item] ~ normal(0, 1), b[item] ~ normal(0, 3), Y ~ bernoulli_logit(exp(la) * (theta - b)))
  e <- ecm(m, list(step_mala(la), step_gibbs(b)))                                    # la has no closed-form step
  for (me in c("newton", "quasi-newton")) {
    f <- mml(m, method = me)
    expect_true(f$converged)
    expect_lt(abs(f$loglik - e$loglik), 1e-5)
    expect_lt(max(abs(f$estimates$est - e$estimates$est)), 1e-3)
    expect_lt(max(abs(f$estimates$se / e$estimates$se - 1)), 1e-3)
  }
  expect_match(al_tex(f), "nlminb", fixed = TRUE)
})

test_that("the Hessian by Louis' formula equals the numerical derivative of the score", {
  d <- sim_2pl_model(N = 300, K = 5, seed = 8)
  m <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(mu, sqrt(tau)),
                   mu ~ normal(0, 1), tau ~ inv_gamma(2, 1), Y ~ bernoulli_logit(a * (theta - b)))
  for (q in c("adaptive", "ba81")) {
    fit <- ecm(m, prior = "map", se = FALSE, quadrature = q)
    E <- bm_setup(m, NULL, NULL, NULL, "map", NULL, quadrature = q); s <- fit$state[E$fixed]; C <- fit$ecm$C
    sc <- bm_score(E, C); g <- function(x) unlist(sc(utils::relist(x, s), NULL)[E$fixed])
    Hn <- calc_hessian(g, unlist(s), h = 1e-5); Hl <- bm_information(E, C)(s, NULL)
    expect_lt(max(abs(Hl - Hn)) / max(abs(Hn)), 1e-6)
  }
})

test_that("ecm(quadrature = \"ba81\") is the Bock-Aitkin (1981) EM", {
  d <- sim_2pl_model(N = 300, K = 5, seed = 6)
  st <- function(s, x) { e <- estep_irt(x$Y, s$a, s$b); s[c("a", "b")] <- mstep_items(e, s$a, s$b); s$.loglik <- e$loglik; s }
  ml <- run_em(list(Y = d$Y), function(x) list(a = rep(1, 5), b = rep(0, 5)), st, loglik = function(s, x) estep_irt(x$Y, s$a, s$b)$loglik, positive = "a")
  fit <- ecm(d$m, prior = "ml", quadrature = "ba81", cycles = 1)                # 41 nodes of N(0, 1), as estep_irt()
  expect_lt(max(abs(fit$estimates$est - ml$estimates$est)), 1e-3)
  expect_lt(abs(fit$loglik - ml$loglik), 1e-4)
  expect_true(attr(check_ecm(d$m, list(step_gibbs(a, b)), prior = "ml", quadrature = "ba81", n_start = 2, seed = 1), "passed"))
})

test_that("vi() uses the same steps, raises the ELBO and agrees with ECM on the means", {
  d <- sim_2pl_model(N = 300, K = 5, seed = 8)
  v <- vi(d$m, list(step_gibbs(theta, a, b), step_shift(theta = 1, b = 1)))
  expect_s3_class(v, "block_model_vi"); expect_true(v$converged)
  expect_equal(unname(v$elbo_drop[["worst"]]), "0")
  expect_true(all(diff(v$elbo[seq_len(min(10, length(v$elbo)))]) >= -1e-8 * abs(v$elbo[1])))
  e <- ecm(d$m, list(step_gibbs(a, b)))
  expect_lt(max(abs(v$summary$mean - e$estimates$est[match(v$summary$parameter, e$estimates$parameter)])), 0.1)
  al <- paste(as.character(algorithm(v, format = "latex")), collapse = "\n")                    # p and q written out, mean-field style
  expect_match(al, "q &= \\\\prod_\\{j,i\\} q\\(\\\\omega_\\{ji\\} \\\\mid \\\\theta_\\{j\\}\\) \\\\prod_\\{j\\} q\\(\\\\theta_\\{j\\}\\)")
  expect_match(al, "q\\(a_\\{i\\}\\) &= \\\\mathcal\\{N\\}\\^\\{\\+\\}")
  expect_output(print(v), "SDs of item parameters are too small")
  expect_error(vi(d$m, list(step_mala(a, b))), "no CAVI update")
  m2 <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), la[item] ~ normal(0, 1), b[item] ~ normal(0, 3), Y ~ bernoulli_logit(exp(la) * (theta - b)))
  expect_error(vi(m2, list(step_gibbs(b))), "not multilinear|no step")
  r <- check_vi(d$m, list(step_gibbs(a, b)), n_start = 2, sweeps = 10, seed = 1)
  expect_true(attr(r, "passed")); expect_output(print(r), "CAVI check")
})

test_that("a cut leaves a data part out of a Gibbs step", {
  set.seed(2); N <- 60; K <- 4
  Y <- matrix(rbinom(N * K, 1, 0.5), N); L <- matrix(rnorm(N * K, 3), N)
  m <- block_model(Y = Y, logT = L, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
                   lambda[item] ~ normal(3, 1), rho[item] ~ normal(0, 1), s2[item] ~ inv_gamma(3, 1),
                   Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - rho * theta, sqrt(s2)))
  st <- list(step_gibbs(a, b), step_gibbs(theta, cut = "logT"), step_gibbs(rho, lambda, s2))
  expect_match(al_tex(m, steps = st), "cut: logT left out", fixed = TRUE)
  expect_error(gibbs(m, list(step_gibbs(theta, cut = "T"), step_gibbs(a, b, rho, lambda, s2))), "not observed data")
  expect_error(ecm(m, st), "sampling scheme")
  f <- gibbs(m, st, n_iter = 200, n_chain = 1)
  expect_s3_class(f, "block_fit")
})

sim_rt_model <- function(seed = 5, N = 150, K = 5, var_prior = "inv_gamma") {
  set.seed(seed); cond <- set_cond(n_subj = N, n_item = K); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  LT <- unname(sim$log_t); LT[sample(length(LT), 20)] <- NA; Y <- unname(sim$Y); Y[sample(length(Y), 20)] <- NA
  m <- if (var_prior == "inv_gamma") block_model(Y = Y, logT = LT, theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, 1),
      a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
      sigma2[item] ~ inv_gamma(2, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
    else block_model(Y = Y, logT = LT, theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, sqrt(v)), v ~ half_t(3, 1),
      a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 10), rho[item] ~ normal(0, 1),
      sigma2[item] ~ half_t(3, 1), Y ~ bernoulli_logit(a * (theta - b)), logT ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  list(m = m, init = list(lambda = colMeans(LT, na.rm = TRUE)))
}

test_that("the Bock-Aitkin collapse gives the same E-step, Q, M-step and score as evaluating every person", {
  d <- sim_rt_model()
  E <- bm_setup(d$m, NULL, NULL, c(theta = 7, zeta = 7), "map", d$init, analytic = FALSE, quadrature = "ba81")
  expect_setequal(vapply(E$collapse, function(id) d$m$factors[[id]]$y, ""), c("Y", "logT"))   # a logistic and a normal part
  E2 <- E; E2$collapse <- integer()
  s <- E$s; C <- bm_centres(E, s, NULL); e1 <- bm_estep(E, s, C); e2 <- bm_estep(E2, s, C)
  expect_equal(e1$obj, e2$obj, tolerance = 1e-10)
  expect_equal(bm_qvalue(E, s, e1), bm_qvalue(E2, s, e2), tolerance = 1e-10)
  s1 <- bm_step(E, C)(s, NULL); s2 <- bm_step(E2, C)(s, NULL)
  for (v in E$fixed) expect_equal(s1[[v]], s2[[v]], tolerance = 1e-8)
  expect_equal(unlist(bm_score(E, C)(s, NULL)), unlist(bm_score(E2, C)(s, NULL)), tolerance = 1e-8)
})

test_that("ECM and Newton reach the same maximum with variance parameters (inverse-gamma and half-t priors)", {
  for (vp in c("inv_gamma", "half_t")) {
    d <- sim_rt_model(seed = 6, var_prior = vp)
    e <- ecm(d$m, init = d$init); f <- mml(d$m, init = d$init)
    expect_lt(abs(f$loglik - e$loglik), 1e-5)
    expect_lt(max(abs(f$estimates$est - e$estimates$est)), 1e-3)
    expect_lt(max(abs(f$estimates$se / e$estimates$se - 1)), 1e-3)
  }
})

test_that("the R evaluation (expressions that are not compiled) gives the same ECM fit as the compiled one", {
  d <- sim_rt_model(seed = 7, var_prior = "half_t")
  m2 <- d$m; for (k in seq_along(m2$factors)) m2$factors[[k]]$ir <- NULL
  set.seed(1); e1 <- ecm(d$m, init = d$init); set.seed(1); e2 <- ecm(m2, init = d$init)   # same starting values
  expect_lt(abs(e1$loglik - e2$loglik), 1e-8)
  expect_lt(max(abs(e1$estimates$est - e2$estimates$est)), 1e-6)
  expect_lt(max(abs(e1$estimates$se / e2$estimates$se - 1)), 1e-6)
})

test_that("CAVI with variances and an exact latent passes check_vi() and agrees with ECM", {
  d <- sim_rt_model(seed = 8, N = 200, var_prior = "half_t")
  v <- vi(d$m, init = d$init); e <- ecm(d$m, init = d$init)
  i <- grep("^(lambda|rho)", v$summary$parameter)
  expect_lt(max(abs(v$summary$mean[i] - e$estimates$est[match(v$summary$parameter[i], e$estimates$parameter)])), 0.05)
  r <- check_vi(d$m, n_start = 2, sweeps = 8, seed = 1)
  expect_true(attr(r, "passed"))
})

test_that("derivatives of compiled expressions with powers match stats::D()", {
  for (ex in list(quote(a * x^2), quote(a * (x - b)^2 + c * x), quote(x^3 * a - 2 * x))) {
    pk <- ml_deriv(ml_compile(ex), "x"); ref <- stats::D(ex, "x")
    vals <- list(a = 1.3, b = -0.4, c = 0.7, x = 0.9)
    got <- sum(vapply(pk, function(mn) mn$coef * prod(unlist(vals[mn$syms])), 0))
    expect_equal(got, eval(ref, vals), tolerance = 1e-12)
  }
})

test_that("fixes from the independent review: truncation, person-level data, ml hyperparameters, ba81 accuracy, bounds", {
  d <- sim_2pl_model(N = 300, K = 6, seed = 11)
  expect_error(block_model(Y = d$Y, theta[person] ~ normal(0, 1), mu ~ normal(0, 1), a[item] ~ normal(mu, 1, lower = 0),
                           b[item] ~ normal(0, 3), Y ~ bernoulli_logit(a * (theta - b))), "truncated prior needs a mean and sd without parameters")
  set.seed(2); th <- rowMeans(d$Y) * 3 - 1.5; g <- rbinom(300, 1, plogis(0.5 + th)); x <- th + rnorm(300, 0, 0.7)
  m <- block_model(Y = d$Y, g = g, x = x, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3),
                   c0 ~ normal(0, 3), c1 ~ normal(0, 3), d1 ~ normal(0, 3), s2 ~ inv_gamma(2, 1),
                   Y ~ bernoulli_logit(a * (theta - b)), g ~ bernoulli_logit(c0 + c1 * theta), x ~ normal(d1 * theta, sqrt(s2)))
  set.seed(1); e <- ecm(m); set.seed(1); f <- mml(m)                # person-level data (vectors)
  expect_lt(abs(f$loglik - e$loglik), 1e-5); expect_lt(max(abs(f$estimates$est - e$estimates$est)), 1e-3)
  expect_true(all(is.finite(e$estimates$se)))
  mh <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(mu_b, 1), mu_b ~ normal(0, 1),
                    Y ~ bernoulli_logit(a * (theta - b)))
  expect_error(ecm(mh, prior = "ml"), "mu_b appears only in priors")
  mb <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), a[item] ~ normal(1, 1, lower = 1.5), b[item] ~ normal(0, 3), Y ~ bernoulli_logit(a * (theta - b)))
  for (me in c("ecm", "newton")) {
    expect_warning(fb <- if (me == "ecm") ecm(mb, prior = "ml") else mml(mb, prior = "ml"), "at the bound of its prior")
    expect_true(any(is.na(fb$estimates$se[grep("^a", fb$estimates$parameter)])))
  }
})

test_that("ecm(quadrature = \"ba81\") warns when the common nodes are too coarse", {
  skip_on_cran()
  set.seed(12); cond <- set_cond(n_subj = 60, n_item = 5); sim <- sim_data(cond, sim_para(cond, "cross"), "cross")
  m <- block_model(Y = unname(sim$Y), L = unname(sim$log_t), theta[person] ~ normal(0, 1), zeta[person] ~ normal(0, sqrt(tau2)),
    a[item] ~ normal(1, 1, lower = 0), b[item] ~ normal(0, 3), lambda[item] ~ normal(0, 3), rho[item] ~ normal(0, 1),
    sigma2[item] ~ half_t(3, 1), tau2 ~ inv_gamma(2, 1), Y ~ bernoulli_logit(a * (theta - b)), L ~ normal(lambda - zeta - rho * theta, sqrt(sigma2)))
  init <- list(lambda = colMeans(sim$log_t))
  expect_warning(fb <- ecm(m, init = init, quadrature = "ba81", analytic = FALSE, se = FALSE), "too coarse.*analytic = TRUE")
  expect_gt(abs(fb$adaptive_gap), 0.5)
  f <- suppressWarnings(ecm(m, init = init, quadrature = "ba81", se = FALSE))                # speed exact, theta on 41 common nodes
  expect_lt(abs(f$adaptive_gap), 0.5)
})

test_that("ecm() is ECM only; mml() maximizes directly", {
  d <- sim_2pl_model(N = 100, K = 4, seed = 3)
  m <- block_model(Y = d$Y, theta[person] ~ normal(0, 1), la[item] ~ normal(0, 1), b[item] ~ normal(0, 3), Y ~ bernoulli_logit(exp(la) * (theta - b)))
  expect_warning(ecm(m, method = "newton"), "unused arguments: method")
  expect_error(mml(rtirt_null(sim_data(set_cond(n_subj = 50, n_item = 3), sim_para(set_cond(n_subj = 50, n_item = 3), "null"), "null"))), "block_model")
})

test_that("starting values from the data: intercepts match the data; ECM needs no init", {
  d <- sim_rt_model(seed = 6, var_prior = "half_t")
  s <- birtRcpp:::init_model_state(d$m, list())
  f <- Filter(function(z) z$observed && z$y == "logT", d$m$factors)[[1]]
  mk <- f$mask; logT <- d$m$data0$logT                                                       # missing cells: 0 in the data, 0 in the mask
  expect_equal(unname(s$lambda), unname(colSums(logT * mk) / colSums(mk)), tolerance = 1e-8)  # lambda - zeta - rho * theta at zeta = theta = 0
  s0 <- birtRcpp:::init_model_state(d$m, list(), from_data = FALSE)
  expect_lt(max(abs(s0$lambda)), 1)                                                         # the prior means (0) without the data
  e0 <- ecm(d$m); e1 <- ecm(d$m, init = d$init)
  expect_lt(abs(e0$loglik - e1$loglik), 1e-6)
  expect_lt(max(abs(e0$estimates$est - e1$estimates$est)), 1e-3)
})
