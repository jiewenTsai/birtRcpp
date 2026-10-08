test_that("check_sampler() passes a correct sampler and catches a wrong prior", {
  skip_on_cran()
  set.seed(1); cond <- set_cond(n_subj = 50, n_item = 3)
  m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
  narrow <- rtirt_priors(lambda = c(mean = 0, sd = 1))          # fast mixing of the joint chain
  g <- check_sampler(m, n_iter = 1e5, priors = narrow, seed = 2)
  expect_s3_class(g, "rtirt_check"); expect_true(attr(g, "passed"))
  expect_output(print(g), "Geweke")
  bad <- check_sampler(m, n_iter = 2e4, priors = narrow, sampler_priors = rtirt_priors(lambda = c(mean = 0, sd = 1), a = c(mean = 1, sd = 0.3)), seed = 2)
  expect_false(attr(bad, "passed"))
  expect_true(all(grepl("^a\\[", head(bad$parameter[order(bad$p.value)], 3))))
  s <- check_sampler(m, method = "sbc", n_rep = 60, n_iter = 4000, priors = narrow, sampler_priors = rtirt_priors(lambda = c(mean = 0, sd = 1), b = c(mean = 1, sd = 1)), seed = 3)
  expect_false(attr(s, "passed")); expect_equal(dim(attr(s, "ranks")), c(60L, 16L))
  expect_error(check_sampler(m, priors = rtirt_priors_julia()), "do not depend on the data")
})

test_that("rank-normalized diagnostics (Vehtari et al., 2021) are available as an option", {
  skip_if_not_installed("posterior")
  set.seed(4); cond <- set_cond(n_subj = 100, n_item = 4)
  m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
  g <- gibbs(m, n_iter = 4000, n_chain = 2, seed = 1, verbose = FALSE, fit_indices = FALSE)
  e0 <- estimates(g); e1 <- estimates(g, diagnostics = "rank")
  expect_equal(e0$est, e1$est); expect_true(all(c("ess", "ess_tail") %in% names(e1)))
  expect_true(all(e1$rhat > 0.99 & e1$rhat < 1.1, na.rm = TRUE))
  # with very short chains the rank-normalized R-hat is the stricter one (it sees the skewed variances)
  s <- gibbs(m, n_iter = 600, n_chain = 2, seed = 1, verbose = FALSE, fit_indices = FALSE)
  expect_gt(max(estimates(s, diagnostics = "rank")$rhat), max(estimates(s)$rhat))
  expect_output(convergence(g, diagnostics = "rank"), "bulk and tail")
})

test_that("person_fit() and ci_test() work and behave under the model", {
  set.seed(5); cond <- set_cond(n_subj = 300, n_item = 6)
  d <- sim_data(cond, sim_para(cond, "null"), "null")                  # data of the fitted model
  f <- gibbs(rtirt_null(d), n_iter = 800, n_chain = 1, seed = 1, verbose = FALSE, fit_indices = FALSE)
  pf <- person_fit(f)
  expect_equal(nrow(pf), 300); expect_true(all(pf$p_flag >= 0 & pf$p_flag <= 1)); expect_lt(attr(pf, "share"), 0.05)
  expect_error(person_fit(gibbs(rtirt_cross(d, quantile = 0.5), n_iter = 100, n_chain = 1, verbose = FALSE, fit_indices = FALSE)), "normal log times")
  ci <- ci_test(ecm(rtirt_null(d), se = FALSE))
  expect_s3_class(ci, "rtirt_ci_test"); expect_equal(nrow(ci), 6); expect_gt(attr(ci, "overall")[["p.value"]], 0.001)
  expect_output(print(ci), "van der Linden")
  # a shift of the log times of correct responses is detected
  d2 <- d; d2$log_t <- d$log_t + 0.3 * d$Y
  expect_lt(attr(ci_test(ecm(rtirt_null(d2), se = FALSE)), "overall")[["p.value"]], 1e-6)
})

test_that("preknowledge_test() has the right size and detects faster responses on compromised items", {
  set.seed(6); cond <- set_cond(n_subj = 600, n_item = 20)
  d <- sim_data(cond, sim_para(cond, "null"), "null")
  pk <- preknowledge_test(ecm(rtirt_null(d), se = FALSE), compromised = 1:6)
  expect_equal(nrow(pk), 600); expect_lt(mean(pk$p_Lambda_s < 0.05), 0.09); expect_lt(mean(pk$p_chi_pf < 0.05), 0.12)
  d2 <- d; d2$log_t[1:60, 1:6] <- d2$log_t[1:60, 1:6] - 1                 # 60 test takers much faster on the compromised items
  pk2 <- preknowledge_test(ecm(rtirt_null(d2), se = FALSE), compromised = 1:6)
  expect_gte(mean(pk2$p_Lambda_s[1:60] < 0.01), 0.7)
  expect_error(preknowledge_test(ecm(rtirt_cross(d), se = FALSE), 1:6), "rtirt_null")
})

test_that("check_ecm() checks the EM of ecm() for every built-in normal model and catches a wrong M-step", {
  skip_on_cran()
  set.seed(1); cond <- set_cond(n_subj = 200, n_item = 5, n_feat = 2)
  mods <- list(mlirt = function(d) mlirt(d), null = function(d) rtirt_null(d, speed_var = "fixed"), latreg = function(d) rtirt_latreg(d),
               latent = function(d) rtirt_latent(d), cross = function(d) rtirt_cross(d), cross_1pl = function(d) rtirt_cross(d, itemtype = "1pl"))
  sims <- c(mlirt = "mlirt", null = "null", latreg = "latreg", latent = "latent", cross = "cross", cross_1pl = "cross")
  for (n in names(mods)) {
    m <- mods[[n]](sim_data(cond, sim_para(cond, sims[[n]]), sims[[n]]))
    for (pr in c("ml", "map")) {
      r <- check_ecm(m, prior = pr, n_start = 2, seed = 1)
      expect_s3_class(r, "em_check"); expect_true(attr(r, "passed"), label = paste(n, pr)); expect_equal(nrow(r), 4L)
    }
  }
  expect_output(print(r), "all checks passed")
  ms <- birtRcpp:::em_mstep; sc <- birtRcpp:::em_score
  wrong <- list(function(p, E, es) { n <- ms(p, E, es); n$sigma2t <- 1.1 * n$sigma2t; n },     # a variance step off its maximum
                function(p, E, es) { n <- ms(p, E, es); n$rho <- n$rho + 0.05; n })
  for (w in wrong) local({
    local_mocked_bindings(em_mstep = w)
    r <- check_ecm(m, prior = "ml", n_start = 2, seed = 1)
    expect_false(attr(r, "passed")); expect_false(r$passed[r$check == "ECM stops at a maximum: largest |gradient|"])
  })
  local({                                                                    # a wrong score breaks Fisher's identity
    local_mocked_bindings(em_score = function(x, E) { s <- sc(x, E); s[1] <- 1.05 * s[1]; s })
    r <- check_ecm(m, prior = "ml", n_start = 2, seed = 1)
    expect_false(r$passed[grepl("Fisher", r$check)])
  })
  expect_error(check_ecm(rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"), quantile = 0.3)), "quantile models")
})
