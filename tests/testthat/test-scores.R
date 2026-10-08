test_that("casewise scores are the gradients of the per-person marginal log-likelihoods", {
  set.seed(5)
  cond <- set_cond(n_subj = 200, n_item = 5, n_feat = 1)
  for (m in list(rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross")),
                 rtirt_latreg(sim_data(cond, sim_para(cond, "latreg"), "latreg")),
                 rtirt_latent(sim_data(cond, sim_para(cond, "latent"), "latent")),
                 mlirt(sim_data(cond, sim_para(cond, "mlirt"), "mlirt")))) {
    f <- suppressWarnings(ecm(m, se = FALSE)); S <- estfun.rtirt(f)
    E <- birtRcpp:::em_setup(f, f$ecm$nodes, FALSE); x <- f$ecm$x
    ll <- function(x) birtRcpp:::em_estep(birtRcpp:::em_unpack(x, E), E)$ll
    G <- sapply(seq_along(x), function(j) { h <- 1e-5 * max(1, abs(x[j])); e <- replace(numeric(length(x)), j, h); (ll(x + e) - ll(x - e)) / (2 * h) })
    expect_lt(max(abs(S - G)) / max(1, abs(G)), 1e-4)        # mlirt: GH nodes move with beta (discretisation ~3e-5)
    expect_equal(unname(colSums(S)), unname(birtRcpp:::em_score(x, E)), tolerance = 1e-8)
  }
})

test_that("score_test() runs, refuses what it cannot test and handles clusters", {
  skip_if_not_installed("strucchange")
  set.seed(6)
  cond <- set_cond(n_subj = 300, n_item = 6)
  d <- sim_data(cond, sim_para(cond, "cross"), "cross"); f <- ecm(rtirt_cross(d), se = FALSE)
  r <- score_test(f, rowSums(d$log_t), by_item = TRUE)
  expect_s3_class(r, "rtirt_score_test"); expect_equal(nrow(r), 7); expect_true(all(r$p.value >= 0 & r$p.value <= 1))
  expect_output(print(r), "DM test")
  expect_equal(attr(score_test(f, factor(rep(1:2, 150))), "test"), "LM")
  expect_error(score_test(f, 1:10), "persons")
  expect_error(score_test(f, rnorm(300), cluster = rep(1:30, each = 10)), "constant within clusters")
  r2 <- score_test(f, rep(rnorm(30), each = 10), cluster = rep(1:30, each = 10))
  expect_equal(attr(r2, "clusters"), 30)
  expect_error(estfun.rtirt(ecm(rtirt_cross(d), se = FALSE, prior = "map")), "maximum likelihood")
})

