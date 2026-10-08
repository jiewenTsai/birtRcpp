test_that("equations() for every model, Unicode and LaTeX", {
  cd <- set_cond(n_subj = 60, n_item = 4, n_feat = 2)
  for (ty in c("mlirt", "latreg", "latent", "cross")) {
    d <- sim_data(cd, sim_para(cd, ty), ty)
    m <- switch(ty, mlirt = mlirt(d), latreg = rtirt_latreg(d), latent = rtirt_latent(d), cross = rtirt_cross(d))
    u <- as.character(equations(m)); l <- as.character(equations(m, "latex"))
    expect_true(any(grepl("logit P(y", u, fixed = TRUE)))
    expect_true(any(grepl("Priors", u)))
    expect_equal(l[1], "\\begin{aligned}")
    expect_false(any(grepl("[^\\x01-\\x7f]", l, perl = TRUE)))
  }
  d <- sim_data(cd, sim_para(cd, "cross"), "cross")
  expect_true(any(grepl("ALD", as.character(equations(rtirt_cross(d, quantile = 0.25))))))
  expect_true(any(grepl("Quantiles", as.character(equations(rtirt_cross(d, quantile = c(0.25, 0.5)))))))
  e <- ecm(rtirt_cross(d), se = FALSE)
  expect_false(any(grepl("Priors", as.character(equations(e)))))
})
