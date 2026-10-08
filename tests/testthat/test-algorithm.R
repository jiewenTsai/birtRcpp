test_that("algorithm() writes every sampler and the EM steps in Unicode and LaTeX", {
  set.seed(8)
  cond <- set_cond(n_subj = 120, n_item = 5)
  d <- sim_data(cond, sim_para(cond, "cross"), "cross")
  X <- matrix(rnorm(120), ncol = 1, dimnames = list(NULL, "x"))
  dc <- input_data(resp = d$Y, log_time = d$log_t, cov = X)
  models <- list(mlirt(dc), mlirt(dc, itemtype = "1pl"), rtirt_null(d), rtirt_null(d, speed_var = "fixed"), rtirt_latreg(dc),
                 rtirt_latent(dc), rtirt_latent(dc, quantile = 0.2), rtirt_cross(d), rtirt_cross(d, quantile = 0.7),
                 rtirt_cross(d, speed_var = "fixed"))
  for (m in models) for (co in c(FALSE, TRUE)) {
    u <- as.character(algorithm(m, collapse = co))
    expect_false(any(grepl("\\\\", u)), info = paste(class(m)[1], co))     # no LaTeX left in the Unicode text
    expect_true(any(grepl("draw_omega()", u, fixed = TRUE)))
    expect_equal(any(grepl("draw_items_joint", u)), co && !identical(m$settings$itemtype, "1pl"))
    l <- as.character(algorithm(m, collapse = co, format = "latex"))
    expect_equal(l[c(1, length(l))], c("\\begin{aligned}", "\\end{aligned}"))
    if (!isTRUE(m$qr)) expect_true(any(grepl("SQUAREM", as.character(algorithm(m, method = "ecm")))))
  }
  u <- as.character(algorithm(rtirt_cross(d)))
  expect_true(any(grepl("Cross-relation move", u))); expect_true(any(grepl("Huang & Wand", u)))
  expect_true(any(grepl("IG(0.001 + N/2", as.character(algorithm(gibbs(rtirt_cross(d), n_iter = 20, n_chain = 1, verbose = FALSE,
    fit_indices = FALSE, priors = rtirt_priors_julia()))), fixed = TRUE)))
  expect_true(any(grepl("collapse = TRUE", as.character(algorithm(gibbs(rtirt_null(d), n_iter = 20, n_chain = 1, verbose = FALSE,
    fit_indices = FALSE, collapse = TRUE))), fixed = TRUE)))
  expect_true(any(grepl("E-step", as.character(algorithm(ecm(rtirt_null(d), se = FALSE))))))
  expect_error(algorithm(rtirt_cross(d, quantile = 0.5), method = "ecm"), "quantile")
  expect_s3_class(algorithm(rtirt_cross(d, quantile = c(0.2, 0.8))), "rtirt_algorithm")
})

test_that("tex_to_unicode() converts the LaTeX subset", {
  f <- birtRcpp:::tex_to_unicode
  expect_equal(f("\\theta_i \\mid \\cdot"), "θᵢ | ·")
  expect_equal(f("\\sigma_j^{-2}"), "σⱼ⁻²")
  expect_equal(f("x_i^\\top\\beta_1"), "xᵢᵀβ₁")
  expect_equal(f("\\sqrt{\\psi_j/\\chi_{ij}}"), "√(ψⱼ/χᵢⱼ)")
  expect_equal(f("\\text{half-t} - c"), "half-t − c")
  expect_equal(f("Q_{\\theta\\zeta}"), "Q_θζ")
})
