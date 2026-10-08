quiet <- function(expr) { utils::capture.output(x <- expr); x }

# small simulations: every model runs, methods work, estimates track the truth
fit_sim <- function(model, cond, quantile = NULL, ...) {
  d <- sim_data(cond, sim_para(cond, model), model)
  m <- switch(model, cross = rtirt_cross(d, quantile, cond = cond), latent = rtirt_latent(d, quantile, cond = cond),
              latreg = rtirt_latreg(d, cond = cond), null = rtirt_null(d, cond = cond), mlirt = mlirt(d, cond = cond))
  gibbs(m, verbose = FALSE, ...)
}

test_that("cross-relation models recover rho and item parameters", {
  set.seed(10)
  cond <- set_cond(n_subj = 800, n_item = 8, n_iter = 1500, n_chain = 2)
  m <- fit_sim("cross", cond)
  expect_s3_class(m, "rtirt_cross"); expect_s3_class(m, "rtirt")
  expect_gt(compare_para(m, "rho")$corr, 0.9)
  expect_gt(compare_para(m, "a")$corr, 0.7)
  expect_gt(compare_para(m, "lambda")$corr, 0.99)
  expect_lt(abs(compare_para(m, "rho")$bias), 0.05)
  expect_equal(nrow(dic(m)), 1)
  expect_s3_class(coda::as.mcmc.list(m), "mcmc.list")
  qr <- fit_sim("cross", cond, quantile = 0.5)
  expect_true(qr$qr); expect_equal(qr$q, 0.5)
  expect_gt(compare_para(qr, "rho")$corr, 0.85)               # median model, symmetric errors
})

test_that("latent regression of speed recovers beta and the ability effect", {
  set.seed(11)
  cond <- set_cond(n_subj = 1000, n_item = 10, n_feat = 2, n_iter = 1500, n_chain = 2)
  d <- sim_data(cond, sim_para(cond, "latent"), "latent"); tp <- attr(d, "true_para")
  m <- gibbs(rtirt_latent(d), n_iter = 1500, n_chain = 2, verbose = FALSE)
  expect_identical(m$true_para, tp)                            # picked up from sim_data()
  expect_lt(max(abs(m$post$mean$beta - tp$beta)), 0.06)
  expect_lt(abs(m$post$mean$b_theta - tp$b_theta), 0.06)
  expect_equal(compare_para(m, "rho")[, -1], compare_para(m, "b_ability")[, -1])   # Julia name
  expect_equal(unname(m$post$mean$sigma_p[2, 2]), tp$sigma_p[2, 2], tolerance = 0.25)
  q <- gibbs(rtirt_latent(d, quantile = 0.5), n_iter = 1500, n_chain = 2, verbose = FALSE)
  expect_lt(max(abs(q$post$mean$beta - tp$beta)), 0.08)     # median = mean for normal errors
})

test_that("joint latent regression model recovers the regressions and the correlation", {
  set.seed(12)
  cond <- set_cond(n_subj = 1000, n_item = 10, n_feat = 2, n_iter = 1500, n_chain = 2)
  d <- sim_data(cond, sim_para(cond, "latreg", true_corr = 0.4), "latreg"); tp <- attr(d, "true_para")
  m <- gibbs(rtirt_latreg(d, cond = cond), verbose = FALSE)        # speed_var = "free" by default
  # against the regressions of the realized theta / zeta on X (what this sample can tell), within 3 SD
  ols <- cbind(coef(lm(tp$theta ~ d$X - 1)), coef(lm(tp$zeta ~ d$X - 1)))
  sdb <- matrix(estimates(m, "^beta")$sd, ncol = 2)
  expect_true(all(abs(m$post$mean$beta - ols) < 3 * sdb))
  S <- m$post$mean$sigma_p
  expect_equal(dimnames(S), list(c("ability", "speed"), c("ability", "speed")))
  expect_equal(S[1, 2] / sqrt(S[2, 2]), 0.4, tolerance = 0.25)
  expect_equal(S[2, 2], tp$sigma_p[2, 2], tolerance = 0.2)
  n <- gibbs(suppressWarnings(rtirt_null(d, speed_var = "fixed", cond = cond)), verbose = FALSE)
  expect_equal(unname(n$post$mean$sigma_p[2, 2]), 1)
  p <- estimates(m)
  expect_true(all(c("est", "sd", "q025", "q975", "rhat", "ess", "sig") %in% names(p)))
  expect_true(all(c("cor_ability_speed", "var_speed") %in% p$parameter))
})

test_that("accessors: coef, summary, scores, reliability, fit_indices, plot", {
  set.seed(14)
  cond <- set_cond(n_subj = 300, n_item = 6)
  d <- sim_data(cond, sim_para(cond, "cross"), "cross")
  m <- gibbs(rtirt_cross(d), n_iter = 600, n_chain = 2, seed = 3, verbose = FALSE)
  m2 <- gibbs(rtirt_cross(d), n_iter = 600, n_chain = 2, seed = 3, verbose = FALSE)
  expect_identical(coef(m), coef(m2))                          # seed
  cf <- coef(m)
  expect_true(is.numeric(cf) && !is.null(names(cf)))
  expect_true(all(c("a[item01]", "rho[item06]", "var_speed") %in% names(cf)))
  s <- scores(m)
  expect_equal(names(s), c("id", "ability", "ability_psd", "speed", "speed_psd")); expect_equal(nrow(s), 300)
  r <- reliability(m); expect_equal(names(r), c("ability", "speed")); expect_true(all(r > 0 & r < 1))
  expect_equal(fit_indices(m)$DIC_complete, dic(m)$dic)
  s <- summary(m); expect_output(print(s), "[$]items"); expect_equal(nrow(s$items), ncol(m$data$Y))
  expect_true(all(nzchar(as.data.frame(s)$prior[!as.data.frame(s)$fixed])))
  expect_output(print(m), "Reliability")
  expect_output(cv <- convergence(m), "parameters with ESS")
  expect_true(cv$share >= 0 && cv$share <= 1)                       # was NaN without tail ESS
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(m)); expect_silent(plot(m, "trace"))
  qs <- gibbs(rtirt_cross(d, quantile = c(0.25, 0.75)), n_iter = 400, n_chain = 2, verbose = FALSE)
  expect_s3_class(qs, "rtirt_qset")
  expect_equal(dim(coef(qs)), c(length(cf), 2))
  expect_silent(plot(qs))
})

test_that("input_data: column names, ids, and checks", {
  set.seed(15)
  raw <- data.frame(pid = 101:150, i1 = rbinom(50, 1, .5), i2 = rbinom(50, 1, .5), t1 = rlnorm(50, 3), t2 = rlnorm(50, 3), x = rnorm(50))
  d <- input_data(resp = c("i1", "i2"), time = c("t1", "t2"), cov = "x", id = "pid", data = raw)
  expect_equal(d$items, c("i1", "i2")); expect_equal(d$id, 101:150); expect_equal(colnames(d$X), "x")
  expect_equal(d$log_t, log(as.matrix(raw[c("t1", "t2")])), ignore_attr = TRUE)
  expect_warning(input_data(raw[c("i1", "i2")], log(raw[c("t1", "t2")])), "already log times")
  expect_warning(input_data(raw[c("i1", "i2")], log_time = raw[c("t1", "t2")]), "raw seconds")
  expect_warning(input_data(raw[c("i1", "i2")], cov = raw$x + 5), "not centred")
  expect_error(input_data(raw[c("i1", "i2")] + 1), "recode it to 0/1")
  r2 <- raw; r2$i1[3] <- NA
  expect_error(input_data(c("i1", "i2"), data = r2), "missing values in 1 persons")
  expect_error(input_data(c("i1", "zz"), data = raw), "zz")
  expect_warning(rtirt_cross(d), "ignored")
  expect_error(rtirt_latent(input_data(raw[c("i1", "i2")], raw[c("t1", "t2")])), "needs covariates")
  expect_error(rtirt_cross(set_cond(), d), "the data come first")
  expect_error(suppressWarnings(rtirt_cross(d, quantile = 1.2)), "quantile")
})

test_that("multilevel IRT and dimension checks", {
  set.seed(13)
  cond <- set_cond(n_subj = 800, n_item = 10, n_feat = 2, n_iter = 1000, n_chain = 2)
  m <- fit_sim("mlirt", cond, itemtype = "2pl")
  expect_lt(max(abs(m$post$mean$beta - m$true_para$beta)), 0.15)
  expect_error(rtirt_cross(input_data(matrix(c(0, 1, 1, 0), 4, 5), matrix(exp(rnorm(20, 3)), 4)), cond = cond), "persons")
  expect_error(input_data(matrix(2, 3, 3)), "0/1")
  expect_error(estimates(mlirt(sim_data(cond, sim_para(cond, "mlirt"), "mlirt"))), "gibbs")
})

test_that("Julia aliases give the same objects", {
  Cond <- setCond(nSubj = 200, nItem = 5, nFeat = 0, nIter = 300, nChain = 2, qRt = 0.25)
  expect_identical(Cond, set_cond(n_subj = 200, n_item = 5, n_feat = 0, n_iter = 300, n_chain = 2, q_rt = 0.25))
  set.seed(1); tp <- setTrueParaRtIrtCross(Cond); Data <- setDataRtIrtCross(Cond, tp)
  M <- GibbsRtIrtCrossQr(Cond, Data, attr(Data, "true_para"))
  expect_s3_class(M, "rtirt_cross"); expect_true(M$qr); expect_equal(M$q, 0.25)
  expect_equal(M$settings$speed_var, "fixed")                    # Julia default cov2one = true
  M <- sampleGibbs(M, verbose = FALSE)
  expect_equal(unname(M$post$mean$sigma_p[2, 2]), 1)
  expect_equal(getDic(M), dic(M))
  expect_equal(comparePara(M, "rho"), compare_para(M, "rho"))
  expect_s3_class(GibbsRtIrt(setCond(nSubj = 50, nItem = 3, nFeat = 1),
                             InputData(matrix(rbinom(150, 1, .5), 50), matrix(exp(rnorm(150, 3)), 50), matrix(rnorm(50)))),
                  "rtirt_latreg")
})

test_that("items without variation are refused", {
  Y <- matrix(c(0, 1, 1, 0, 1), 5, 4); Y[, 2] <- 1
  expect_error(input_data(Y, matrix(exp(rnorm(20, 3)), 5)), "without variation")
})

test_that("gibbs() takes one set of argument names and flags unknown ones; cov2one only in sampleGibbs()", {
  cond <- set_cond(n_subj = 100, n_item = 4)
  m <- rtirt_null(sim_data(cond, sim_para(cond, "null"), "null"))
  f <- gibbs(m, n_iter = 200, n_chain = 1, n_burnin = 50, n_thin = 2, seed = 1, verbose = FALSE)
  expect_equal(f$cond$n_chain, 1); expect_equal(f$cond$n_burnin, 50); expect_equal(f$cond$n_thin, 2)
  expect_warning(gibbs(m, n_iter = 100, n_chain = 1, n_iters = 5, verbose = FALSE), "unused arguments: n_iters")
  expect_warning(gibbs(m, n_iter = 100, n_chain = 1, n_chains = 2, verbose = FALSE), "unused arguments: n_chains")
  g <- sampleGibbs(m, cov2one = TRUE, n_iter = 100, n_chain = 1, verbose = FALSE)
  expect_equal(g$settings$speed_var, "fixed")
})

test_that("collapse = TRUE samples the same posterior as the original sampler", {
  set.seed(7)
  cond <- set_cond(n_subj = 400, n_item = 8, n_feat = 2)
  data <- list(mlirt = sim_data(cond, sim_para(cond, "mlirt"), "mlirt"),
               latreg = sim_data(cond, sim_para(cond, "latreg"), "latreg"),
               latent = sim_data(cond, sim_para(cond, "latent"), "latent"),
               cross = sim_data(cond, sim_para(cond, "cross"), "cross"))
  models <- list(mlirt(data$mlirt), rtirt_latreg(data$latreg), rtirt_latreg(data$latreg, speed_var = "fixed"),
                 rtirt_latent(data$latent), rtirt_latent(data$latent, quantile = 0.25),
                 rtirt_cross(data$cross), rtirt_cross(data$cross, quantile = 0.25, speed_var = "fixed"))
  for (m in models) {
    s <- gibbs(m, n_iter = 1500, n_chain = 2, seed = 1, verbose = FALSE, fit_indices = FALSE)
    c <- gibbs(m, n_iter = 1500, n_chain = 2, seed = 2, verbose = FALSE, fit_indices = FALSE, collapse = TRUE)
    expect_true(c$post$settings$collapse)
    expect_true(c$post$ab_accept > 0.5 && c$post$ab_accept <= 1)
    es <- estimates(s); ec <- estimates(c)
    ok <- es$sd > 0
    z <- (ec$est - es$est)[ok] / sqrt(es$sd^2 / pmax(es$ess, 50) + ec$sd^2 / pmax(ec$ess, 50))[ok]
    expect_lt(max(abs(z)), 5)
    expect_gt(stats::cor(s$post$person$theta$mean, c$post$person$theta$mean), 0.995)
  }
  expect_error(gibbs(models[[1]], n_iter = 100, collapse = NA, verbose = FALSE), "collapse")
})

test_that("every camelCase (Julia) alias calls only functions that exist", {
  ns <- asNamespace("birtRcpp")
  aliases <- grep("^([A-Z]|[a-z]+[A-Z])", getNamespaceExports("birtRcpp"), value = TRUE)
  aliases <- aliases[!grepl("\\.", aliases)]
  heads <- function(e) {                                   # names of the functions a body calls
    if (!is.call(e)) return(character())
    h <- if (is.symbol(e[[1]])) as.character(e[[1]]) else character()
    c(h, unlist(lapply(as.list(e)[-1], heads)))
  }
  expect_gt(length(aliases), 20)
  for (a in aliases) {
    f <- get(a, envir = ns); expect_true(is.function(f), info = a)
    calls <- setdiff(unique(heads(body(f))), c("function", "<-", "{", "(", "if", "for", "return", "::", ":::", "$", "[[", "["))
    calls <- setdiff(calls, names(formals(f)))                     # arguments that are functions (e.g. funcGibbs)
    missing <- calls[!vapply(calls, exists, logical(1), envir = ns, mode = "function")]
    expect_identical(missing, character(), info = a)
  }
})
