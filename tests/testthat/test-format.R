# Printed numbers: 3 decimals in every print() and summary(), full precision in the objects

decimals_of <- function(text) {
  numbers <- unlist(regmatches(text, gregexpr("-?[0-9]+\\.[0-9]+", text)))
  unique(nchar(sub("^[^.]*\\.", "", numbers)))
}

format_fit <- local({
  fits <- NULL
  function() {
    if (is.null(fits)) {
      set.seed(5)
      cond <- set_cond(n_subj = 300, n_item = 5, n_feat = 1)
      dat <- sim_data(cond, sim_para(cond, "latreg"), "latreg")
      fits <<- list(null = ecm(dat, model = "null"), latreg = ecm(dat, model = "latreg"), dat = dat)
    }
    fits
  }
})

test_that("the formatting helpers follow the rules", {
  expect_identical(fmt_fixed(c(1.23456, -0.5, NA, 1234567.1)), c("1.235", "-0.500", "", "1234567.100"))
  expect_identical(fmt_fixed(1e-12), "0.000")                      # no scientific notation
  expect_identical(fmt_p(c(0.5, 0.0004, 2.2e-16, NA)), c("0.500", "<0.001", "<0.001", ""))
  expect_identical(p_text(1e-20), "p < 0.001")
  expect_identical(p_text(0.0456), "p = 0.046")
  expect_identical(fmt_column(c(300, 200), "N"), c("300", "200"))   # counts without decimals
  expect_identical(fmt_column(c(1.5, 2), "N"), c("1.500", "2.000")) # not whole: decimals
  expect_identical(fmt_column(c(3L, NA), "anything"), c("3", ""))
  expect_identical(fmt_column(c(0.01, 1e-5), "p.value"), c("0.010", "<0.001"))
  expect_identical(fmt_column(c("x", NA)), c("x", ""))
})

test_that("estimates() prints 3 decimals and keeps full precision", {
  fit <- format_fit()$latreg
  est <- estimates(fit)
  expect_s3_class(est, c("rtirt_tbl", "tbl_df"))
  text <- capture.output(print(est))
  expect_identical(decimals_of(text), 3L)
  expect_match(text[1], "^# A tibble: [0-9]+ x 7$")
  expect_false(any(grepl("e-[0-9]", text)))
  expect_length(text, nrow(est) + 2)                                # all rows of a short table
  # stored values are not rounded
  expect_false(isTRUE(all.equal(est$est, round(est$est, 3))))
  expect_identical(est$est, fit$ecm$estimates$est)
  # tibble, dplyr and data frame use
  plain <- as.data.frame(est)
  expect_identical(class(plain), "data.frame")
  expect_identical(plain$est, est$est)
  filtered <- dplyr::filter(est, grepl("^b\\[", parameter))
  expect_s3_class(filtered, "rtirt_tbl")
  expect_identical(filtered$est, est$est[grepl("^b\\[", est$parameter)])
  expect_identical(decimals_of(capture.output(print(dplyr::mutate(est, twice = 2 * est)))), 3L)
  expect_length(capture.output(print(est, n = 3)), 3 + 3)
})

test_that("long tables show the first rows and say how many more", {
  fit <- format_fit()$null
  sc <- scores(fit)
  text <- capture.output(print(sc))
  expect_length(text, 1 + 1 + 20 + 1)
  expect_match(text[length(text)], "280 more rows")
  expect_identical(decimals_of(text), 3L)
  expect_length(capture.output(print(sc, n = Inf)), 2 + nrow(sc))
  expect_identical(sc$ability, fit$ecm$person$theta$mean)
})

test_that("print(), summary(), fit_indices() and anova() show 3 decimals", {
  fits <- format_fit()
  printed <- capture.output(print(fits$latreg))
  expect_match(printed, "AIC [0-9]+\\.[0-9]{3}, BIC [0-9]+\\.[0-9]{3}$", all = FALSE)
  s <- capture.output(print(summary(fits$latreg)))
  expect_identical(decimals_of(s[!grepl("ECM, ", s)]), 3L)          # the run time keeps 1 decimal
  expect_match(s, "<0.001", fixed = TRUE, all = FALSE)               # P(>|z|) of the structural parameters
  fi <- capture.output(print(fit_indices(fits$latreg)))
  expect_match(fi[3], "^ECM +-[0-9]+\\.[0-9]{3} +[0-9]+ +300 ")       # npar and N without decimals
  a <- anova(fits$null, fits$latreg)
  text <- capture.output(print(a))
  expect_match(text, "<0.001  \\*\\*\\*$", all = FALSE)               # p-value and stars
  expect_false(any(grepl("e-[0-9]", text)))
  expect_identical(decimals_of(text[!startsWith(text, "Signif. codes")]), 3L)
  expect_identical(a$logLik, c(fits$null$ecm$logLik, fits$latreg$ecm$logLik))
})

test_that("coef(), reliability() and the tests print 3 decimals", {
  fits <- format_fit()
  cf <- coef(fits$latreg)
  expect_identical(decimals_of(capture.output(print(cf))), 3L)
  expect_identical(unclass(cf)[["b[item01]"]], fits$latreg$ecm$estimates$est[6])
  expect_identical(decimals_of(capture.output(print(cf[1:2]))), 3L)
  rel <- reliability(fits$latreg)
  expect_identical(decimals_of(capture.output(print(rel))), 3L)
  expect_false(isTRUE(all.equal(as.numeric(rel), round(as.numeric(rel), 3))))
  ci <- capture.output(print(ci_test(fits$null)))
  expect_identical(decimals_of(ci), 3L)
  pk <- preknowledge_test(fits$null, compromised = 1:2)
  expect_identical(decimals_of(capture.output(print(pk))), 3L)
  expect_s3_class(pk, "data.frame")
})

test_that("knitr documents print the tables as at the console (not with kable())", {
  skip_if_not_installed("knitr")
  fits <- format_fit()
  est <- estimates(fits$latreg)
  expect_identical(capture.output(knitr::knit_print(est)), capture.output(print(est)))
  a <- anova(fits$null, fits$latreg)
  expect_identical(capture.output(knitr::knit_print(a)), capture.output(print(a)))
})
