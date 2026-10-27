# Regression test of the refactor: every case of helper-regression.R against the baseline in
# fixtures/regression.rds (made by fixtures/make-regression.R on the same computer). Cases with
# `fit_new` are checked twice: through the old interface and through ecm(data, model = ...) /
# gibbs(data, model = ...) of task 3.1; both must reproduce the same baseline.

test_that("the regression baseline exists", {
  skip_on_cran()
  skip_on_ci()
  expect_true(file.exists(regression_fixture()))
})

baseline <- if (file.exists(regression_fixture())) readRDS(regression_fixture())

expect_snapshot_match <- function(new, old, kind, label) {
  if (kind == "ecm") {
    expect_table_close(new$estimates, old$estimates, 1e-6, paste(label, "estimates"))
    expect_lte(abs(new$logLik - old$logLik), 1e-8)
    expect_table_close(new$scores, old$scores, 1e-6, paste(label, "scores"))
    for (tab in names(old$tables)) {
      expect_table_close(new$tables[[tab]], old$tables[[tab]], 1e-6, paste(label, "summary", tab))
    }
    expect_identical(new$equations, old$equations)
    expect_identical(new$algorithm, old$algorithm)
    expect_identical(drop_hints(new$print), drop_hints(old$print))
  }

  if (kind == "gibbs") {
    how <- expect_draws_match(new, old, label)
    if (identical(how, "identical")) {
      expect_table_close(new$scores, old$scores, 1e-10, paste(label, "scores"))
      for (tab in names(old$tables)) {
        expect_table_close(new$tables[[tab]], old$tables[[tab]], 1e-8, paste(label, "summary", tab))
      }
    }
    expect_identical(new$equations, old$equations)
    expect_identical(new$algorithm, old$algorithm)
  }

  if (kind == "block_ecm") {
    expect_table_close(new$estimates, old$estimates, 1e-6, paste(label, "estimates"))
    expect_lte(abs(new$logLik - old$logLik), 1e-8)
    expect_identical(new$algorithm, old$algorithm)
    expect_identical(new$equations, old$equations)
  }

  if (kind == "block_gibbs") {
    expect_draws_match(new, old, label)
    expect_identical(new$algorithm, old$algorithm)
  }

  if (kind == "block_vi") {
    expect_table_close(new$estimates, old$estimates, 1e-6, paste(label, "estimates"))
    expect_lte(abs(new$elbo - old$elbo), 1e-8)
    expect_identical(new$algorithm, old$algorithm)
  }
}

for (case_name in names(regression_cases())) {
  test_that(paste("regression:", case_name), {
    skip_on_cran()
    skip_on_ci()
    skip_if(is.null(baseline), "no regression baseline")
    case <- regression_cases()[[case_name]]
    old <- baseline$cases[[case_name]]

    via_old <- suppressWarnings(regression_snapshot(case, "old"))
    expect_snapshot_match(via_old, old, case$kind, paste(case_name, "(old interface)"))

    if (!is.null(case$fit_new)) {
      via_new <- regression_snapshot(case, "new")
      expect_snapshot_match(via_new, old, case$kind, paste(case_name, "(new interface)"))
    }
  })
}
