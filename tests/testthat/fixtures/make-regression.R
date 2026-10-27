# Make the regression baseline for test-regression.R.
# Run from the package root, on the computer where the tests will be run, before any change:
#   LANG=en_US.UTF-8 Rscript tests/testthat/fixtures/make-regression.R

pkgload::load_all(".", quiet = TRUE)
invisible(testthat::source_test_helpers("tests/testthat", env = environment()))

commit <- system("git rev-parse --short HEAD", intern = TRUE)
cases <- regression_cases()
snapshots <- list()
for (case_name in names(cases)) {
  message("snapshot: ", case_name)
  snapshots[[case_name]] <- suppressWarnings(regression_snapshot(cases[[case_name]]))
}

baseline <- list(
  commit = commit,
  made = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  r_version = R.version.string,
  platform = R.version$platform,
  rng = RNGkind(),
  packages = vapply(c("Rcpp", "RcppArmadillo", "coda"), function(p) as.character(utils::packageVersion(p)), ""),
  cases = snapshots)

saveRDS(baseline, file.path("tests", "testthat", "fixtures", "regression.rds"))
message("saved ", length(snapshots), " cases (commit ", commit, ")")
