#' birtRcpp: BIvariate Response Time Modeling with 'Rcpp'
#'
#' Built-in models of response accuracy and response times (ports of ExtendedRtIrtModeling.jl),
#' named by `model` and fitted with [ecm()] (the main estimator) or [gibbs()]; [algorithm()] and
#' [equations()] show the steps and the model, [check_ecm()] and [check_sampler()] check them:
#' ```
#' dat <- input_data(resp = items, time = times, cov = covs, id = "IDSTUD", data = raw)
#' fit_ecm <- ecm(dat, model = "cross")
#' fit <- gibbs(dat, model = "cross", init = fit_ecm, n_iter = 10000, n_chain = 3, seed = 1)
#' summary(fit)
#' estimates(fit)
#' scores(fit)
#' ```
#' The declarative [block_model()], [vi()] and [mml()] are deprecated (they still work).
#' The accessors have the names of the birt package. Julia names are aliases ([julia_aliases]).
#' The generics shared with birt (`estimates`, `scores`, `reliability`, `fit_indices`, `equations`,
#' `convergence`, `score_test`, `algorithm`) dispatch to the methods of both packages whichever was attached last: when
#' both are loaded, each package's methods are registered in the other's generics.
#' @useDynLib birtRcpp, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @rawNamespace export(estimates.rtirt, estimates.rtirt_qset, scores.rtirt, reliability.rtirt, fit_indices.rtirt, convergence.rtirt, score_test.rtirt)
#' @keywords internal
"_PACKAGE"

shared_generics <- c("estimates", "scores", "reliability", "fit_indices", "convergence", "equations", "score_test", "algorithm")

# register the methods of birtRcpp in birt's generics and vice versa (S3 lookup from the global
# environment skips attached packages, so exported methods alone are not found)
bridge_birt <- function(...) {
  if (!isNamespaceLoaded("birt")) return(invisible())
  nb <- asNamespace("birt"); me <- asNamespace("birtRcpp")
  tb <- get(".__S3MethodsTable__.", envir = nb); tm <- get(".__S3MethodsTable__.", envir = me)
  copy <- function(g, from, to, env) for (n in ls(from, pattern = paste0("^", g, "\\.")))
    if (!exists(n, envir = to, inherits = FALSE)) registerS3method(g, substring(n, nchar(g) + 2), get(n, envir = from), envir = env)
  for (g in intersect(shared_generics, getNamespaceExports("birt"))) { copy(g, tm, tb, nb); copy(g, tb, tm, me) }
  invisible()
}

.onLoad <- function(libname, pkgname) {
  setHook(packageEvent("birt", "onLoad"), bridge_birt)
  bridge_birt()
}
