# Julia names (ExtendedRtIrtModeling.jl) as aliases of the snake_case functions, so code written
# for the Julia package can be translated line by line. Object fields are snake_case
# (fit$cond, fit$post$mean, ...).

#' Julia-compatible aliases
#'
#' The function names of ExtendedRtIrtModeling.jl, mapped to the snake_case functions:
#'
#' | Julia | birtRcpp |
#' |---|---|
#' | `setCond(nSubj, nItem, nFeat, nIter, nChain, nBurnin, nThin, nRep, qRa, qRt)` | [set_cond()] |
#' | `InputData(Y, T, X)`, `InputData4R(Y, kappa, logT, X)` | [input_data()] (`resp`, `time`, `cov`, `log_time`) |
#' | `InputPara(...)` | [input_para()] |
#' | `GibbsMlIrt`, `GibbsRtIrtNull`, `GibbsRtIrt` | [mlirt()], [rtirt_null()], [rtirt_latreg()] |
#' | `GibbsRtIrtLatent`, `GibbsRtIrtLatentQr` | [rtirt_latent()], `quantile = NULL` / `Cond$q_rt` |
#' | `GibbsRtIrtCross`, `GibbsRtIrtCrossQr` | [rtirt_cross()], `quantile = NULL` / `Cond$q_rt` |
#' | `sample!(MCMC)` | `sampleGibbs(MCMC)` = [gibbs()] |
#' | `coef(MCMC)`, `precis(MCMC)` | [summary()][summary.rtirt], [precis()] |
#' | `getDic`, `checkConvergence`, `comparePara` | [dic()], [convergence()], [compare_para()] |
#' | `getBias`, `getRmse`, `getCorr` | [get_bias()], [get_rmse()], [get_corr()] |
#' | `setTruePara<Model>`, `setData<Model>` | [sim_para()], [sim_data()] |
#' | `runSimulation(Cond, truePara; funcData, funcGibbs)` | [run_simulation()] |
#'
#' The constructors take `(Cond, Data, truePara, Para)` as in Julia and keep the Julia
#' defaults: the speed variance is fixed at 1 (`speed_var = "fixed"`) in `GibbsRtIrtNull`,
#' `GibbsRtIrt` and `GibbsRtIrtCross` (`sampleGibbs(M, cov2one = FALSE)` frees it), and the
#' `Qr` versions use the quantile level `Cond$q_rt`. `setCond()` keeps the Julia default
#' dimensions (2000 persons, 15 items, 3 covariates), which must match the data.
#' @param Cond,Data,truePara,Para,MCMC,par,est,true,... As in the Julia package.
#' @param nSubj,nItem,nFeat,nIter,nChain,nBurnin,nThin,nRep,qRa,qRt As in `setCond()`.
#' @param cov2one Julia's setting: `TRUE` fixes the variance of speed at 1 (`speed_var = "fixed"`), `FALSE` frees it.
#' @param Y,T,X,kappa,logT Data, as in `InputData()` / `InputData4R()`.
#' @param trueCorr,type,funcData,funcGibbs As in the Julia functions.
#' @name julia_aliases
NULL

#' @rdname julia_aliases
#' @export
setCond <- function(nSubj = 2000, nItem = 15, nFeat = 3, nIter = 5000, nChain = 4, nBurnin = NULL, nThin = 1,
                    nRep = 10, qRa = 0.5, qRt = 0.5)
  set_cond(nSubj, nItem, nFeat, nIter, nChain, nBurnin, nThin, nRep, qRa, qRt)
#' @rdname julia_aliases
#' @export
InputData <- function(Y, T = NULL, X = NULL) input_data(resp = Y, time = T, cov = X)
#' @rdname julia_aliases
#' @export
InputData4R <- function(Y, kappa = Y - 0.5, logT = NULL, X = NULL) input_data(resp = Y, cov = X, log_time = logT)
#' @rdname julia_aliases
#' @export
InputPara <- function(...) input_para(...)

tp_of <- function(truePara, Data) if (is.null(truePara)) attr(Data, "true_para") else truePara
#' @rdname julia_aliases
#' @export
GibbsMlIrt <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("mlirt", Data, Cond, tp_of(truePara, Data), Para)
#' @rdname julia_aliases
#' @export
GibbsRtIrtNull <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("rtirt_null", Data, Cond, tp_of(truePara, Data), Para, speed_var = "fixed")
#' @rdname julia_aliases
#' @export
GibbsRtIrt <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("rtirt_latreg", Data, Cond, tp_of(truePara, Data), Para, speed_var = "fixed")
#' @rdname julia_aliases
#' @export
GibbsRtIrtLatent <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("rtirt_latent", Data, Cond, tp_of(truePara, Data), Para)
#' @rdname julia_aliases
#' @export
GibbsRtIrtLatentQr <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("rtirt_latent", Data, Cond, tp_of(truePara, Data), Para, quantile = Cond$q_rt)
#' @rdname julia_aliases
#' @export
GibbsRtIrtCross <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("rtirt_cross", Data, Cond, tp_of(truePara, Data), Para, speed_var = "fixed")
#' @rdname julia_aliases
#' @export
GibbsRtIrtCrossQr <- function(Cond, Data, truePara = NULL, Para = NULL)
  new_model("rtirt_cross", Data, Cond, tp_of(truePara, Data), Para, quantile = Cond$q_rt, speed_var = "fixed")

#' @rdname julia_aliases
#' @export
sampleGibbs <- function(MCMC, cov2one = NULL, ...) {               # Julia: cov2one = TRUE fixes Var(speed) = 1
  if (!is.null(cov2one)) MCMC$settings$speed_var <- if (cov2one) "fixed" else "free"
  gibbs(MCMC, ...)
}
#' @rdname julia_aliases
#' @export
getDic <- function(MCMC) dic(MCMC)
#' @rdname julia_aliases
#' @export
checkConvergence <- function(MCMC, ...) convergence(MCMC, ...)
#' @rdname julia_aliases
#' @export
comparePara <- function(MCMC, par = "a") compare_para(MCMC, par)
#' @rdname julia_aliases
#' @export
getBias <- function(est, true) get_bias(est, true)
#' @rdname julia_aliases
#' @export
getRmse <- function(est, true) get_rmse(est, true)
#' @rdname julia_aliases
#' @export
getCorr <- function(est, true) get_corr(est, true)

#' @rdname julia_aliases
#' @export
setTrueParaMlIrt <- function(Cond) sim_para(Cond, "mlirt")
#' @rdname julia_aliases
#' @export
setTrueParaRtIrt <- function(Cond, trueCorr = 0.3) sim_para(Cond, "latreg", true_corr = trueCorr)
#' @rdname julia_aliases
#' @export
setTrueParaRtIrtLatent <- function(Cond) sim_para(Cond, "latent")
#' @rdname julia_aliases
#' @export
setTrueParaRtIrtCross <- function(Cond) sim_para(Cond, "cross")
#' @rdname julia_aliases
#' @export
setDataMlIrt <- function(Cond, truePara) sim_data(Cond, truePara, "mlirt")
#' @rdname julia_aliases
#' @export
setDataRtIrt <- function(Cond, truePara) sim_data(Cond, truePara, "latreg")
#' @rdname julia_aliases
#' @export
setDataRtIrtNull <- function(Cond, truePara) { truePara$beta <- NULL; sim_data(Cond, truePara, "null") }
#' @rdname julia_aliases
#' @export
setDataRtIrtLatent <- function(Cond, truePara, type = "norm") sim_data(Cond, truePara, "latent", type)
#' @rdname julia_aliases
#' @export
setDataRtIrtCross <- function(Cond, truePara, type = "norm") sim_data(Cond, truePara, "cross", type)

#' @rdname julia_aliases
#' @export
runSimulation <- function(Cond, truePara, funcData = setDataRtIrt, funcGibbs = GibbsRtIrt, par = c("a", "b", "theta"), ...) {
  # the Julia data generator and sampler -> run_simulation(model, quantile, speed_var)
  jl <- list(list(setDataMlIrt, GibbsMlIrt, "mlirt", FALSE, "free"), list(setDataRtIrtNull, GibbsRtIrtNull, "null", FALSE, "fixed"),
             list(setDataRtIrt, GibbsRtIrt, "latreg", FALSE, "fixed"), list(setDataRtIrtLatent, GibbsRtIrtLatent, "latent", FALSE, "free"),
             list(setDataRtIrtLatent, GibbsRtIrtLatentQr, "latent", TRUE, "free"), list(setDataRtIrtCross, GibbsRtIrtCross, "cross", FALSE, "fixed"),
             list(setDataRtIrtCross, GibbsRtIrtCrossQr, "cross", TRUE, "fixed"))
  hit <- Filter(function(r) identical(funcData, r[[1]]) && identical(funcGibbs, r[[2]]), jl)
  if (!length(hit)) stop("runSimulation(): funcData and funcGibbs must be a matching pair of the Julia functions (e.g. setDataRtIrtCross, GibbsRtIrtCross); for other designs use run_simulation()", call. = FALSE)
  r <- hit[[1]]
  run_simulation(Cond, truePara, r[[3]], quantile = if (r[[4]]) Cond$q_rt, pars = par, speed_var = r[[5]], ...)
}
