#' Title
#'
#' @param LifeHistory
#' @param Exploitation
#' @param opt_type
#' @param utilpow
#'
#' @return
#' @export
#'
calculate_optimal_fishing <- function(LifeHistory, Exploitation, opt_type=1, utilpow=0.3) {

  data <- list(model='optF',
               rec_pattern=LifeHistory$R0_m,
               opt_type=opt_type,
               utilpow=utilpow,
               h=LifeHistory$steepness,
               Wght_Age=LifeHistory$Weight_Age_Mean,
               Mat_at_Age=LifeHistory$Maturity_at_Age,
               M_at_Age=LifeHistory$M_at_Age,
               PSM_at_Age=LifeHistory$Post_Spawning_Mortalit,
               selA= 1/(1+exp(-log(19)*((LifeHistory$Ages-Exploitation$SA50)/(Exploitation$SA95-Exploitation$SA50)))))

  parameters <- list(logF_m=rep(log(0.01),12))
  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE)

  starts <- obj$par
  opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr))
  opt_rep <- obj$report()
  opt_rep

}


