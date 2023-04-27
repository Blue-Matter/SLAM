#' Calculate the optimal pattern in fishing mortality
#'
#' Uses a HARA model to calculate the optimal exploitation pattern
#'
#' @param R0_m The expected fraction of recruitment in each calendar month
#' @param steepness The assumed steepness of the Beverton-Holt stock-recruitment function
#' @param Weight_Age_Mean The mean weight-at-age
#' @param Maturit_at_Age The maturity-at-age schedule
#' @param M_at_Age The natural mortality-at-age schedule
#' @param Post_Spawing_Mortality The post-spawning mortality at age schedule
#' @param sel_at_age The selectivity-at-age schedule
#' @param opt_type 0 to maximize yield, 1 to maximize HARA utility
#' @param utilpow Exponent of the HARA model
#'
#' @return A list
#' @export
#'
calculate_optimal_fishing <- function(R0_m, steepness,
                                      Weight_Age_Mean,
                                      Maturity_at_Age,
                                      M_at_Age,
                                      Post_Spawning_Mortality,
                                      sel_at_age,
                                      opt_type=1, utilpow=0.3) {

  data <- list(model='optF',
               rec_pattern=R0_m,
               opt_type=opt_type,
               utilpow=utilpow,
               h=steepness,
               Wght_Age=Weight_Age_Mean,
               Mat_at_Age=Maturity_at_Age,
               M_at_Age=M_at_Age,
               PSM_at_Age=Post_Spawning_Mortality,
               selA=sel_at_age)

  parameters <- list(logF_m=rep(log(0.01),12))
  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE)

  starts <- obj$par
  opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr))
  opt_rep <- obj$report()
  opt_rep

}


