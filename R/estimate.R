#' Title
#'
#' @param data
#' @param parameters
#' @param random
#'
#' @return
#' @export
#'
#' @examples
EstModel <- function(data, parameters, random, starts, bounds=NULL) {

  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, random = random)

  step.min <- 1
  step.max <- 1
  control <- list(eval.max=1e4, iter.max=1e4,
                  step.min=step.min, step.max=step.max,
                  trace=0, abs.tol=1e-20)

  if (is.null(bounds)) {
    opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr,
                                   control=control))
  } else {
    opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr,
                                   lower=bounds[[1]], upper=bounds[[2]],
                                   control=control))
  }

 list(obj=obj, opt=opt)
}
