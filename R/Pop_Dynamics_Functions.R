

#' Title
#'
#' @param R0
#' @param h
#' @param SB
#' @param SBpR
#'
#' @return
#' @export
BH_SRR <- function(R0, h, SB, SBpR) {
  (4*R0*h*SB)/(SBpR*R0*(1-h)+(5*h-1)*SB)
}

