# Generate monthly recruitment pattern
GenMonthlyRec <- function(mu=6, sigma=2) {
  tt <- dnorm(-5:6, 0, sigma)
  tt <- tt/sum(tt)
  t2 <- rep(tt, 5)
  rec_pattern <- rep(0, 12)
  n <- 12-mu
  rec_pattern[mu:(mu+n)] <- t2[18:(18+n)]
  if (mu>1)
    rec_pattern[1:(mu-1)] <- t2[(18-mu+1):17]
  rec_pattern
}

# Calculate selectivity at length
calSelL <- function(L, LF5, LFS, Vmaxlen, Linf) {

  sigma_asc <- (LFS-LF5)/sqrt(-log(0.05,2))
  sigma_dec <- (Linf-LFS)/sqrt(-log(Vmaxlen,2))

  sl <- rep(0, length(L))
  ind1 <- which(L<=LFS)
  ind2 <- which(L>LFS)
  sl[ind1] <-2^(-((L[ind1]-LFS)/sigma_asc)^2)
  sl[ind2] <- 2^(-((L[ind2]-LFS)/sigma_dec)^2)
  sl
}


#' Title
#'
#' @param scen
#'
#' @return
#' @export
#'
#' @examples
Load_Scenario <- function(scen=1) {
  Scenario_Parameters <- as.matrix(Scenario_Parameters)
  vars <- Scenario_Parameters[,1]
  vals <- Scenario_Parameters[,scen+1]
  List <- list()
  for (i in 1:length(vars)) {
    List[[vars[i]]] <- as.numeric(vals[i])
  }
  List
}


h2CR <- function(h) {
  (4*h)/(1-h)
}

CR2h <- function(CR) {
  CR/(CR+4)
}

BH_SRR <- function(R0, h, SB, SBpR) {
  (4*R0*h*SB)/(SBpR*R0*(1-h)+(5*h-1)*SB)
}





