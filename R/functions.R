# Generate monthly recruitment pattern
GenMonthlyRec <- function(mu=6, sigma=2) {
  tt <- dnorm(-5:6, 0, sigma)
  tt <- tt/sum(tt)
  t2 <- rep(tt, 5)
  rec_pattern <- rep(0, 12)
  n <- 12-mu
  rec_pattern[mu:(mu+n)] <- t2[18:(18+n)]
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
  vars <- Scenario_Parameters[,1]
  vals <- Scenario_Parameters[,scen+1]
  List <- list()
  for (i in 1:length(vars)) {
    List[[vars[i]]] <- as.numeric(vals[i])
  }
  List
}