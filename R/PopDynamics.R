



#' Title
#'
#' @param Pars
#'
#' @return
#' @export
#'
#' @examples
Simulate <- function(Pars) {

  Ages <- 0:Pars$maxage
  nAge <- length(Ages)
  LatAge <- Pars$Linf*(1-exp(-Pars$K*(Ages-Pars$t0)))
  LatAge[LatAge==0] <- 0.001
  LenSD <- LatAge * Pars$LenCV

  WatAge <- Pars$wa*LatAge^Pars$wb

  # Probability of length-at-age
  Len_Bins <- seq(0, to=Pars$BinMax, by=Pars$BinWidth)
  nbins <- length(Len_Bins)-1
  By <- Len_Bins[2]-Len_Bins[1]
  Len_Mids <- seq(Len_Bins[1]+0.5*By, by=By, length.out=nbins)
  ALK <- matrix(0, nrow=nAge, ncol=nbins)
  ALK[,1] <- pnorm((Len_Bins[2] - LatAge)/LenSD, 0, 1) # probability of length-at-age
  for (i in 2:(nbins-1)) {
    ALK[,i] <- pnorm((Len_Bins[i+1] - LatAge)/LenSD, 0, 1) -
      pnorm((Len_Bins[i] - LatAge)/LenSD, 0, 1)
  }
  ALK[,nbins] <- 1 - pnorm((Len_Bins[nbins] - LatAge)/LenSD, 0, 1)

  # Maturity
  pA <- 1/(1+exp(-log(19)*((Ages-Pars$A50)/(Pars$A95-Pars$A50))))
  # pL <- 1/(1+exp(-log(19)*((Len_Mids-Pars$L50)/(Pars$L95-Pars$L50))))

  # Selectivity - double-normal
  sL <- calSelL(Len_Mids, Pars$SL5, Pars$SFS, Pars$Vmaxlen, Pars$Linf)
  sA <- apply(ALK%*%sL, 1, sum)

  # Set up arrays
  nTS <- Pars$nyears * 12 # total number of time-steps
  N <- matrix(0, nrow=nAge, ncol=nTS)
  C <- N
  SB <- rep(0, nTS)

  # Generate recruitment deviations
  rec_sd <- Pars$sigmaR
  rec_devs <- exp(rnorm(nTS, -0.5*rec_sd^2, rec_sd))
  # monthly pattern in recruitment
  if (!is.null(Pars$rec_pattern)) {
    rec_pattern <- rep(Pars$rec_pattern, nTS)
  } else {
    rec_pattern <- GenMonthlyRec(mu=Pars$rec_mu, sigma=Pars$rec_sd)
  }


  Rbar <- Pars$Rbar # mean annual recruitment

  M_at_age <- rep(Pars$M, nAge)
  phi_at_age <- pA # proportion that die after spawning - assumed to follow maturity curve

  # Fishing mortality
  if (is.null(Pars$Fvector)) {
    # generate F vector
    Fvector <- Pars$F_mu * exp(rnorm(nTS, -0.5*Pars$F_sd^2, Pars$F_sd))
  } else {
    Fvector <- Pars$Fvector
  }
  Fvector[1] <- 0

  F_at_age <-  t(replicate(nAge, Fvector)) * sA
  Z_at_age <- M_at_age + F_at_age

  # SPR
  theta_0 <- rep(1, nAge)
  theta_F <- matrix(0, nrow=nAge, ncol=nTS)
  theta_F[1,] <- 1
  for (a in 2:nAge) {
    theta_0[a] <- theta_0[a-1]*exp(-M_at_age[a-1])*(1-phi_at_age[a-1])
  }
  for (t in 1:nTS) {
    for (a in 2:nAge) {
      theta_F[a,t] <- theta_F[a-1,t]*exp(-Z_at_age[a-1,t])*(1-phi_at_age[a-1])
    }
  }

  E0 <- sum(theta_0*pA*WatAge)
  Ef <- apply(theta_F*pA*WatAge, 2, sum)
  SPR <- Ef/E0

  # Recruitment
  if (!is.null(Pars$h)) {
    h <- Pars$h
  } else {
    h <- 0.999
  }

  # Initial time-step - unfished
  for (t in 1:36) {
    month <- t%%12
    if (month==0) month <-12
    for (a in 0:Pars$maxage) {
      if (a ==0)
        N[a+1,month] <- rec_pattern[month]*Rbar* rec_devs[t]
      if (a>0) {
        if (month==1) {
          N[a+1, month] <- N[a,12]*exp(-M_at_age[a])*(1-phi_at_age[a])
        } else {
          N[a+1, month] <- N[a,month-1]*exp(-M_at_age[a])*(1-phi_at_age[a])
        }
      }
    }
    SB[month] <- sum(N[,month]*WatAge*pA)*exp(-Fvector[month]/2)
  }

  # Loop over time-steps
  for (t in 13:nTS) {
    month <- t%%12
    if (month==0) month <-12
    for (a in 1:Pars$maxage) {
      N[a+1, t] <- N[a,t-1]*exp(-Z_at_age[a,t-1])*(1-phi_at_age[a])
    }
    #SB
    SB[t] <- sum(N[,t]*WatAge*pA)*exp(-Fvector[t]/2)
    # recruits
    N[1,t] <- BH_SRR(R0=Pars$Rbar*rec_pattern[month], h, SB=SB[t], SBpR=E0) * rec_devs[t] # rec_pattern[month]*rec_devs[t] *  Rbar
    C[,t] <- N[,t]*((1-pA)*exp(-M_at_age/2)+pA*exp(-phi_at_age/2))*(1-exp(-F_at_age[,t]))
  }

  B <- apply(N * replicate(nTS,WatAge), 2, sum)
  CB <- apply(C * replicate(nTS,WatAge), 2, sum)


  # Index of Abundance
  I <- B * exp(rnorm(nTS, -0.5*Pars$sigmaI^2, Pars$sigmaI))
  I <- I/mean(I)

  # Size structure of catch
  ALK_C <- array(0, dim=dim(ALK))
  for (a in 1:nAge) {
    ALK_C[a,] <- ALK[a,] * sL
  }
  ALK_C <- ALK_C/apply(ALK_C,1, sum)

  CAL_exp <- matrix(0, nrow=nbins, ncol=nTS)
  for (t in 1:nTS) {
    CAL_exp[,t] <- apply(ALK_C * C[,t], 2, sum)
  }
  CAL_samp <- matrix(0, nrow=nbins, ncol=nTS)
  for (t in 1:nTS) {
    if (sum(CAL_exp[,t])>0)
      CAL_samp[,t] <- Pars$CAL_nsamp * (rmultinom(1, size=Pars$CAL_ESS, prob=CAL_exp[,t]))/Pars$CAL_ESS
  }



  Effort <- Fvector/mean(Fvector) *  exp(rnorm(nTS, -0.5*Pars$sigmaE^2, Pars$sigmaE))

  out <- list()
  out$E0 <- E0
  out$Ef <- Ef
  out$Number <- N
  out$Biomass <- B
  out$SB <- SB
  out$Catch <- C
  out$Catch_Biomass <- CB
  out$CAL_exp <- CAL_exp
  out$CAL_samp <- CAL_samp
  out$Index <- I
  out$SPR <- SPR
  out$Len_Bins <- Len_Bins
  out$Len_Mids <- Len_Mids
  out$M_at_Age <- M_at_age
  out$Len_at_Age <- LatAge
  out$SD_Len_at_Age <- LenSD
  out$Wt_at_Age <- WatAge
  out$Mat_at_Age <- pA
  out$phi_at_Age <- phi_at_age
  out$F_m <- Fvector
  out$Effort <- Effort
  out$Rec_Pattern <- rec_pattern
  out$Sel_at_Age <- sA
  out$Sel_at_Length <- sL
  out$Z_at_age <- Z_at_age
  out
}




Predict <- function() {

}


