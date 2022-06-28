


#' Simulate the fishery
#'
#' @param Pars A list of parameters
#'
#' @return A dataframe with simulated fishery dynamics
#' @export
#'
#' @examples
Simulate <- function(Pars, sim=1) {

  Ages <- 0:Pars$maxage
  nAge <- length(Ages)
  Weight_Age <- Pars$Weight_Age
  Weight_Age_SD <- Pars$Weight_Age_SD
  Weight_Bins <- Pars$WghtBins
  Weight_Mids <- Pars$WghtMids

  # Probability of weight-at-age
  nbins <- length(Weight_Mids)
  AWK <- matrix(0, nrow=nAge, ncol=nbins)
  mu <- log(Weight_Age) -0.5*Weight_Age_SD^2
  AWK[,1] <- plnorm(Weight_Bins[2], mu, Weight_Age_SD) # probability of weight-at-age
  for (i in 2:(nbins-1)) {
    AWK[,i] <- plnorm(Weight_Bins[i+1], mu, Weight_Age_SD) -
      plnorm(Weight_Bins[i], mu, Weight_Age_SD)
  }
  AWK[,nbins] <- 1 - plnorm(Weight_Bins[nbins], mu, Weight_Age_SD)

  # Natural mortality
  M_at_Age <- rep(Pars$M, nAge)

  # Post-spawning mortality
  PSM_at_Age <- Pars$PSM_at_Age

  # Maturity
  Mat_at_Age <- Pars$Mat_at_Age

  # Selectivity - double-normal
  Sel_at_Age <-  1/(1+exp(-log(19)*((Ages-Pars$sA50)/(Pars$sA95-Pars$sA50))))

  # recruitment deviations
  rec_devs <- Pars$rec_devs[,sim]
  rec_pattern <- Pars$rec_mu
  Rbar <- Pars$Rbar # mean annual recruitment

  # Fishing Effort and Mortality
  Eff_DF <- Pars$Effort %>% dplyr::filter(Sim==sim)
  nts <- nrow(Eff_DF)
  Effort <- Eff_DF$Effort

  Fvector <- Pars$q * Effort # fishing mortality by month


  # Set up arrays
  N <- matrix(0, nrow=nAge, ncol=nts)
  C <- N
  SB <- rep(0, nts)
  F_at_Age <-  t(replicate(nAge, Fvector)) * Sel_at_Age
  Z_at_Age <- M_at_Age + F_at_Age

  # SPR
  theta_0 <- rep(1, nAge)
  theta_F <- matrix(0, nrow=nAge, ncol=nts)
  theta_F[1,] <- 1
  for (a in 2:nAge) {
    theta_0[a] <- theta_0[a-1]*exp(-M_at_Age[a-1])*(1-PSM_at_Age[a-1])
  }
  for (t in 1:nts) {
    for (a in 2:nAge) {
      theta_F[a,t] <- theta_F[a-1,t]*exp(-Z_at_Age[a-1,t])*(1-PSM_at_Age[a-1])
    }
  }

  E0 <- sum(theta_0*pA*Weight_Age)
  Ef <- apply(theta_F*pA*Weight_Age, 2, sum)
  SPR <- Ef/E0

  # Recruitment
  if (!is.null(Pars$h)) {
    h <- Pars$h
  } else {
    h <- 0.999
  }

  # Unfished population
  N_unfished <- matrix(0, nAge, 12)

  # Initial time-step - unfished
  for (t in 1:60) {
    month <- t%%12
    if (month==0) month <-12
    for (a in 0:Pars$maxage) {
      if (a ==0)
        N_unfished[a+1,month] <- rec_pattern[month]*Rbar
      if (a>0) {
        if (month==1) {
          N_unfished[a+1, month] <- N_unfished[a,12]*exp(-M_at_Age[a])*(1-PSM_at_Age[a])
        } else {
          N_unfished[a+1, month] <- N_unfished[a,month-1]*exp(-M_at_Age[a])*(1-PSM_at_Age[a])
        }
      }
    }
  }

  SB0 <- apply(N_unfished * Weight_Age*Mat_at_Age, 2, sum)

  # Fished population - Loop over time-steps
  for (t in 1:nts) {
    month <- t%%12
    # recruits
    if (month==0) month <-12
    if (t ==1) {
      for (a in 1:Pars$maxage) {
        N[a+1, t] <- N_unfished[a,12]*exp(-Z_at_Age[a,1])*(1-PSM_at_Age[a])
      }
    } else {
      for (a in 1:Pars$maxage) {
        N[a+1, t] <- N[a,t-1]*exp(-Z_at_Age[a,t-1])*(1-PSM_at_Age[a])
      }
    }
    #SB
    SB[t] <- sum(N[,t]*Weight_Age*Mat_at_Age)*exp(-Fvector[t]/2)
    C[,t] <- N[,t]*((1-Mat_at_Age)*exp(-M_at_Age/2)+Mat_at_Age*exp(-PSM_at_Age/2))*(1-exp(-F_at_Age[,t]))
    N[1,t] <- BH_SRR(R0=Pars$Rbar*rec_pattern[month], h, SB=SB[t], SBpR=E0) * rec_devs[t] # rec_pattern[month]*rec_devs[t] *  Rbar
  }

  B <- apply(N * replicate(nts,Weight_Age), 2, sum)
  CB <- apply(C * replicate(nts,Weight_Age), 2, sum)

  # Generate Data
  data_nts <- Pars$nyears * 12
  data_ts <- (nts-data_nts+1):nts # time-steps for generating data

  # Index of Abundance
  cpue <- CB[data_ts]/Effort[data_ts] * exp(rnorm(data_nts, -0.5*Pars$sigmaI^2, Pars$sigmaI))
  cpue <- cpue/mean(cpue)

  # Index of Effort
  Eff_ind <- Effort[data_ts] * exp(rnorm(data_nts, -0.5*Pars$sigmaE^2, Pars$sigmaE))
  Eff_ind <- Eff_ind/mean(Eff_ind)

  # Catch-at-Weight
  CAW_exp <- matrix(0, nrow=nbins, ncol=nts)
  for (t in 1:nts) {
    CAW_exp[,t] <- apply(AWK * C[,t], 2, sum)
  }
  CAW_samp <- matrix(0, nrow=nbins, ncol=nts)
  for (t in 1:nts) {
    if (sum(CAW_exp[,t])>0)
      CAW_samp[,t] <- Pars$CAW_nsamp * (rmultinom(1, size=Pars$CAW_ESS, prob=CAW_exp[,t]))/Pars$CAW_ESS
  }

  out <- list()
  out$time_steps <- Pars$Effort %>% select(t, Year, Month, Date)
  out$E0 <- E0
  out$Ef <- Ef
  out$Number <- N
  out$Biomass <- B
  out$SB <- SB
  out$Catch <- C
  out$Catch_Biomass <- CB
  out$CAW_exp <- CAW_exp
  out$CAW_samp <- CAW_samp
  out$Index <- cpue
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
  out$Eff_ind <- Eff_ind
  out$Rec_Pattern <- rec_pattern
  out$Sel_at_Age <- Sel_at_Age
  out$Pars <- Pars
  out
}

Project <- function() {

}


Predict <- function() {

}


# ---- Calculate Optimal Fishing Pattern ----
Optimize <- function(Data, Rec_Pattern, selA, opt_type=1, utilpow=0.3, assumed_h=0.6) {

  data <- list(model='optF',
               rec_pattern=Rec_Pattern,
               opt_type=opt_type,
               utilpow=utilpow,
               h=assumed_h,
               Wght_Age=Data$Wt_at_Age,
               Mat_at_Age=Data$Mat_at_Age,
               M_at_Age=Data$M_at_Age,
               PSM_at_Age=Data$phi_at_Age,
               selA=selA)

  parameters <- list(logF_m=rep(log(0.01),12))
  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE)
  starts <- obj$par
  opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr))
  opt_rep <- obj$report()
  opt_rep
}
