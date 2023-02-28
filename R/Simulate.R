
#' Simulate a fishery for a short-lived semelparous species
#'
#' @param LifeHistory A `list` object with the life-history parameters
#' @param Explotation A `list` object with the exploitation parameters
#' @param Data A `list` object with the parameters for generating data
#' @param nsim The number of stochastic simulations
#' @param seed The seed for the random number generator
#' @param currentYr The current (most recent) year for the simulations
#'
#' @export
Simulate <- function(LifeHistory, Exploitation, Data, nsim=3, seed=101,
                     currentYr=format(Sys.Date(), "%Y")) {

  set.seed(seed)

  # Import Life-History Parameters
  Ages <- LifeHistory$Ages
  maxage <- max(Ages)
  nAge <- length(Ages)
  Weight_Age_Mean <- LifeHistory$Weight_Age_Mean
  Weight_Age_SD <- LifeHistory$Weight_Age_SD

  M_at_Age <- LifeHistory$M_at_Age
  Maturity_at_Age <- LifeHistory$Maturity_at_Age
  Post_Spawning_Mortality <- LifeHistory$Post_Spawning_Mortality

  R0_bar <- LifeHistory$R0_bar
  R0_m <- LifeHistory$R0_m
  R0_m <- R0_m/sum(R0_m)
  sigmaR <- LifeHistory$sigmaR
  h <- LifeHistory$steepness
  if (is.null(h)) h <- 0.999
  if (is.na(h)) h <- 0.999

  # Import Exploitation Parameters
  SA50 <- Exploitation$SA50
  SA95 <- Exploitation$SA95
  Effort_Annual_Mean <- Exploitation$Effort_Annual_Mean
  Effort_Annual_SD <- Exploitation$Effort_Annual_SD
  Effort_Month_Mean <- Exploitation$Effort_Month_Mean
  Effort_Month_SD <- Exploitation$Effort_Month_SD
  nyears <- Exploitation$nyears
  nts <- nyears * 12
  q <- Exploitation$q

  # Import Data Parameters
  n_recent_months <- Data$n_recent_months
  Rel_Sample_Month <- Data$Rel_Sample_Month
  CPUE_CV <- Data$CPUE_CV
  Catch_CV <- Data$Catch_CV
  Effort_CV <- Data$Effort_CV

  Weight_Bin_Max <- Data$Weight_Bin_Max
  Weight_Bin_Width <- Data$Weight_Bin_Width
  Weight_Bins <- seq(0, Weight_Bin_Max, Weight_Bin_Width)
  nBins <- length(Weight_Bins)-1
  Weight_Mids <- seq(0.5*Weight_Bin_Width, by=Weight_Bin_Width, length.out=nBins)

  Data$Weight_Bins <- Weight_Bins
  Data$Weight_Mids <- Weight_Mids

  CAW_Annual_Sample_Size <- Data$CAW_Annual_Sample_Size
  CAW_Annual_ESS <- Data$CAW_Annual_ESS

  # Parameter Checks
  #TODO

  # ---- Model Fishery Dynamics ----

  # Generate Monthly Recruitment Deviations
  Rec_Devs <- exp(rnorm((nts+maxage)*nsim, -0.5*sigmaR^2, sigmaR))
  Rec_Devs <- matrix(Rec_Devs, nrow=nts+maxage, ncol=nsim)
  Rec_Devs <- t(Rec_Devs)

  # Generate Fishing Effort and Mortality
  Effort_Annual_Deviations <- rlnorm(nyears*nsim, -0.5*Effort_Annual_SD^2, Effort_Annual_SD)
  Effort_Annual_Deviations <- matrix(Effort_Annual_Deviations, nrow=nyears, ncol=nsim)
  Effort_Annual <- Effort_Annual_Mean * Effort_Annual_Deviations

  Effort_Month_Deviations <- rlnorm(nts*nsim, -0.5*Effort_Month_SD^2, Effort_Month_SD)
  Effort_Month_Deviations <- matrix(Effort_Month_Deviations, nrow=nts, ncol=nsim)

  Effort_Month_Mean <- matrix(rep(Exploitation$Effort_Month_Mean, nyears),
                              nrow=nts, ncol=nsim, byrow=FALSE)

  Effort_Month <- matrix(rep(Effort_Annual, each=12), nrow=nts, ncol=nsim) * Effort_Month_Mean * Effort_Month_Deviations
  Effort_Month <- t(Effort_Month)
  F_Month <- Effort_Month * q

  # Set up Arrays
  M_at_Age <- replicate(nsim, LifeHistory$M_at_Age)
  M_at_Age <- replicate(nts, M_at_Age)
  M_at_Age <- aperm(M_at_Age, c(2,1,3))
  F_at_Age <- Z_at_Age <- array(0, dim=dim(M_at_Age))

  N_Age_fished <- array(0, dim=c(nsim, nAge, nts))
  N_Age_unfished <- N_Age_unfished_eq <- SPR <- SB_Age_fished <- SB_Age_unfished <-
    SB_Age_unfished_eq <- Catch_Age <- N_Age_fished
  SB_fished <- SB_unfished <- SB_unfished_eq <- matrix(0, nrow=nsim, ncol=nts)

  Sel_at_Age <-  1/(1+exp(-log(19)*((Ages-SA50)/(SA95-SA50))))
  Sel_at_Age <- replicate(nsim, Sel_at_Age)
  Sel_at_Age <- replicate(nts, Sel_at_Age)
  Sel_at_Age <- aperm(Sel_at_Age, c(2,1,3))

  ind <- as.matrix(expand.grid(1:nsim, 1:nAge, 1:nts))
  F_at_Age[ind] <- F_Month[ind[,c(1,3)]] * Sel_at_Age[ind]
  Z_at_Age <- M_at_Age + F_at_Age

  # Calculate Spawning Potential Ratio
  theta_0 <- array(1, dim=c(nsim, nAge))
  theta_F <- array(0, dim=c(nsim, nAge,nts))
  theta_F[,1,] <- 1
  for (a in 2:nAge) {
    theta_0[,a] <- theta_0[,a-1]*exp(-M_at_Age[a-1])*(1-Post_Spawning_Mortality[a-1])
  }
  for (t in 1:nts) {
    for (a in 2:nAge) {
      theta_F[,a,t] <- theta_F[,a-1,t]*exp(-(Z_at_Age[,a-1,t]))*(1-Post_Spawning_Mortality[a-1])
    }
  }

  Maturity_at_Age_array <- t(replicate(nsim, Maturity_at_Age))
  Weight_Age_Mean_array <- t(replicate(nsim, Weight_Age_Mean))
  Post_Spawning_Mortality_array <- t(replicate(nsim, Post_Spawning_Mortality))
  E0 <- apply(theta_0*Maturity_at_Age_array*Weight_Age_Mean_array, 1, sum)
  E0 <- replicate(nts, E0)

  Maturity_at_Age_array <- replicate(nts, Maturity_at_Age_array)
  Weight_Age_Mean_array <- replicate(nts, Weight_Age_Mean_array)
  Post_Spawning_Mortality_array <- replicate(nts, Post_Spawning_Mortality_array)
  Ef <- apply(theta_F*Maturity_at_Age_array*Weight_Age_Mean_array,c(1,3), sum)
  SPR <- Ef/E0

  # Initialize unfished population
  N_Age_unfished[,1,1] <-  R0_bar * R0_m[1] * Rec_Devs[,1]
  N_Age_unfished_eq[,1,1] <-  R0_bar * R0_m[1]
  for (a in 2:nAge) {
    age <- a -1
    if (age %% 12==0) {
      m <- -(age%%12)+1
    } else {
      m <- -(age %% 12)+13
    }

    cumM <- apply(M_at_Age[,1:(a-1),1, drop=FALSE], 1, sum)
    cumPSM <- prod(1-Post_Spawning_Mortality[1:(a-1)])

    N_Age_unfished[,a,1] <-  R0_bar * R0_m[m]* exp(-cumM) * cumPSM  * Rec_Devs[,a]
    N_Age_unfished_eq[,a,1] <-  R0_bar * R0_m[m]* exp(-cumM) * cumPSM
  }

  N_Age_fished <- N_Age_unfished
  SB_Age_fished[,,1] <- N_Age_fished[,,1] * Weight_Age_Mean_array[,,1] * Maturity_at_Age_array[,,1] * exp(-F_Month[,1]/2)
  SB_Age_unfished[,,1] <- N_Age_unfished_eq[,,1] * Weight_Age_Mean_array[,,1] * Maturity_at_Age_array[,,1]
  SB_Age_unfished_eq[,,1] <- N_Age_unfished_eq[,,1] * Weight_Age_Mean_array[,,1] * Maturity_at_Age_array[,,1]
  SB_fished[,1] <- apply(SB_Age_fished[,,1], 1, sum)
  SB_unfished[,1] <- apply(SB_Age_unfished[,,1], 1, sum)
  SB_unfished_eq[,1] <- apply(SB_Age_unfished_eq[,,1], 1, sum)
  Catch_Age[,,1] <- N_Age_fished[,,1]*((1-Maturity_at_Age_array[,,1])*exp(-M_at_Age[,,1]/2)+Maturity_at_Age_array[,,1]*exp(-Post_Spawning_Mortality_array[,,1]/2))*(1-exp(-F_at_Age[,,1]))

  # Fished population - Loop over remaining time-steps
  for (t in 2:nts) {
    month <- t%%12
    if (month==0) month <-12
    for (a in 1:maxage) {
      N_Age_fished[,a+1, t] <- N_Age_fished[,a,t-1]*exp(-Z_at_Age[,a,t-1])*(1-Post_Spawning_Mortality_array[,a,t-1])
      N_Age_unfished[,a+1, t] <- N_Age_unfished[,a,t-1]*exp(-M_at_Age[,a,t-1])*(1-Post_Spawning_Mortality_array[,a,t-1])
      N_Age_unfished_eq[,a+1, t] <- N_Age_unfished_eq[,a,t-1]*exp(-M_at_Age[,a,t-1])*(1-Post_Spawning_Mortality_array[,a,t-1])
    }

    # Spawning Biomass mid time-step
    SB_Age_fished[,,t] <- N_Age_fished[,,t] * Weight_Age_Mean_array[,,t] * Maturity_at_Age_array[,,t] * exp(-F_Month[,t]/2)
    SB_Age_unfished[,,t] <- N_Age_unfished[,,t] * Weight_Age_Mean_array[,,t] * Maturity_at_Age_array[,,t]
    SB_Age_unfished_eq[,,t] <- N_Age_unfished_eq[,,t] * Weight_Age_Mean_array[,,t] * Maturity_at_Age_array[,,t]
    SB_fished[,t] <- apply(SB_Age_fished[,,t], 1, sum)
    SB_unfished[,t] <- apply(SB_Age_unfished[,,t], 1, sum)
    SB_unfished_eq[,t] <- apply(SB_Age_unfished_eq[,,t], 1, sum)

    # Catch during time-step
    Catch_Age[,,t] <- N_Age_fished[,,t]*((1-Maturity_at_Age_array[,,t])*exp(-M_at_Age[,,t]/2)+Maturity_at_Age_array[,,t]*exp(-Post_Spawning_Mortality_array[,,t]/2))*(1-exp(-F_at_Age[,,t]))

    # Recruitment to age-0 this time-step
    N_Age_fished[,1, t] <- sapply(1:nsim, function(x)
      BH_SRR(R0=R0_bar*R0_m[month], h=h, SB=SB_fished[x,t], SBpR=E0[x,t]) * Rec_Devs[x,t])

    N_Age_unfished[,1, t] <- sapply(1:nsim, function(x)
      BH_SRR(R0=R0_bar*R0_m[month], h=h, SB=SB_unfished[x,t], SBpR=E0[x,t]) * Rec_Devs[x,t])

    N_Age_unfished_eq[,1, t] <- sapply(1:nsim, function(x)
      BH_SRR(R0=R0_bar*R0_m[month], h=h, SB=SB_unfished_eq[x,t], SBpR=E0[x,t]))
  }

  # Calculate Monthly Time-series
  B_unfished_Month <- apply(N_Age_unfished * Weight_Age_Mean_array, c(1,3), sum)
  SB_unfished_Month <- apply(N_Age_unfished * Weight_Age_Mean_array * Maturity_at_Age_array, c(1,3), sum)

  B_Month <- apply(N_Age_fished * Weight_Age_Mean_array, c(1,3), sum)
  Catch_N_Month <- apply(Catch_Age, c(1,3), sum)
  Catch_B_Month <- apply(Catch_Age * Weight_Age_Mean_array, c(1,3), sum)

  N_unfished <- apply(N_Age_unfished, c(1,3), sum)
  N_fished <- apply(N_Age_fished, c(1,3), sum)

  # ---- Generate Data (for all years) ----

  # Calculate Age-Weight Key (assume log-normal variability in weight-at-age)
  AWK <- matrix(0, nrow=nAge, ncol=nBins)
  mu <- log(Weight_Age_Mean) -0.5*Weight_Age_SD^2
  AWK[,1] <- plnorm(Weight_Bins[2], mu, Weight_Age_SD) # probability of weight-at-age
  for (i in 2:(nBins-1)) {
    AWK[,i] <- plnorm(Weight_Bins[i+1], mu, Weight_Age_SD) -
      plnorm(Weight_Bins[i], mu, Weight_Age_SD)
  }
  AWK[,nBins] <- 1 - plnorm(Weight_Bins[nBins], mu, Weight_Age_SD)

  # Catch-at-Weight - expected
  CAW_exp <- array(0, dim=c(nsim, nBins, nts))
  for (ts in 1:nts) {
    CAW_exp[,,ts] <- t(sapply(1:nsim, function(x)
                 apply(AWK * Catch_Age[x,,ts], 2, sum)))
  }

  # Catch (absolute)
  Catch_Sample <- Catch_B_Month * exp(rnorm(nts*nsim, -0.5*Catch_CV^2, Catch_CV))

  # CPUE
  CPUE_Sample <- B_Month/apply(B_Month, 1, mean) * exp(rnorm(nts*nsim, -0.5*CPUE_CV^2, CPUE_CV))

  # Effort (Index)
  Effort_Sample <- Effort_Month/apply(Effort_Month, 1, mean) * exp(rnorm(nts*nsim, -0.5*Effort_CV^2, Effort_CV))


  # ---- Generate Monthly Samples for n_recent_months ----
  sample_ts <- (nts-n_recent_months+1):nts
  n_sample_ts <- length(sample_ts)
  Catch_Sample <- Catch_Sample[,sample_ts, drop=FALSE]
  CPUE_Sample <- CPUE_Sample[,sample_ts, drop=FALSE]
  CPUE_Sample <- CPUE_Sample/apply(CPUE_Sample, 1, mean)
  Effort_Sample <- Effort_Sample[,sample_ts, drop=FALSE]
  Effort_Sample <- Effort_Sample/apply(Effort_Sample, 1, mean)

  CAW_Sample <- array(0, dim=c(nsim, nBins, n_sample_ts))
  Rel_Sample_Month <- Rel_Sample_Month/sum(Rel_Sample_Month)

  for (i in seq_along(sample_ts)) {
    ts <- sample_ts[i]
    month <- ts %%12
    if (month==0) month <- 12
    ts_sample_size <- CAW_Annual_Sample_Size * Rel_Sample_Month[month]
    ts_ess <- CAW_Annual_ESS * Rel_Sample_Month[month]

    CAW_Sample[,,i] <- t(sapply(1:nsim, function(x)
      if (sum(CAW_exp[x,,ts])>0) {
        val <-  ts_sample_size * (rmultinom(1, size=ts_ess, prob=CAW_exp[x,,ts]))/ts_ess
      } else {
        val <- rep(NA, nBins)
      }
    ))
  }


  CAA_Sample <- array(0, dim=c(nsim, nAge, n_sample_ts))

  for (i in seq_along(sample_ts)) {
    ts <- sample_ts[i]
    month <- ts %%12
    if (month==0) month <- 12
    ts_sample_size <- CAW_Annual_Sample_Size * Rel_Sample_Month[month]
    ts_ess <- CAW_Annual_ESS * Rel_Sample_Month[month]

    CAA_Sample[,,i] <- t(sapply(1:nsim, function(x)
      if (sum(Catch_Age[x,,ts])>0) {
        val <-  ts_sample_size * (rmultinom(1, size=ts_ess, prob=Catch_Age[x,,ts]))/ts_ess
      } else {
        val <- rep(NA, nAge)
      }
    ))
  }


  # Return info
  currentYr <- as.numeric(currentYr)
  Years <- (currentYr-nyears+1):currentYr
  Years <- rep(Years, each=12)
  Months <- rep(month.abb, nyears)

  OM_DF <- data.frame(Sim=1:nsim,
                      Year=rep(Years, each=nsim),
                      Month=rep(Months, each=nsim),
                      Month_ind=rep(1:nts, each=nsim),
                      B_unfished=as.vector(B_unfished_Month),
                      B_fished=as.vector(B_Month),
                      SB_unfished=as.vector(SB_unfished_Month),
                      SB_unfished_eq=as.vector(SB_unfished_eq),
                      SB_fished=as.vector(SB_fished),
                      N_unfished=as.vector(N_unfished),
                      N_fished=as.vector(N_fished),
                      Recruits=as.vector(N_Age_fished[,1,]),
                      Catch=as.vector(Catch_B_Month),
                      SPR=as.vector(SPR),
                      Effort=as.vector(Effort_Month),
                      F_mort=as.vector(F_Month),
                      Rec_Devs=as.vector(Rec_Devs[,(maxage+1):ncol(Rec_Devs)]))

  # Data
  month_ind <- sample_ts %%12
  month_ind[month_ind==0] <- 12

  Months <- month.abb[month_ind]
  n_recent_years <- Months %>% table() %>% max()

  Years <- (currentYr-n_recent_years+1):currentYr

  Data_TS_DF <- data.frame(Sim=1:nsim,
                         Year=rep(Years, each=nsim),
                         Month=rep(Months, each=nsim),
                         Month_ind=rep(1:n_sample_ts, each=nsim),
                         Catch=as.vector(Catch_Sample),
                         CPUE=as.vector(CPUE_Sample),
                         Effort=as.vector(Effort_Sample))

  Data_CAW_DF <- data.frame(Sim=1:nsim,
                            Weight=rep(Weight_Mids, each=nsim),
                            Year=rep(Years, each=nsim*nBins),
                            Month=rep(Months, each=nsim*nBins),
                            Month_ind=rep(1:n_sample_ts, each=nsim*nBins),
                            Count=as.vector(CAW_Sample))

  Data_CAA_DF <- data.frame(Sim=1:nsim,
                            Weight=rep(Weight_Mids, each=nsim),
                            Year=rep(Years, each=nsim*nAge),
                            Month=rep(Months, each=nsim*nAge),
                            Month_ind=rep(1:n_sample_ts, each=nsim*nAge),
                            Count=as.vector(CAA_Sample))

  list(LifeHistory=LifeHistory,
       Exploitation=Exploitation,
       Data=Data,
       OM_DF=OM_DF,
       Data_TS_DF=Data_TS_DF,
       Data_CAW_DF=Data_CAW_DF,
       Data_CAA_DF=Data_CAA_DF)

}







