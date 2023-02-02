Initialize_Parameters_OM <- function(SimMod,
                                     log_sigmaF=log(0.05),
                                     log_sigmaR0=log(0.1),
                                     log_sigmaR=log(0.3)) {
  parameters <- list()

  parameters$ls50 <- log(SimMod$Exploitation$SA50)
  parameters$lsdelta <- log(SimMod$Exploitation$SA95-SimMod$Exploitation$SA50)

  n_age <- SimMod$LifeHistory$maxage+1
  parameters$logF_minit <- rep(log(0.05), n_age)


  Data_Y_M <- SimMod$Data_TS_DF %>% filter(Sim==1) %>% distinct(Year, Month)
  OM_sub <-SimMod$OM_DF %>% filter(Sim==1, Year%in%Data_Y_M$Year, Month%in%Data_Y_M$Month)

  parameters$logF_m <- log(OM_sub$F_mort)
  parameters$log_sigmaF <- log_sigmaF # standard deviation for random walk penalty for F
  parameters$logR0_m_est <- log(SimMod$LifeHistory$R0_m[2:12]/mean(SimMod$LifeHistory$R0_m[2:12]))
  parameters$log_sigmaR0 <- log_sigmaR0 # sd for random walk penalty for monthly recruitment
  parameters$logRec_Devs <- log(OM_sub$Rec_Devs)
  parameters$log_sigmaR  <- log(SimMod$LifeHistory$sigmaR) # monthly rec dev sd (usually fixed)

  parameters
}

Construct_Data_OM <- function(SimMod, sim=1,
                           CAW_Monthly_ESS=100,
                           Effort_CV=0.2,
                           CPUE_CV=0.2) {

  data <- list()
  # Assumed life-history parameters
  data$Weight_Age <- SimMod$LifeHistory$Weight_Age_Mean
  data$Weight_Age_SD <- SimMod$LifeHistory$Weight_Age_SD
  data$Mat_at_Age <- SimMod$LifeHistory$Maturity_at_Age
  data$M_at_Age <- SimMod$LifeHistory$M_at_Age
  data$PSM_at_Age <- SimMod$LifeHistory$Post_Spawning_Mortality
  data$h <- SimMod$LifeHistory$steepness

  # CAW Data
  data$WghtBins <- SimMod$Data$Weight_Bins
  data$WghtMids <- SimMod$Data$Weight_Mids

  CAW_DF <- SimMod$Data_CAW_DF %>% filter(Sim==sim) %>% select(Month_ind, Weight, Count)
  nBins <- length(unique(CAW_DF$Weight))
  nMonths <- length(unique(CAW_DF$Month_ind))

  CAW <- matrix(CAW_DF$Count, nrow=nBins, nMonths)
  data$CAW <- CAW

  data$CAW_ESS <- rep(CAW_Monthly_ESS, nMonths)

  # Effort Index
  Effort_DF <- SimMod$Data_TS_DF %>% filter(Sim==sim) %>% select(Effort)
  data$Effort <- Effort_DF$Effort
  data$Effort_SD <- rep(Effort_CV, nMonths)

  # CPUE Index
  CPUE_DF <- SimMod$Data_TS_DF %>% filter(Sim==sim) %>% select(CPUE)
  data$CPUE <- CPUE_DF$CPUE
  data$CPUE_SD <- rep(CPUE_CV, nMonths)

  # Priors and penalties
  data$F_meanprior <- 0

  data$Fit_Effort <- 1
  data$Fit_CPUE <- 1
  data$use_Frwpen <- 1
  data$use_R0rwpen <- 1
  data$use_Fmeanprior <- 0

  data$model <- 'SLAM'
  data$currentYr <- max(SimMod$OM_DF$Year)

  data
}
