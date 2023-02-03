#' Iniitialize Parameters for SLAM Estimation
#'
#' @param data A data object
#' @param as50 initial value for age at 50% selection
#' @param as95 initial value for age at 95% selection
#' @param Feq_init
#' @param sigmaR
#' @param log_sigmaF
#' @param log_sigmaR0
#' @param log_sigmaR
#'
#' @return
#' @export
#'
Initialize_Parameters <- function(data,
                                  as50=4, as95=6,
                                  Feq_init=0.05,
                                  sigmaR=0.3,
                                  log_sigmaF=log(0.05),
                                  log_sigmaR0=log(0.1),
                                  log_sigmaR=log(0.3)) {
  parameters <- list()

  parameters$ls50 <- log(as50)
  parameters$lsdelta <- log(as95-as50)

  n_age <- length(data$Weight_Age)
  parameters$logF_minit <- log(0.05)

  n_ts <- length(data$Effort)

  parameters$logF_m <- rep(log(mean(data$M_at_Age)), n_ts)
  parameters$log_sigmaF <- log_sigmaF # standard deviation for random walk penalty for F
  parameters$logR0_m_est <- rep(1/12, 12)
  parameters$log_sigmaR0 <- log_sigmaR0 # sd for random walk penalty for monthly recruitment
  parameters$logRec_Devs <- rep(log(1), n_ts)
  parameters$log_sigmaR  <- log(sigmaR) # monthly rec dev sd (usually fixed)

  parameters
}


#' Title
#'
#' @param SimMod
#' @param log_sigmaF
#' @param log_sigmaR0
#' @param log_sigmaR
#'
#' @return
#' @export
Initialize_Parameters <- function(SimMod,
                                  log_sigmaF=log(0.05),
                                  log_sigmaR0=log(0.1),
                                  log_sigmaR=log(0.3)) {
  parameters <- list()

  parameters$ls50 <- log(SimMod$Exploitation$SA50)
  parameters$lsdelta <- log(SimMod$Exploitation$SA95-SimMod$Exploitation$SA50)

  n_age <- SimMod$LifeHistory$maxage+1
  parameters$logF_minit <- log(0.05)


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

#' Title
#'
#' @param SimMod
#' @param sim
#' @param CAW_Monthly_ESS
#' @param Effort_CV
#' @param CPUE_CV
#' @param Fit_Effort
#' @param Fit_CPUE
#' @param Fit_CAW
#' @param use_Frwpen
#' @param use_R0rwpen
#' @param use_Fmeanprior
#'
#' @return
#' @export
Construct_Data_OM <- function(SimMod, sim=1,
                           CAW_Monthly_ESS=100,
                           Effort_CV=0.2,
                           CPUE_CV=0.2,
                           Fit_Effort=1,
                           Fit_CPUE=1,
                           Fit_CAW=1,
                           use_Frwpen=0,
                           use_R0rwpen=0,
                           use_Fmeanprior=0) {

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

  data$Fit_Effort <- Fit_Effort
  data$Fit_CPUE <- Fit_CPUE
  data$Fit_CAW <- Fit_CAW
  data$use_Frwpen <- use_Frwpen
  data$use_R0rwpen <- use_R0rwpen
  data$use_Fmeanprior <- use_Fmeanprior

  data$model <- 'SLAM'
  data$currentYr <- max(SimMod$OM_DF$Year)

  data
}
