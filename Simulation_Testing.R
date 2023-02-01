library(SLAM)
library(dplyr)
library(ggplot2)

LifeHistory <- Import_LifeHistory(dir='inst')
Data <- Import_Data(dir='inst')
Exploitation <- Import_Exploitation(dir='inst')


MyMod <- Simulate(LifeHistory, Exploitation, Data)



mod <- MyMod$OM_DF %>% filter(Sim==1)

sim <- 3
# Fit Assessment Model

# ---- data ----
data <- list()
data$Weight_Age <- MyMod$LifeHistory$Weight_Age_Mean
data$Weight_Age_SD <- MyMod$LifeHistory$Weight_Age_SD
data$Mat_at_Age <- MyMod$LifeHistory$Maturity_at_Age
data$M_at_Age <- MyMod$LifeHistory$M_at_Age
data$PSM_at_Age <- MyMod$LifeHistory$Post_Spawning_Mortality

data$WghtBins <- MyMod$Data$Weight_Bins
data$WghtMids <- MyMod$Data$Weight_Mids

CAW_DF <- MyMod$Data_CAW_DF %>% filter(Sim==sim) %>% select(Month_ind, Weight, Count)
nBins <- length(unique(CAW_DF$Weight))
nMonths <- length(unique(CAW_DF$Month_ind))

CAW <- matrix(CAW_DF$Count, nrow=nBins, nMonths)

data$CAW <- CAW

data$CAW_ESS <- rep(20, nMonths)

Effort_DF <- MyMod$Data_TS_DF %>% filter(Sim==sim) %>% select(Effort)
data$Effort <- Effort_DF$Effort
data$Effort_SD <- rep(0.1, nMonths)

CPUE_DF <- MyMod$Data_TS_DF %>% filter(Sim==sim) %>% select(CPUE)
data$CPUE <- CPUE_DF$CPUE
data$CPUE_SD <- rep(0.1, nMonths)

data$h <- MyMod$LifeHistory$steepness

data$F_meanprior <- 0

data$Fit_Effort <- 1
data$Fit_CPUE <- 1
data$use_Frwpen <- 1
data$use_R0rwpen <- 1
data$use_Fmeanprior <- 0

# ---- parameters ----
log_sigmaF=log(0.1)
log_sigmaR=log(0.3)
log_sigmaR0=log(0.4)

map=list(log_sigmaF=factor(NA),
         log_sigmaR=factor(NA),
         log_sigmaR0=factor(NA))

# Starting parameters
ls50 <- log(5)
lsdelta <- log(0.1)
logR0_m_est <- rep(log(1), 11)

n_age <- length(data$Weight_Age)
nts <- length(data$Effort)
logF_m <- rep(log(0.2), nts)
logF_minit <- rep(log(0.25), n_age)
logRec_Devs <- rep(0, nts-1)
parameters <- list(ls50=ls50,
                   lsdelta=lsdelta,
                   logF_minit=logF_minit,
                   logF_m=logF_m,
                   log_sigmaF=log_sigmaF,
                   logR0_m_est=logR0_m_est,
                   log_sigmaR0=log_sigmaR0,
                   logRec_Devs=logRec_Devs,
                   log_sigmaR=log_sigmaR)


if (!is.null(map$log_sigmaR)) {
  Random <- NULL
} else {
  Random <- 'logRec_Devs'
}


obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, map=map, random=Random)


lapply(parameters, class)
lapply(data, class)



Assess <- function() {

}



plot(mod$Month_ind, mod$Effort, type='l')

plot(mod$Month_ind, mod$Catch, type='l')
plot(mod$Month_ind, mod$SB_fished/mod$SB_unfished, type='l', ylim=c(0,1))

lines(mod$Month_ind, mod$SPR, type='l', col='blue')

plot(mod$SB_fished/mod$SB_unfished, mod$Recruits, type='l', ylim=c(0,max(mod$Recruits)))
plot( mod$Recruits, type='l', ylim=c(0,max(mod$Recruits)))
