remotes::install_github('Blue-Matter/SLAM')

library(SLAM)
library(dplyr)
library(ggplot2)
library(ggthemes)

source('R/functions.R')


data <- SLAM::casestudydata

nyears <- 10
nts <- nyears * 12
nsim <- 100
set.seed(1001)

Effort_DF <- data.frame(t=1:nts,
                        Effort=1,
                        Sim=rep(1:nsim, each=nts),
                        Year=rep(1:nyears, each=12))
Effort_DF$M <- Effort_DF$t %% 12
Effort_DF$M[Effort_DF$M==0] <- 12
Effort_DF$Month <- month.abb[Effort_DF$M]
Effort_DF$Year <- Effort_DF$Year + 2011
Effort_DF$Date <- lubridate::my(paste(Effort_DF$Month, Effort_DF$Year, sep="-"))

Pars <- list()
Pars$Name <- 'Test'
Pars$maxage <- 14
Pars$M <- 0.15
Pars$Weight_Age <- data$Weight_Age
Pars$Weight_Age_SD <- data$Weight_Age_SD
Pars$Mat_at_Age <- data$Mat_at_Age
Pars$PSM_at_Age <- data$PSM_at_Age
Pars$h <- data$h
Pars$WghtBins <-  data$WghtBins
Pars$WghtMids <-  data$WghtMids

# recruitment pattern
mulist <- list(1, 6.5, 6.5, c(3,9))
sdlist <- list(100, 1, 3, c(1.5,1.5))
reclist <- list()
for (i in 1:length(mulist)) {
  reclist[[i]] <- GenMonthlyRec(mulist[[i]], sdlist[[i]])
}

Pars$rec_mu <- rep(1/12, 12)
Pars$Rbar <- 1E6

# recruitment deviations
rec_sd <- 0.000
rec_devs <- exp(rnorm(nts*nsim, -0.5*rec_sd^2, rec_sd))
Pars$rec_devs <- matrix(rec_devs, nrow=nts, ncol=nsim)

# selectivity
Pars$sA50 <- 5
Pars$sA95 <- 6

# effort & fishing mortality
Pars$Effort <- Effort_DF
Pars$q <- 0.2

Pars$Effort$Effort[Pars$Effort$Sim==1] <- 1

# Observation process
Pars$nyears <- 5 # number of years of data
Pars$sigmaI <- 0.02
Pars$sigmaE <- 0.02

Pars$CAW_nsamp <- 2000
Pars$CAW_ESS <- 2000
SimPop <- Simulate(Pars, sim=1)

# assess
sim_data <- list()
nts <- length(SimPop$Biomass)
ndata_ts <- length(SimPop$Eff_ind)
data_ts <- (nts-ndata_ts+1):nts
sim_data$Weight_Age <- Pars$Weight_Age
sim_data$Weight_Age_SD <- Pars$Weight_Age_SD
sim_data$M_at_Age <- SimPop$M_at_Age
sim_data$Mat_at_Age <- Pars$Mat_at_Age
sim_data$PSM_at_Age <- Pars$PSM_at_Age
sim_data$h  <- Pars$h
sim_data$Effort <- SimPop$Eff_ind
sim_data$Effort_SD <- rep(0.1, ndata_ts)
sim_data$CPUE <- SimPop$Index
sim_data$CPUE_SD <- rep(0.1, ndata_ts)
sim_data$WghtBins <- Pars$WghtBins
sim_data$WghtMids <- Pars$WghtMids
sim_data$CAW <- SimPop$CAW_samp
sim_data$CAW_ESS <- rep(Pars$CAW_ESS, ndata_ts)
sim_data$model <- 'SLAM'
sim_data$Fit_Effort <- 1
sim_data$Fit_CPUE <- 1
sim_data$use_Frwpen <- 1
sim_data$use_R0rwpen <- 1
sim_data$use_Fmeanprior <- 0
sim_data$F_meanprior <- 1

# Estimation
data <- sim_data
# Starting parameters
ls50 <- log(5)
lsdelta <- log(1)
logR0_m_est <- rep(log(1), 11)

log_sigmaF=log(0.01)
log_sigmaR=log(0.01)
log_sigmaR0=log(0.01)

nts <- length(data$Effort)
histF <- Pars$Effort$Effort[Pars$Effort$Sim==1]*Pars$q

ts <- (length(SimPop$Biomass)-nts+1):length(SimPop$Biomass)
logF_m <- log(histF[ts])
logF_minit <- log(0.01)
logRec_Devs <- rep(0, nts-1)
parameters <- list(ls50=ls50,
                   lsdelta=lsdelta,
                   logR0_m_est=logR0_m_est,
                   log_sigmaR0=log_sigmaR0,
                   logF_m=logF_m,
                   logF_minit=logF_minit,
                   logRec_Devs=logRec_Devs,
                   log_sigmaF=log_sigmaF,
                   log_sigmaR=log_sigmaR)

control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20)

map=list(
         lsdelta=factor(NA),
         log_sigmaF=factor(NA),
         log_sigmaR=factor(NA),
         log_sigmaR0=factor(NA),
         logF_m=factor(rep(NA, length(logF_m))),
         logR0_m_est=factor(rep(NA, 11)),
         logRec_Devs=factor(rep(NA,nts-1)))

Random <- NULL

data$use_Fmeanprior <- 0
data$F_meanprior <- 1

if (is.null(map[["logF_m"]])) {
  map$logF_m <- 1:nts
  for (i in 2:nts) {
    if(i %in% df$t) {
      map$logF_m[i] <- map$logF_m[i-1]
    }
  }
  # don't estimate F for time-steps before data
  if (sum(map$logF_m ==1)>1) {
    ind <- which(map$logF_m==1)
    map$logF_m[ind[1:(length(ind)-1)]] <- NA
  }

  map$logF_m <- factor(map$logF_m)
}

if (is.null(map[["logRec_Devs"]])) {
  map$logRec_Devs <- 1:(nts-1)
  for (i in 2:nts) {
    if(i %in% df$t) {
      map$logRec_Devs[i] <- map$logRec_Devs[i-1]
    }
  }
  map$logRec_Devs <- factor(map$logRec_Devs)
}

obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, map=map, random=Random)

starts <- obj$par
opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control)),silent=TRUE)

rep <- obj$report(obj$env$last.par.best)

#### TO DO ####
# get the initial year right

plot(SimPop$Number[,ts[1]]/Pars$Rbar, type='l')
lines(rep$N_m[,1], col='blue')


rep$F_minit






# plot model fits
y <- 1
plot(SimPop$CAW_samp[,y]/sum(SimPop$CAW_samp[,y]), type='l')
lines(rep$predCAW[,y], col='blue')

plot(SimPop$Index, type='l', ylim=c(0,2))
lines(rep$stpredCPUE, col='blue')



plot(simn[,1]/sum(simn[,1]))
lines(estn[,1]/sum(estn[,1]))
























# compare assessment and simulation

do_assess$rep$S50
do_assess$rep$S95






