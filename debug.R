
devtools::install_github('Blue-Matter/SLAM',
                         auth_token='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')



# Compare simulation and estimation
library(SLAM)
library(dplyr)
source("R/functions.r")
source('R/PopDynamics.R')
nts <- 120
q <- 0.2

WghtBins <- seq(0, 6, by=0.1)
WghtMids <- seq(0.05, by=0.1, length.out=length(WghtBins))

data <- SLAM::casestudydata
Pars <- list()
Pars$maxage <- 14
Pars$M <- 0.15
Pars$Weight_Age <- data$Weight_Age
Pars$Weight_Age_SD <- data$Weight_Age_SD
Pars$Mat_at_Age <- data$Mat_at_Age
Pars$PSM_at_Age <- data$PSM_at_Age
Pars$h <- 0.99 #data$h
Pars$WghtBins <-  WghtBins
Pars$WghtMids <-  WghtMids

Pars$rec_mu <- rep(1/12,12)
Pars$Rbar <- 1E5

rec_sd <- 0.001
rec_devs <- exp(rnorm(nts, -0.5*rec_sd^2, rec_sd))
Pars$rec_devs <- matrix(rec_devs, ncol=1)

# selectivity
Pars$sW50 <- 0.9
Pars$sW95 <- 1

# effort & fishing mortality
Pars$Effort <- data.frame(Effort=rep(1, nts), Sim=1, t=1:nts, Year=1, Month=1, Date=1)
Pars$q <- q

# Observation process
Pars$nyears <- 5 # number of years of data
Pars$sigmaI <- 0.01
Pars$sigmaE <- 0.01

Pars$CAW_nsamp <- 1000
Pars$CAW_ESS <- 5000

Pars$Rbar <- 1
SimPop <- Simulate(Pars)

data <- list()
nts <- length(SimPop$Biomass)
ndata_ts <- length(SimPop$Eff_ind)
data_ts <- (nts-ndata_ts+1):nts
data$Weight_Age <- Pars$Weight_Age
data$Weight_Age_SD <- Pars$Weight_Age_SD
data$M_at_Age <- SimPop$M_at_Age
data$Mat_at_Age <- Pars$Mat_at_Age
data$PSM_at_Age <- Pars$PSM_at_Age
data$h  <- Pars$h
data$Effort <- SimPop$Eff_ind
data$Effort_SD <- rep(0.01, ndata_ts)
data$CPUE <- SimPop$Index
data$CPUE_SD <- rep(0.01, ndata_ts)
data$WghtBins <- Pars$WghtBins
data$WghtMids <- Pars$WghtMids
data$CAW <- SimPop$CAW_samp
data$CAW_ESS <- rep(Pars$CAW_ESS, ndata_ts)
data$model <- 'SLAM'
data$Fit_Effort <- 1
data$Fit_CPUE <- 1
data$use_Frwpen <- 1
data$use_R0rwpen <- 1
data$use_Fmeanprior <- 0
data$F_meanprior <- 1



run <- Assess(data)
run$sdreport
run$chk

# profile over ls50
ls50vec <- seq(5.5, 6.5, by=0.01)
nll <- rep(NA, length(ls50vec))
for (i in seq_along(ls50vec)) {
  message(i, '/', length(ls50vec))
  parameters2 <- parameters

  map2 <- map
  map2$ls50 <- factor(NA)

  parameters2$ls50 <- log(ls50vec[i])
  map2$lsdelta <- factor(NA)
  parameters2$lsdelta <- log(0.5)

  obj <- TMB::MakeADFun(data=data, parameters2, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, map=map2, random=Random)
  # lower <- c(log(0.1), rep(-Inf, length(obj$par)-1))
  # upper <- c(log(2), rep(Inf, length(obj$par)-1))
  starts <- obj$par
  # opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr,
  #                                    control = control,
  #                                    upper=upper, lower=lower)),silent=TRUE)

  opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr,
                                     control = control)),silent=TRUE)

  if (class(opt)=='list') {
    nll[i] <- opt$objective
  }

}


plot(ls50vec, nll, type='b')

map2 <- map
# map2$lsdelta <- factor(NA)
# parameters2$lsdelta <- log(0.4)
obj <- TMB::MakeADFun(data=data, parameters2, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, map=map2, random=Random)
starts <- obj$par
opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr,
                                   control = control)),silent=TRUE)
rep <- obj$report()
rep$S50
rep$S95
sd <- TMB::sdreport(obj, obj$env$last.par.best)
sd
sd$gradient.fixed[1:3]


comparefits <- function(s50) {
  parameters2 <- parameters

  map2 <- map
  map2$ls50 <- factor(NA)
  parameters2$ls50 <- log(s50)

  upper <- c(log(12), c(3), rep(length(parameters2)-2))
  lower <- c(log(0.5), c(0.1), rep(length(parameters2)-2))
  obj <- TMB::MakeADFun(data=data, parameters2, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, map=map2, random=Random)

  rep <- obj$report(obj$env$last.par.best)
  out <- list()
  out$selA <- rep$selA
  out$CPUE <- rep$stpredCPUE
  out$CAW <- rep$predCAW
  out$Effort <-rep$StEffort
  out$nll <- rep$nll_joint
  out
}

run1 <- comparefits(4.7)
run2 <- comparefits(4.9)
run3 <- comparefits(5)
run4 <- comparefits(5.1)

S50 <- 4.9
Sdelta <- 0.1
a <- 0:14
selA1= 1 / (1 + exp(-log((19))*((a - 4.7)/Sdelta)));
selA2= 1 / (1 + exp(-log((19))*((a - 4.9)/Sdelta)));
cbind(selA1, selA2)
plot(selA1, type='l')
lines(selA2, col='blue')
run1$nll %>% sum()
run2$nll %>% sum()
run3$nll %>% sum()
run4$nll %>% sum()

plot(run1$selA, type='l')
lines(run2$selA, col='blue')
lines(run3$selA, col='red')
lines(run4$selA, col='green')


sdreport$cov[1:3,1:3]

rep$S50
rep$S95
plot(rep$selA, type='l')

plot(2:6, data$Weight_Age[3:7], type='b')
abline(h=data$WghtMids, lty=3)
abline(v=2:6, lty=3)


obj <- TMB::MakeADFun(data=data, parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, map=map, random=Random)

starts <- obj$par
lower <- c(log(0.1), rep(-Inf, length(obj$par)-1))
upper <- c(log(0.8*max(data$Weight_Age)), rep(Inf, length(obj$par)-1))
opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr,
                                   control = control)),silent=TRUE)

obj$report()$S50
obj$report()$S95
sd <- sdreport(obj ,obj$env$last.par.best)
plot(abs(sd$gradient.fixed[1,]))
plot(abs(sd$cov[1,]))

Pars$sW50
Pars$sW95
plot(SimPop$Sel_at_Age, type='l')
lines(rep$selA, col='blue')

plot(data$Weight_Age, SimPop$Sel_at_Age, type='l')
lines(data$Weight_Age, rep$selA, col='blue')

par(mfrow=c(4,5), bty='n', mar=rep(0,4), oma=rep(0,4))
for (i in 1:20) {
  plot(data$WghtMids, (data$CAW[,101:120]/apply(data$CAW[,101:120], 2, sum))[,i] ,
       type='l', xlab='',     ylab='')
  lines(data$WghtMids, (rep$predCAW[,41:60])[,i], col='blue')

}




options=list()
control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20)
map=list(log_sigmaF=factor(NA),
         log_sigmaR=factor(NA),
         log_sigmaR0=factor(NA))
Fit_Effort=1
Fit_CPUE=1
log_sigmaF=log(0.5)
log_sigmaR=log(0.9)
log_sigmaR0=log(0.6)

  # Starting parameters
  ls50 <- log(0.5)
  lsdelta <- log(0.15)
  logR0_m_est <- rep(log(1), 11)

  nts <- length(data$Effort)
  logF_m <- rep(log(0.01), nts)
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

  if (!is.null(map$log_sigmaR)) {
    Random <- NULL
  } else {
    Random <- 'logRec_Devs'
  }

  data$use_Fmeanprior <- 0
  data$F_meanprior <- 1

  if (!Fit_Effort) data$Fit_Effort <- 0
  if (!Fit_CPUE) data$Fit_CPUE <- 0

  # Group estimates from sequential time-steps with no data

  # Find time-steps with no data
  df <- data.frame(t=1:nts, CAW=apply(data$CAW, 2, sum)==0,
                   Effort=is.na(data$Effort),
                   CPUE=is.na(data$CPUE)) %>%
    filter(CAW==TRUE, Effort==TRUE, CPUE==TRUE)

  # if no data in first time-step, don't estimate F_minit
  df1 <- df %>% dplyr::filter(t==1)
  if (nrow(df1)>0) {
    if (all(df1[1,]))
      map$logF_minit <- factor(NA)
  }

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


  # parameters$log_sigmaR  <- log(0.05)

  # map$logR0_m_est <- factor(rep(NA, 11))
  # map$logRec_Devs <- factor(rep(NA, length(logRec_Devs)))
  # map$ls50 <- factor(NA)
  # map$lsdelta <- factor(NA)


  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, map=map, random=Random)

  starts <- obj$par
  lower <- c(log(0.1), log(0.01), rep(-Inf, length(obj$par)-2))
  upper <- c(log(0.8*max(data$Weight_Age)),
             log(2), rep(Inf, length(obj$par)-2))
  opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control
                                     )),silent=TRUE)


  rep <- obj$report(obj$env$last.par.best)
  sdreport <- TMB::sdreport(obj, obj$env$last.par.best)

  ind <- which.max(sdreport$gradient.fixed)
  print(ind)
  obj$par[1:ind]

  mean(rep$logRec_Devs)
  exp(rep$logRec_Devs)

  plot(rep$selA)
  lines(SimPop$Sel_at_Age)

  plot(data$Weight_Age, rep$selA)
  lines(data$Weight_Age, SimPop$Sel_at_Age)

  plot(SimPop$Rec_Pattern, ylim=c(0,0.5))
  lines(rep$R0_m)

  plot(SimPop$F_m, pch=16, type='l')
  lines(rep$F_m, col='blue')

  # N at age - compare
  t <- 10
  plot(SimPop$Number[,t], type='l')
  lines(rep$N_m[,t], col='blue')

  cbind(SimPop$Number[,t], rep$N_m[,t])
  a <- 2:15
  a2 <- a-1
  dn <- SimPop$Number[a, t]/SimPop$Number[a2, t-1]

  dn2 <- rep$N_m[a, t]/rep$N_m[a2, t-1]
  cbind(  -log(dn) - 0.15,-log(dn2) - 0.15)

  rep$F_m[t]

  SimPop$Number[1,]
  rep$N_m[1,]

  plot(SimPop$Eff_ind, type='l', ylim=c(0,1.5))
  lines(rep$StEffort, col='blue')

  plot(SimPop$Index, type='l', ylim=c(0,1.5))
  lines(rep$stpredCPUE, col='blue')

  plot(SimPop$Catch_Biomass[61:120]/mean(SimPop$Catch_Biomass[61:120]), type='l', ylim=c(0,1))
  lines(rep$predCB/mean(rep$predCB), col='blue')

  plot(data$WghtMids, data$CAW[,t+60]/sum(data$CAW[,t+60]), type="l")
  lines(data$WghtMids, rep$predCAW[,t])



