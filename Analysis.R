

# Application to Case Study ----









map=list(log_sigmaF=factor(NA),
         log_sigmaR=factor(NA),
         log_sigmaR0=factor(NA),
         logR0_m_est=factor(rep(NA,11)))

data <- SLAM::casestudydata
data$use_R0rwpen <- 1

# data$CAW_ESS <- rep(500, length(data$CAW_ESS))
Base_Case <- Assess(data)


Base_Case$rep$S50
Base_Case$rep$S95

calcMonthlyMean <- function(dat, offset=0) {

  nts <- length(dat)
  df <- data.frame(t=1:nts, n=dat)
  df$n[df$n==0] <- NA
  df$M <- (df$t+ offset) %% 12
  df$M[df$M ==0] <- 12
  df$Month <- month.abb[df$M ]
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)
  df %>% group_by(Month) %>%
    summarize(n=mean(n, na.rm=TRUE))
}

tt <- calcMonthlyMean(Base_Case$rep$N_m[1,])
small_catch <- calcMonthlyMean(data$CAW[3,], 5)

lines(tt$Month, tt$n/sum(tt$n), col='green')
nn <- Base_Case$rep$R0_m+tt$n/sum(tt$n)
nn <- nn/sum(nn)
plot(tt$Month, tt$n/sum(tt$n), ylim=c(0,0.5), type='l')
lines(tt$Month, nn, col='red')
lines(small_catch$Month, small_catch$n/sum(small_catch$n), col='blue')


Base_Case$rep$N_m[1,]


# Plot Fit to CAW
plotCAW <- function(assess, first_yr=2017) {
  wght_mids <- assess$obj$env$data$WghtMids
  caw_data <- assess$obj$env$data$CAW
  d <- dim(caw_data)
  nbins <- d[1]
  nts <- d[2]
  df <- data.frame(t=rep(1:nts, each=nbins),
                   Bin=wght_mids, each=nts,
                   Obs=as.numeric(caw_data),
                   Pred=as.numeric(assess$rep$predCAW))
  df$M <- df$t %% 12
  df$M[df$M ==0] <- 12
  df$Month <- month.abb[df$M]
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)

  df$Year <- rep(Year, each=12*nbins)[1:nrow(df)]
  df <- df %>% group_by(t) %>%
    mutate(Obs_p=Obs/sum(Obs), Pred_p=Pred/sum(Pred))

  ggplot(df, aes(x=Bin)) +
    facet_grid(Month~Year) +
    geom_bar(stat='identity', aes(y=Obs_p)) +
    geom_line(aes(y=Pred_p))

}






# Base Case and Sensitivity Tests
data <- SLAM::casestudydata
data_HigherM <- data_LowerH <- data_HigherH <- data
data_HigherM$M_at_Age <- rep(0.2, length(data_HigherM$M_at_Age))
data_HigherM_LowerH <- data_HigherM_HigherH <- data_HigherM

data_LowerH$h <- 0.6
data_HigherH$h <- 0.9

data_HigherM_LowerH$h <- 0.6
data_HigherM_HigherH$h <- 0.9

data$use_R0rwpen <- 0
Base_Case <- Assess(data)

Base_Case_optY <- Optimize(data, Rec_Pattern=Base_Case$rep$R0_m, selA=Base_Case$rep$selA, opt_type = 0)
Base_Case_optU <- Optimize(data, Base_Case$rep$R0_m, Base_Case$rep$selA)

par(mfrow=c(2,2))
plot(Base_Case$rep$R0_m, type='l')
plot(Base_Case_optY$F_m, type='l')
lines(Base_Case$rep$meanF, col='blue')

plot(Base_Case_optU$F_m, type='l')
lines(Base_Case$rep$meanF, col='blue')

plot(Base_Case$rep$meanSPR, type='l', ylim=c(0,1))
lines(Base_Case_optY$SPR, col='blue')
lines(Base_Case_optU$SPR, col='red')

# UNDERSTAND WHY OPTIMIZING FOR YIELD RESULTS IN THIS PATTERN
# How can SPR be 0???

# Check SPR corresponds with equilibrium SB/SB0

# IS THE ESTIMATED RECRUITMENT PATTERN RELIABLE?
# - try with and without random walk penalty

# Compare with constant recruitment?





HigherM <- Assess(data_HigherM)
LowerH <- Assess(data_LowerH)
HigherH <- Assess(data_HigherH)
HigherM_LowerH <- Assess(data_HigherM_LowerH)
HigherM_HigherH <- Assess(data_HigherM_HigherH)


# Optimize for maximum yield and utility
Optimize <- function(data, Rec_Pattern, selA, opt_type=1, utilpow=0.4) {

  data <- list(model='optF',
               rec_pattern=Rec_Pattern,
               opt_type=opt_type,
               utilpow=utilpow,
               h=data$h,
               Wght_Age=data$Weight_Age,
               Mat_at_Age=data$Mat_at_Age,
               M_at_Age=data$M_at_Age,
               PSM_at_Age=data$PSM_at_Age,
               selA=selA)

  parameters <- list(logF_m=rep(log(0.01),12))
  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE)
  starts <- obj$par
  opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr))
  obj$report(obj$env$last.par.best)
}

Base_Case_optY <- Optimize(data, Rec_Pattern=Base_Case$rep$R0_m, selA=Base_Case$rep$selA, opt_type = 0)

par(mfrow=c(2,2))
plot(Base_Case$rep$R0_m, type='l')
plot(Base_Case_optY$F_m, type='l')
plot(Base_Case_optY$predCB, type='l')

# UNDERSTAND WHY OPTIMIZING FOR YIELD

Base_Case$rep$selA


Base_Case_optU <- Optimize(data, Base_Case$rep$R0_m, Base_Case$rep$selA)

HigherM_optY <- Optimize(data, HigherM$rep$R0_m, HigherM$rep$selA, opt_type = 0)
HigherM_optU <- Optimize(data_HigherM, HigherM$rep$R0_m, HigherM$rep$selA)


plot(Base_Case_optY$F_m, type='l')
lines(HigherM_optY$F_m, col='blue')

plot(Base_Case_optU$F_m, type='l')
lines(HigherM_optU$F_m, col='blue')

plot(Base_Case$rep$R0_m, type='l')
lines(HigherM$rep$R0_m, col='blue')
lines(LowerH$rep$R0_m, col='red')
lines(HigherH$rep$R0_m, col='green')
lines(HigherM_LowerH$rep$R0_m, col='orange')
lines(HigherM_HigherH$rep$R0_m, col='grey')

nts <- length(Base_Case$rep$F_m)
df <- data.frame(t=1:nts, F=Base_Case$rep$F_m)
df$M <- df$t %% 12
df$M[df$M ==0] <- 12
df$Month <- month.abb[df$M]
df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)

df2 <- df %>% group_by(Month) %>% summarize(Mean=mean(F))

maxY <- max(c(Base_Case_optY$F_m, Base_Case_optU$F_m, df2$Mean))
plot(df2$Mean, type='l', ylim=c(0, maxY))
lines(Base_Case_optY$F_m,  col='blue')
lines(Base_Case_optU$F_m, col='red')

# Check

# plot fits

# Plot model fits
par(mfrow=c(5,6))

plot(data$CPUE, pch=16, type='l', ylim=c(0, 1.2))
lines(rep$stpredCPUE, col='blue')

plot(data$Effort, pch=16, type='l')
lines(rep$StEffort, col='blue')

plot(rep$R0_m,type='l', ylim=c(0, max(rep$R0_m)))
plot(rep$F_m,  type='l')
plot(rep$SPR,  type='l', ylim=c(0,1))

for (i in 1:24) {
  caw <-  data$CAW[,i]
  plot(data$WghtMids, caw, type='l')
  lines(data$WghtMids, sum(caw)*Base_Case$rep$predCAW[,i], col='blue')
}

plot(rep$selA, type='l')

# ---- Assessment Figures -----


plot_CPUE <- function(data, results) {

}


# Things to follow up on:
# - why super high F in time-step 2?
# - doesn't work well without random walk penalties? add description to report
# - why does F go to zero around ts 40?
# - modify optF to take selA
# - make HTML report
# - report Effort and CPUE on same scale as Data

par(mfrow=c(2,2))
plot(rep$StEffort,type='l', pch=16)
plot(rep$F_m,type='l', pch=16, ylim=c(0, 1))
plot(rep$SPR,type='l', pch=16, ylim=c(0, 1))








lower <- rep(-Inf, length(starts))
upper <- rep(Inf, length(starts))
lower[1] <- log(0.05)
upper[1] <- log(0.8)
lower[2] <- log(0.1)

opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control,
                               lower=lower, upper=upper))

opt$objective

rep <- obj$report()
sdreport <- TMB::sdreport(obj)

# - test model

# - add BH SRR

# - add option to fit to length comps















# Formalize processing of octopus data

# get model to reliably fit to empirical data
# generate assessment figures
# write-up

# TO DO TODAY:
# - modify model to fit to weight composition DONE
# - process Indo data and write up
# - modify model to fit to relative trend in effort in each year

# - fit to Indo data
# - re-write report following new structure



Octopus_Data <- Process_Data(IndoData)


n_ages <- length(mean_length)
first_yr <- min(Octopus_Data$Catch$Year)
last_yr <- max(Octopus_Data$Catch$Year)
years <- first_yr:last_yr
n_year <- length(years)
n_months <- 12 * n_year

# CPUE
CPUE <- CPUE_SD <- rep(NA, n_months)
t <- 0
for (y in years) {
  for (m in 1:12) {
    t <- t+1
    tt <- Octopus_Data$CPUE %>% filter(Year==y, Month==m)
    if (nrow(tt)>0) {
      CPUE[t] <- tt$Mean
      CPUE_SD[t] <- tt$SD
    }
  }
}

# Effort
Effort <- Effort_SD <- rep(NA, n_months)
t <- 0
for (y in years) {
  for (m in 1:12) {
    t <- t+1
    tt <- Octopus_Data$Effort %>% filter(Year==y, Month==m)
    if (nrow(tt)>0) {
      Effort[t] <- tt$Mean
      Effort_SD[t] <- tt$SD
    }
  }
}

# CAL
nbins <- length(Octopus_Data$Length_Mids)
CAL <- matrix(0, nrow=nbins, ncol=n_months)
t <- 0
for (y in years) {
  for (m in 1:12) {
    t <- t+1
    tt <- Octopus_Data$Length_Comp %>% filter(Year==y, Month==m, Sex=='Betina')
    if(nrow(tt>0))
      CAL[,t] <- tt$Count
  }
}

CAL_ESS <- apply(CAL, 2, sum)
CAL_ESS[CAL_ESS>200] <- 200

# Last available month
tt <- Octopus_Data$Effort %>% filter(Year==last_yr)
last_month <- max(tt$Month)

first_yr <- min(Octopus_Data$Catch$Year)
last_yr <- max(Octopus_Data$Catch$Year)
years <- first_yr:last_yr
n_year <- length(years)-1
n_months <- (12 * n_year) + last_month

Effort <- Effort[1:n_months]
Effort_SD <- Effort_SD[1:n_months]
CPUE <- CPUE[1:n_months]
CPUE_SD <- CPUE_SD[1:n_months]
CAL <- CAL[,1:n_months]

# # impute missing effort and cpue
first <- FALSE
for (t in 1:length(Effort)) {
  chk <- is.na(Effort[t])
  if (!chk) first <- TRUE
  if (first) {
    if (chk) {
      Effort[t] <- Effort[t-1]
      Effort_SD[t] <- Effort_SD[t-1]
      CPUE[t] <- CPUE[t-1]
      CPUE_SD[t] <- CPUE_SD[t-1]
    }
  }

}



Data <- list()
Data$Len_Age <- mean_length
Data$SD_Len_at_Age <- 0.1 * mean_length
Data$Wt_at_Age <- weight2
Data$Mat_at_Age <- c(rep(0,11), 0.9, rep(1, 3))
Data$M_at_Age <- rep(0.1, n_ages)
Data$phi_at_Age <- Data$Mat_at_Age
Data$Len_Bins <- Octopus_Data$Length_Bins
Data$Len_Mids <- Octopus_Data$Length_Mids
Data$CAL <- CAL
Data$CAL_ESS <- CAL_ESS
Data$Effort <- Effort/mean(Effort, na.rm = TRUE)
Data$Effort_SD <- Effort_SD
Data$CPUE <- CPUE/mean(CPUE, na.rm = TRUE)
Data$CPUE_SD <- CPUE_SD
Data$F_meanprior <- c(0.1,0.1)
Data$Fit_Effort <- 1
Data$Fit_CPUE <- 1
Data$use_Fmeanprior <- 0
Data$use_Frwpen <- 0
Data$use_R0rwpen <- 0
Data$Vmaxlen <- 1
Data$maxL=Data$Len_Age[length(Data$Len_Age)]




control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20)
map=list(log_sigmaF=factor(NA),
         log_sigmaR=factor(NA),
         log_sigmaR0=factor(NA))

# Assess <- function(Data, options=list(),
#                    control=list(eval.max=1E4, iter.max=1E4),
#                    map=list(log_sigmaF=factor(NA),
#                             log_sigmaR=factor(NA),
#                             log_sigmaR0=factor(NA))) {


  # --- Checks on data ---
  #TODO
  data <- list(model='SLAM',
               Len_Age=Data$Len_Age,
               SD_Len_Age=Data$SD_Len_at_Age,
               Wght_Age=Data$Wt_at_Age,
               Mat_at_Age=Data$Mat_at_Age,
               M_at_Age=Data$M_at_Age,
               PSM_at_Age=Data$phi_at_Age,
               LenBins=Data$Len_Bins,
               LenMids=Data$Len_Mids,
               CAL=Data$CAL,
               CAL_ESS=Data$CAL_ESS,
               Effort=Data$Effort,
               Effort_SD=Data$Effort_SD,
               CPUE=Data$CPUE,
               CPUE_SD=Data$CPUE_SD,
               F_meanprior=Data$F_meanprior,
               Fit_Effort=Data$Fit_Effort,
               Fit_CPUE=Data$Fit_CPUE,
               use_Fmeanprior=Data$use_Fmeanprior,
               use_Frwpen=Data$use_Frwpen,
               use_R0rwpen=Data$use_R0rwpen,
               Vmaxlen=1,
               maxL=Data$Len_Age[length(Data$Len_Age)]
  )

  # starting values for parameters
  # selectivity - 2 parameters
  Lencumulative <- apply(data$CAL, 2, cumsum)
  ind <- apply(Lencumulative, 2, sum)
  ind2 <- which(ind>0)
  Lencumulative <- Lencumulative[,ind2]/matrix(Lencumulative[nrow(Lencumulative),ind2], nrow=nrow(Lencumulative), ncol=length(ind2),byrow=TRUE)

  Lencumulative <- apply(Lencumulative, 1, mean, na.rm=TRUE)
  maxL <- max(data$Len_Age)
  SL50_start <- data$LenMids[min(which(Lencumulative>0.3))]
  S95_start <- data$LenMids[min(which(Lencumulative>0.5))]
  lrelSL50 <- log(SL50_start/maxL)

  lSLdelta <- log(S95_start/SL50_start)

  # L5_start <- data$LenMids[min(which(Lencumulative>0.05))]
  # LFS_start <- data$LenMids[min(which(Lencumulative>0.5))]
  # t_sl5 <- logit(L5_start/maxL)
  # t_slfint <- logit((LFS_start-L5_start)/maxL)

  nts <- dim(Data$CAL)[2]
  logR0_m_est <- rep(log(1), 11)
  log_sigmaR0 <- log(0.6)
  logF_m <- rep(log(0.1), nts)
  logF_minit <- log(0.1)
  logRec_Devs <- rep(0, nts-1)
  log_sigmaF <- log(0.5)
  log_sigmaR <- log(0.6)
  parameters <- list(lrelSL50=lrelSL50,
                     lSLdelta=lSLdelta,
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

  data$use_Frwpen <- 1
  data$use_R0rwpen <- 1
  data$Fit_CPUE <- 1
  data$Fit_Effort <- 1

  control=list(eval.max=1E4, iter.max=1E4, abs.tol=1E-20, rel.tol=1E-6)


  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, map=map, random=Random)

  starts <- obj$par

  lower <- rep(-Inf, length(starts))
  upper <- rep(Inf, length(starts))
  lower[1] <- log(0.05)
  upper[1] <- log(0.8)
  lower[2] <- log(0.1)

  opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control,
                                 lower=lower, upper=upper))

  opt$objective

  rep <- obj$report()
  sdreport <- TMB::sdreport(obj)


  # Plot model fits
  par(mfrow=c(5,6))

  plot(Data$CPUE, pch=16, type='l', ylim=c(0, 1.2))
  lines(rep$stpredCPUE, col='blue')

  plot(Data$Effort, pch=16, type='l')
  lines(rep$StEffort, col='blue')

  plot(rep$R0_m,type='l')
  plot(rep$F_m,  type='l')
  plot(rep$SPR,  type='l', ylim=c(0,1))

  for (i in 1:24) {
    cal <-  Data$CAL[,i]
    plot(Data$Len_Mids,cal, type='l')
    lines(Data$Len_Mids, sum(cal)*rep$predCAL[,i], col='blue')
  }

  plot(Data$Len_Mids, rep$selL, type='l')







  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, map=map, random=Random)

  opt2 <-TMBhelper::fit_tmb(obj, obj$fn, obj$gr,
                           control=list(trace=0, eval.max=20000, iter.max=20000))

  opt2$objective




  rep <- obj$report()


  rep$nll_joint
  rep$CAL
  rep$SL50
  rep$SL95
  plot(data$LenMids, rep$selL)
  plot(data$LenMids, rep$predCAL[,3], type='l')


  sdreport <- TMB::sdreport(obj)
  list(rep=rep, sdreport=sdreport, opt=opt)
}


rep$nll_joint
rep$CALnll

rep$predCAL[,7]



plot(rep$F_m, rep$SPR, ylim=c(0,1))






doassess <- Assess(Data)



rep <- doassess$rep

# Plot model fits
par(mfrow=c(5,6))

plot(Data$CPUE, pch=16, type='l', ylim=c(0, 1.2))


plot(Data$Effort, pch=16, type='l')


plot(rep$R0_m,type='l')
plot(rep$F_m,  type='l')

for (i in 1:24) {
  cal <-  Data$CAL[,i]
  plot(Pop$Len_Mids, cal, type='l')
  lines(Pop$Len_Mids, sum(cal)*rep$predCAL[,i], col='blue')
}

plot(Pop$Sel_at_Length)
lines(rep$selL)





SimData <- function(Pars, CAL_ESS=200, Effort_SD=0.2, CPUE_SD=0.2,
                    Fmeanprior=NULL,
                    fitCPUE=1, fitEffort=1,
                    use_Frwpen=0,
                    use_R0rwpen=0,
                    use_years=5)


Pars <- Load_Scenario(1)
Pars$rec_sd <- 1
Pop <- Simulate(Pars)
Sim <- SimData(Pop, use_years = 5)  # simulate fishery
Data <- Sim$Data
Hist <- Sim$Pop

Sim$Data$CAL %>% dim()

doassess <- Assess(Data)
rep <- doassess$rep

assess_freq=12

Project <- function(Hist, MPs=NA, use_years=5, pyears=10) {
  Pars <- Hist$Pars
  dd <- dim(Hist$Number)
  nages <- dd[1]
  nyears <- dd[2]
  Hist$Number[,nyears]

  nts <- pyears * 12

  nMPs <- length(MPs)
  # Status quo future effort
  recentE <- Hist$Effort[(nts-11):nts] # effort in last 12 months
  Effort <- recentE * exp(rnorm(nts, -0.5*Pars$Ecv^2, Pars$Ecv))
  projEffort <- matrix(Effort, nrow=length(Effort), ncol=nMPs)

  for (m in 1:nMPs) {
    # Run initial assessment
    Sim <- SimData(Hist, use_years = use_years)
    doassess <- Assess(Sim$Data)

    # apply MP
    fn <- get(MPs[m])
    eff_adj <- fn(doassess, Sim$Data)

    projEffort[,m] <- recentE * eff_adj

  }



  plot(Effort[1:12], type='l')
  lines(recentE * (opt$F_m/recentF), col='blue')



  # Project forward



  for (t in 1:nts) { # projections
    m_ind <- t %%12 # calendar month
    if (m_ind ==0) m_ind <- 12

    if (m_ind %in% manage_m) {
      # assess and management recommendation


    }

    F <- Effort[t] * Hist$q
    F_at_age <-  F * Hist$Sel_at_Age
    Z_at_age <- Hist$M_at_Age + F_at_age






  }
}




# Management Procedures


MP_1 <- function(assess, Data, assumed_h=0.6, utilpow=0.3) {
  opt <- Optimize(Data, assess$rep$R0_m, assess$rep$selA, utilpow=utilpow, assumed_h = assumed_h)

  nts <- length(assess$rep$F_m)
  recentF <- assess$rep$F_m[(nts-11):nts]

  # effort adjustment for target F
  opt$F_m/recentF
}




# ---- Calculate Optimal Fishing Pattern ----


# Closed-Loop



# Plot model fits
par(mfrow=c(5,6))

plot(Data$CPUE, pch=16, type='l', ylim=c(0, 1.2))
lines(rep$stpredCPUE, col='blue')

plot(Data$Effort, pch=16, type='l')
lines(rep$StEffort, col='blue')

plot(Pop$Rec_Pattern[1:12],type='l', ylim=c(0, max(Pop$Rec_Pattern[1:12])))
lines(rep$R0_m, col='blue')

R <- matrix(rep$N_m[1,], nrow=12, ncol=5)
r <- apply(R, 1, mean)
r <- r/sum(r)
rep$Rec_Pattern <- r
plot(Pop$Rec_Pattern[1:12],type='l', ylim=c(0, max(Pop$Rec_Pattern[1:12])))
lines(r, col='blue')

plot(Pop$F_m[61:120], ylim=c(0, max(c(Pop$F_m, rep$F_m))), type='l')
lines(rep$F_m, col='blue')

for (i in 1:24) {
  cal <-  Data$CAL[,i]
  plot(Pop$Len_Mids, cal, type='l')
  lines(Pop$Len_Mids, sum(cal)*rep$predCAL[,i], col='blue')
}

plot(Pop$Sel_at_Length)
lines(rep$selL)



plot(Pop$SPR[241:360], type='l', ylim=c(0,1))
lines(rep$SPR, col='blue')

# - get model to match with selectivity fixed
# - fix the selectivity

rep$sigmaF
rep$sigmaR


rep$nll_joint


















# checks:
rep$selA/Pop$Sel_at_Age
rep$selL/Pop$Sel_at_Length

rep$ALK-Pop$ALK
rep$ALK_C-Pop$ALK_C

simN <- Pop$Number[,169:180]
estN <- rep$N_m

par(mfrow=c(2,6))
for (m in 1:12) {
  plot(0:15, simN[,m], type='l', main=month.abb[m])
  lines(0:15, estN[,m], col='blue')
}

Pop$Effort[169:180]/rep$StEffort
Pop$Catch_Biomass[169:180]/rep$predCB

plot(data$CPUE)
lines(rep$stpredCPUE)

tt <- rep$predCB/rep$StEffort
tt/mean(tt)



s










Pars <- Load_Scenario(1)

Pars$F_mu <- 0.3
Pars$F_sd <- 0.01
Pars$sigmaR <- 0.01
Pars$h <- 0.999
Pars$rec_mu <- 8
Pars$rec_sd <- 1
Pars$sigmaE <- 0.01
Pars$sigmaI <- 0.01
Pars$CAL_ESS <- 1000
Pars$Rbar <- 1000

# Make selectivity same dome-shaped
# Map the Vmaxlen parameter


# Simulate 1 year of fished data
# Set estimation model with exact parameters
# Confirm that simulation and estimation models are identical

Sim <- SimData(Pars, use_years = 1)
Data <- Sim$Data


Data$CPUE_SD <- Data$Effort_SD <- rep(0.01, length(Data$Effort_SD))
Data$log_sigmaR <- log(Pars$sigmaR)

Pop <- Sim$Pop
rep <- Assess(Data)
rep <- Assess(Data, map=list(log_sigmaR0=as.factor(NA)))

rep$sigmaR0

# Plot model fits
par(mfrow=c(5,6))

plot(Data$CPUE, pch=16, type='l', ylim=c(0, 1.2))
lines(rep$stpredCPUE, col='blue')

plot(Data$Effort, pch=16, type='l')
lines(rep$StEffort, col='blue')

plot(Pop$Rec_Pattern[1:12],type='l', ylim=c(0, max(Pop$Rec_Pattern[1:12])))
lines(rep$R0_m, col='blue')

# R <- matrix(rep$N_m[1,], nrow=12, ncol=5)
# r <- apply(R, 1, mean)
# r <- r/sum(r)
# rep$Rec_Pattern <- r
# plot(Pop$Rec_Pattern[1:12],type='l', ylim=c(0, max(Pop$Rec_Pattern[1:12])))
# lines(r, col='blue')

plot(Pop$F_m[61:120], ylim=c(0, max(c(Pop$F_m, rep$F_m))), type='l')
lines(rep$F_m, col='blue')

for (i in 1:24) {
  cal <-  Data$CAL[,i]
  plot(Pop$Len_Mids, cal, type='l')
  lines(Pop$Len_Mids, sum(cal)*rep$predCAL[,i], col='blue')
}





# Calculate Assessment Performance




opt_rep <- Optimize(Data, rep$Rec_Pattern, rep$selA)
opt_rep2 <- Optimize(Data, Pop$Rec_Pattern, Pop$Sel_at_Age, assumed_h = Pars$h)
plot(opt_rep$F_m, type='l', ylim=c(0, max(c(opt_rep$F_m, opt_rep$F_m2))))
lines(opt_rep2$F_m, col='blue')

plot(opt_rep$SPR, type='l', ylim=c(0,1))
lines(opt_rep2$SPR, col='blue')






lower <- rep(-Inf, 12)
upper <- rep(log(2), 12)
opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr))
# lower=lower, upper=upper))







 # Attempt with using phases

data <- list(model='SLAM',
             Len_Age=Data$Len_Age,
             SD_Len_Age=Data$SD_Len_at_Age,
             Wght_Age=Data$Wt_at_Age,
             Mat_at_Age=Data$Mat_at_Age,
             M_at_Age=Data$M_at_Age,
             PSM_at_Age=Data$phi_at_Age,
             LenBins=Data$Len_Bins,
             LenMids=Data$Len_Mids,
             CAL=Data$CAL,
             CAL_ESS=Data$CAL_ESS,
             Effort=Data$Effort,
             Effort_SD=Data$Effort_SD,
             CPUE=Data$CPUE,
             CPUE_SD=Data$CPUE_SD,
             sigmaRprior=Data$sigmaRprior,
             F_meanprior=Data$F_meanprior,
             Fit_Effort=Data$Fit_Effort,
             Fit_CPUE=Data$Fit_CPUE,
             use_sigmaRprior=Data$use_sigmaRprior,
             use_Fmeanprior=Data$use_Fmeanprior,
             use_Frwpen=Data$use_Frwpen,
             use_R0rwpen=Data$use_R0rwpen,
             log_sigmaR0=Data$log_sigmaR0,
             log_sigmaF=Data$log_sigmaF
)

# starting parameters
# selectivity - 2 parameters
modal_lengths <- data$LenMids[apply(data$CAL, 2, which.max)]
sl50_start <- mean(modal_lengths)
log_sl50 <- log(sl50_start)

modal_lengths <- data$LenMids[apply(data$CAL, 2, which.max)+1]
sl95_start <- mean(modal_lengths)
log_sldelta <- log(sl95_start-sl50_start)

# monthly R0 - assume constant
logR0_m_est <- rep(log(1), 11)

nts <- dim(Data$CAL)[2]
nages <- length(Data$Len_Age)
logF_m <- rep(log(0.1), nts)
logF_minit <- log(0.1)
logRec_Devs <- rep(0, nts)
logsigmaR <- log(0.3)

parameters <- list(log_sl50=log_sl50,
                   log_sldelta=log_sldelta,
                   logR0_m_est=logR0_m_est,
                   logF_m=logF_m,
                   logF_minit=logF_minit,
                   logRec_Devs=logRec_Devs)

random <- 'logRec_Devs'

phases <- list(
  log_sl50=2,
  log_sldelta=2,
  logR0_m_est=1,
  logsigmaR=1,
  logF_m=1,
  logF_minit=1,
  logRec_Devs=3
)

# copied and modified from https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R
fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}

#loop over phases
for (phase_cur in 1:max(unlist(phases))) {

  #work out the map for this phase
  # if phases for parameters is less than the current phase
  # then map will contain a factor filled with NAs
  map_use <- list()
  j <- 0
  for (i in 1:length(parameters)) {
    if (phases[[i]]>phase_cur) {
      j <- j+1
      map_use[[j]] <- fill_vals(parameters[[i]],NA)
      names(map_use)[j] <- names(parameters)[i]
    }
  }

  #remove the random effects if they are not estimated
  random_use <- random[!random%in%names(map_use)]

  # initialize the parameters at values in previous phase
  params_use <- parameters
  if (phase_cur>1) params_use <- obj$env$parList(opt$par)

  # Fit the model
  obj <- TMB::MakeADFun(data,params_use,random=random_use,DLL='SLAM_TMBExports',map=map_use,
                        silent=TRUE)
  # TMB::newtonOption(obj,smartsearch=FALSE)
  control <- list(eval.max=1E4, iter.max=2000, trace=0)
  opt <- nlminb(obj$par,obj$fn,obj$gr, control=control)
  sdrep <- TMB::sdreport(obj)
  rep <- obj$report()
  #close phase loop

}

data$sigmaRprior
rep$sigmaR


# Plot model fits
par(mfrow=c(4,5))

plot(Data$CPUE, pch=16, type='l', ylim=c(0, 1.2))
lines(rep$stpredCPUE, col='blue')

plot(Data$Effort, pch=16, type='l')
lines(rep$StEffort, col='blue')

plot(Pop$Rec_Pattern[1:12],type='l', ylim=c(0, max(Pop$Rec_Pattern[1:12])))
lines(rep$R0_m, col='blue')

plot(Pop$F_m[61:120], ylim=c(0, max(c(Pop$F_m, rep$F_m))), type='l')
lines(rep$F_m, col='blue')


# plot(Pop$Catch_Biomass[years]/mean(Pop$Catch_Biomass[years]), type='l', pch=16)
# lines(rep$predCB/mean(rep$predCB), col='blue')

for (i in 1:12) {
  cal <-  Data$CAL[,1:12][,i]
  plot(Pop$Len_Mids, cal, type='l')
  lines(Pop$Len_Mids, sum(cal)*rep$predCAL[,49:60][,i], col='blue')
}











































# check unfished recruitment
data <- list(model='SLAM',
             Len_Age=Data$Len_Age,
             SD_Len_Age=Data$SD_Len_at_Age,
             Wght_Age=Data$Wt_at_Age,
             Mat_at_Age=Data$Mat_at_Age,
             M_at_Age=Data$M_at_Age,
             PSM_at_Age=Data$phi_at_Age,
             LenBins=Data$Len_Bins,
             LenMids=Data$Len_Mids,
             CAL=Data$CAL,
             CAL_ESS=Data$CAL_ESS,
             Effort=Data$Effort,
             Effort_SD=Data$Effort_SD,
             CPUE=Data$CPUE,
             CPUE_SD=Data$CPUE_SD,
             sigmaRprior=Data$sigmaRprior,
             F_meanprior=Data$F_meanprior,
             Fit_Effort=Data$Fit_Effort,
             Fit_CPUE=Data$Fit_CPUE,
             use_sigmaRprior=Data$use_sigmaRprior,
             use_Fmeanprior=Data$use_Fmeanprior,
             use_Frwpen=Data$use_Frwpen,
             use_R0rwpen=Data$use_R0rwpen,
             log_sigmaR0=Data$log_sigmaR0,
             log_sigmaF=Data$log_sigmaF
)

# starting parameters
# selectivity - 2 parameters
modal_lengths <- data$LenMids[apply(data$CAL, 2, which.max)]
sl50_start <- mean(modal_lengths)
log_sl50 <- log(sl50_start)

modal_lengths <- data$LenMids[apply(data$CAL, 2, which.max)+1]
sl95_start <- mean(modal_lengths)
log_sldelta <- log(sl95_start-sl50_start)

logR0_m_est <- log(Pop$Rec_Pattern[2:12])

nts <- dim(Data$CAL)[2]
nages <- length(Data$Len_Age)
logF_m <- rep(log(0.1), nts)
logF_minit <- log(0.1)
logRec_Devs <- rep(0, nts)
logsigmaR <- log(0.3)

parameters <- list(log_sl50=log_sl50,
                   log_sldelta=log_sldelta,
                   logR0_m_est=logR0_m_est,
                   logsigmaR=logsigmaR,
                   logF_m=logF_m,
                   logF_minit=logF_minit,
                   logRec_Devs=logRec_Devs)

Random <- 'logRec_Devs'

obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, random=Random)

rep <- obj$report()

plot(rep$N_unfished[1,], type='l')
lines(Pop$Rec_Pattern, col='blue')
plot(rep$N_m[1,], type='l')

# To do:
# - clean up assessment model code - options for the priors and penalties
# - write methods for assessment
# - write wrapper functions for assessment and optF
# - simulation test assessment
# - develop MPs
# - closed loop simulation test MPs
# - test with missing data

# Things to follow up on:
# - add option for different selectivity curve
# - make sure first month is January? or make clear in model what the months are

# Scenarios:
# - one set of life-history parameters
# - recruitment: single pulse, spread out, constant
# - available data (1 yr, 3 yrs, 10 yrs):
#   - CAL
#   - CAL and index of effort
#   - CAL, effort, and CPUE

# Assessment performance:
# - recruitment pattern
# - relative error in mean SPR
# - relative error in mean F

# Management Procedures:
# - adjust monthly effort to match predicted
# - adjust overall effort iteratively to match mean F

# Performance:
# - compare yields and utility against perfect management


# Pars$rec_pattern <- rep(1/12)

# Pars$Fvector <- rep(0.1, Pars$nyears*12)


n <- Pop$Number[,years]
n2 <- rep$N_m














rep2 <- obj$report()

plot(rep$F_m[109:120],type='l')
lines(rep2$F_m, col='blue')


par(mfrow=c(2,2))
plot(colSums(rep$N_m),type='b', ylim=c(0, max(colSums(rep$N_m))))
plot(Pop$Rec_Pattern,type='b', ylim=c(0, max(Pop$Rec_Pattern)))
plot(rep$F_m,type='b', ylim=c(0, max(rep$F_m)))
plot(rep$predCB,type='b', ylim=c(0, max(rep$predCB)))

sum(rep$predCB)


Pars$Fvector <- rep(rep$F_m, Pars$nyears)
Pars$Rbar <- 1
Pop <- Simulate(Pars)

lines(Pop$Catch_Biomass[169:180], col='blue')





# Compare with Simulation model
Fmonth <- rep$F_m

Pars$Fvector <- rep(Fmonth, Pars$nyears)
Pars$Rbar <- 1
Pop2 <- Simulate(Pars)

s <- Pop2$SB[169:180]
s2 <- rep$SB

n <- Pop2$Number[,169:180]
n2 <- rep$N_m

plot(n[,1])
lines(n2[,1])

df <- data.frame(s=s, r=n[1,])
df <- df %>% arrange(desc(s))
plot(df, type='l')

df2 <- data.frame(s=s2, r=n2[1,])
df2 <- df2 %>% arrange(desc(s))
lines(df2, col='blue')


# TODO
# - write up - optimizing for yield and utility section
# - set life-history parameters
# - explore properties of simulation model
# - sim test and report results
# - case study






sum(n[,1])
sum(n2[,1])

plot(n[,1])
lines(n2[,1])

plot(Pop2$Catch_Biomass[169:180])
lines(rep$predCB)


par(mfrow=c(2,6))
for (m in 1:12) {
  plot(rep1$N_m[,m], type='l', main=month.abb[m], ylim=c(0, max(rep1$N_m)))
  lines(rep2$N_m[,m], type='l', main=month.abb[m],col='blue')
}

par(mfrow=c(2,6))
for (m in 1:12) {
  plot(rep1$predC_a[,m], type='l', main=month.abb[m], ylim=c(0, max(rep1$predC_a)))
  lines(rep2$predC_a[,m], col='blue')
}



1-sum(rep2$predCB)/sum(rep1$predCB)










# Calculate optimal fishing pattern

# a) maximize yield
# b) maximize utility

# HARA utility - hyperbolic absolute risk aversion

# Write model to calculate best monthly fishing pattern
# given a particular life-history, selectivity pattern, and
# monthly recruitment pattern


