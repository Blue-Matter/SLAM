library(SLAM)


Pars <- Load_Scenario(1)

Pop <- Simulate(Pars)

nyears <- dim(Pop$CAL_samp)[2]
use_years <- 10 # use the last use_years

years <- (nyears-use_years*12+1):nyears
nts <- length(years)



data <- list(model='optFpattern',
             rec_pattern=Pop$Rec_Pattern,
             opt_type=0,
             utilpow=1,
             Len_Age=Pop$Len_at_Age,
             SD_Len_Age=Pop$SD_Len_at_Age,
             Wght_Age=Pop$Wt_at_Age,
             Mat_at_Age=Pop$Mat_at_Age,
             M_at_Age=Pop$M_at_Age,
             PSM_at_Age=Pop$phi_at_Age,
             LenBins=Pop$Len_Bins,
             LenMids=Pop$Len_Mids)

logF_m <- rep(log(0.3), nts)
parameters <- list(logF_m=logF_m)

obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE)







# Quick test
CAL <- Pop$CAL_samp[,years]
CAL_ESS <- rep(100, nts)

Effort <- Pop$Effort[years]/mean(Pop$Effort[years])
Effort_SD <- rep(0.1, nts)
EffExists <- rep(1, nts)
nEffMonths <- nts
CPUE <- Pop$Index[years]/mean(Pop$Index[years])
CPUE_SD <- rep(0.2, nts)

sigmaRprior <- c(0.5, 0.2)
F_meanprior <- c(0.4, 0.3)
Fit_Effort <- 1
Fit_CPUE <- 1

data <- list(model='SLAM',
             Len_Age=Pop$Len_at_Age,
             SD_Len_Age=Pop$SD_Len_at_Age,
             Wght_Age=Pop$Wt_at_Age,
             Mat_at_Age=Pop$Mat_at_Age,
             M_at_Age=Pop$M_at_Age,
             PSM_at_Age=Pop$phi_at_Age,
             LenBins=Pop$Len_Bins,
             LenMids=Pop$Len_Mids,
             CAL=CAL,
             CAL_ESS=CAL_ESS,
             Effort=Effort,
             Effort_SD=Effort_SD,
             EffExists=EffExists,
             nEffMonths=nEffMonths,
             CPUE=CPUE,
             CPUE_SD=CPUE_SD,
             sigmaRprior=sigmaRprior,
             F_meanprior=F_meanprior,
             Fit_Effort=Fit_Effort,
             Fit_CPUE=Fit_CPUE)

log_sl50 <- log(Pars$SL5)
log_sldelta <- log(1)
logR0_m <- rep(log(1), 12)
logsigmaR <- log(0.2)
logF_m <- rep(log(0.3), nts)
logRec_Devs <- rep(0, nts)

log_sigmaF <- log(0.3)
log_sigmaR0 <- log(0.5)

parameters <- list(log_sl50=log_sl50,
                   log_sldelta=log_sldelta,
                   logR0_m=logR0_m,
                   log_sigmaR0=log_sigmaR0,
                   logsigmaR=logsigmaR,
                   logF_m=logF_m,
                   log_sigmaF=log_sigmaF,
                   logRec_Devs=logRec_Devs)

Random <- 'logRec_Devs'

obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, random=Random)

starts <- obj$par
opt <- suppressWarnings(nlminb(starts, obj$fn, obj$gr))


rep <- obj$report()

plot(Pop$Rec_Pattern,type='l', ylim=c(0, max(Pop$Rec_Pattern)))
lines(rep$R0_m, col='blue')


# Calculate optimal fishing pattern

# a) maximize yield
# b) maximize utility

# HARA utility - hyperbolic absolute risk aversion

# Write model to calculate best monthly fishing pattern
# given a particular life-history, selectivity pattern, and
# monthly recruitment pattern


