model$opt <- stats::nlminb(model$par, model$fn, model$gr, control = list(iter.max = 1000, eval.max = 1000))
if(n.newton){ # Take a few extra newton steps
  # print("is n.newton")
  tryCatch(for(i in 1:n.newton) {
    g <- as.numeric(model$gr(model$opt$par))
    h <- stats::optimHess(model$opt$par, model$fn, model$gr)
    model$opt$par <- model$opt$par - solve(h, g)
    model$opt$objective <- model$fn(model$opt$par)
  }, error = function(e) {model$err <<- conditionMessage(e)}) # still want fit_tmb to return model if newton steps error out
}


library(SLAM)
library(dplyr)
library(ggplot2)
library(ggthemes)

source('R/functions.R')


# ---- Make Parameters List ----
nyears <- 10
nts <- nyears * 12
nsim <- 100
set.seed(1001)
mulist <- list(1, 6.5, 6.5, c(3,9))
sdlist <- list(100, 1, 3, c(1.5,1.5))
reclist <- list()
for (i in 1:length(mulist)) {
  reclist[[i]] <- GenMonthlyRec(mulist[[i]], sdlist[[i]])
}
rec_scen_names <- c('Constant', 'Pulse', 'Diffuse', 'Bi-Modal')
rec_scen_names <- factor(rec_scen_names, ordered=TRUE, levels=rec_scen_names)
Scens <- list()

calcMonthlyMean <- function(dat, offset=0) {

  nts <- length(dat)
  df <- data.frame(t=1:nts, n=dat)
  df$n[df$n==0] <- NA
  df$M <- (df$t+ offset) %% 12
  df$M[df$M ==0] <- 12
  df$Month <- month.abb[df$M ]
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)
  df %>% group_by(Month) %>%
    summarize(mu=mean(n, na.rm=TRUE), sd=sd(log(n), na.rm=T))
}
data <- SLAM::casestudydata
mu_sd <- calcMonthlyMean(data$Effort)
mu_sd$sd <- 0.3

Effort <- matrix(rlnorm(nts*nsim, log(mu_sd$mu)-0.5*mu_sd$sd^2, mu_sd$sd), nrow=nts, ncol=nsim)

Effort_DF <- data.frame(t=1:nts,
                        Effort=as.numeric(Effort),
                        Sim=rep(1:nsim, each=nts),
                        Year=rep(1:nyears, each=12))
Effort_DF$M <- Effort_DF$t %% 12
Effort_DF$M[Effort_DF$M==0] <- 12
Effort_DF$Month <- month.abb[Effort_DF$M]
Effort_DF$Year <- Effort_DF$Year + 2011
Effort_DF$Date <- lubridate::my(paste(Effort_DF$Month, Effort_DF$Year, sep="-"))

for (s in 1:length(rec_scen_names)) {
  Pars <- list()
  Pars$Name <- rec_scen_names[s]
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
  ind <- switch(as.character(rec_scen_names[s]), Constant=1, Pulse=2, Diffuse=3, 'Bi-Modal'=4)
  Pars$rec_mu <- reclist[[ind]]
  Pars$Rbar <- 1E5

  # recruitment deviations
  rec_sd <- 0.9
  rec_devs <- exp(rnorm(nts*nsim, -0.5*rec_sd^2, rec_sd))
  Pars$rec_devs <- matrix(rec_devs, nrow=nts, ncol=nsim)

  # selectivity
  Pars$sA50 <- 5.5
  Pars$sA95 <- 8

  # effort & fishing mortality
  Pars$Effort <- Effort_DF
  Pars$q <- 0.2

  # Observation process
  Pars$nyears <- 5 # number of years of data
  Pars$sigmaI <- 0.2
  Pars$sigmaE <- 0.2

  Pars$CAW_nsamp <- 1000
  Pars$CAW_ESS <- 200

  Scens[[s]] <- Pars
}

Pars <- Scens[[1]]
i <- 1
Pars$rec_devs[,i] <- 1
Pars$sigmaE <- 0.01
Pars$sigmaI <- 0.01
Pars$CAW_ESS <- 5000
Pars$Effort$Effort <- rep(1, length(Pars$Effort$Effort))
SimPop <- Simulate(Pars, sim=i)

SimPop$Catch_Biomass

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
sim_data$Effort_SD <- rep(0.3, ndata_ts)
sim_data$CPUE <- SimPop$Index
sim_data$CPUE_SD <- rep(0.3, ndata_ts)
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

data <- sim_data

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
lsdelta <- log(0.5)
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

