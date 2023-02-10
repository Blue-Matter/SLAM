#' Iniitialize Parameters for SLAM Estimation
#'
#' @param data A data object
#' @param as50 initial value for age at 50% selection
#' @param as95 initial value for age at 95% selection
#' @param Feq_init
#' @param sigmaR
#' @param log_sigmaF
#' @param log_sigmaR0
#'
#' @return
#' @export
#'
Initialize_Parameters <- function(data,
                                  as50=4, as95=6,
                                  Feq_init=0.05,
                                  F_mean=0.1,
                                  sigmaR=0.5,
                                  sigmaF_m=0.4,
                                  sigmaR0=0.5) {

  parameters <- list()
  parameters$ls50 <- log(as50)
  parameters$lsdelta <- log(as95-as50)

  n_ts <- length(data$Month_ind)

  parameters$logF_minit <- log(Feq_init)
  parameters$logF_ts <- rep(log(0.1), n_ts)

  parameters$log_sigmaF_m <- log(sigmaF_m)

  parameters$logR0_m_est <- rep(1/12, 11)
  parameters$log_sigmaR0 <- log(sigmaR0) # sd for random walk penalty for monthly recruitment
  parameters$logRec_Devs <- rep(log(1),  n_ts)
  parameters$log_sigmaR  <- log(sigmaR) # monthly rec dev sd (usually fixed)
  parameters

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
Initialize_Parameters_OM <- function(SimMod,
                                     sigmaR=0.5,
                                     sigmaF_m=0.4,
                                     sigmaR0=0.5) {
  parameters <- list()

  parameters$ls50 <- log(SimMod$Exploitation$SA50)
  parameters$lsdelta <- log(SimMod$Exploitation$SA95-SimMod$Exploitation$SA50)

  n_age <- SimMod$LifeHistory$maxage+1

  Data_Y_M <- SimMod$Data_TS_DF %>% filter(Sim==1) %>% distinct(Year, Month)
  OM_sub <- SimMod$OM_DF %>% filter(Sim==1, Year%in%Data_Y_M$Year, Month%in%Data_Y_M$Month)

  parameters$logF_minit <- log(0.05)
  F_mean <- mean(OM_sub$F_mort)
  parameters$logF_mean <- log(F_mean)

  nts <- SimMod$Data_TS_DF %>% filter(Sim==1) %>% nrow()
  parameters$logF_ts_dev <- rep(log(1), nts)

  parameters$log_sigmaF_m <- log(sigmaF_m)

  parameters$logR0_m_est <- log(SimMod$LifeHistory$R0_m[2:12])
  parameters$log_sigmaR0 <- log(sigmaR0) # sd for random walk penalty for monthly recruitment
  parameters$logRec_Devs <- log(OM_sub$Rec_Devs)
  parameters$log_sigmaR  <- log(SimMod$LifeHistory$sigmaR) # monthly rec dev sd (usually fixed)
  parameters
}



#' Title
#'
#' @param sim
#' @param SimMod
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
Construct_Data_OM <- function(sim=1,
                              SimMod,
                              CAW_Monthly_ESS=100,
                              Effort_CV=0.2,
                              CPUE_CV=0.2,
                              Fit_Effort=1,
                              Fit_CPUE=1,
                              Fit_CAW=1,
                              use_Frwpen=1,
                              use_R0rwpen=1) {

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

  # effort_by_year <- split(data$Effort, ceiling(seq_along(data$Effort)/12))
  # data$Effort_y_mean <- as.vector(unlist(lapply(effort_by_year, mean)))

  # Options
  data$Fit_Effort <- Fit_Effort
  data$Fit_CPUE <- Fit_CPUE
  data$Fit_CAW <- Fit_CAW

  # Penalties
  data$use_Frwpen <- use_Frwpen
  data$use_R0rwpen <- use_R0rwpen

  data$model <- 'SLAM'
  data$currentYr <- max(SimMod$OM_DF$Year)
  data$Month_ind <- 1:length(data$Effort)

  data
}



#' Title
#'
#' @param data
#' @param as50
#' @param as95
#' @param Feq_init
#' @param sigmaR
#' @param control
#' @param map
#' @param log_sigmaF
#' @param log_sigmaR0
#'
#' @return
#' @export
#'
Do_Assess <- function(data,
                   as50=4,
                   as95=6,
                   Feq_init=0.05,
                   sigmaR=0.5,
                   control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20),
                   map=list(log_sigmaF=factor(NA),
                            log_sigmaR=factor(NA),
                            log_sigmaR0=factor(NA)),
                   log_sigmaF = log(0.2),
                   log_sigmaR0 = log(0.1)) {

  parameters <- Initialize_Parameters(data,
                                      as50, as95,
                                      Feq_init,
                                      sigmaR,
                                      log_sigmaF,
                                      log_sigmaR0)

  if (data$n_years==1) {
    map$logRec_Devs <- rep(factor(NA), length(parameters$logRec_Devs))
  }


  if (!is.null(map$log_sigmaR)) {
    Random <- NULL
  } else {
    Random <- 'logRec_Devs'
  }
  do_opt <- opt_TMB_model(data, parameters, map, Random, control, restarts=10)
  do_opt
}


#' Title
#'
#' @param data
#' @param parameters
#' @param map
#' @param Random
#' @param control
#' @param restarts
#'
#' @return
#' @export
#'
opt_TMB_model <- function(data, parameters, map, Random, control, restarts=10) {

  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, hessian=FALSE, map=map, random=Random)

  starts <- obj$par
  opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control)),silent=TRUE)

  rerun <- FALSE
  if (inherits(opt, 'list')) {
    rep <- obj$report()
    sdreport <- TMB::sdreport(obj)

    # check convergence, gradient and positive definite
    chk <- data.frame(#pdHess=sdreport$pdHess,
                      conv=opt$convergence ==0,
                      grad=max(abs(sdreport$gradient.fixed)) < 0.1)

    if (any(!chk)) rerun <- TRUE
  } else {
    chk <- opt
    rerun <- TRUE
    sdreport <- NA
    rep <- NA
  }

  if (rerun & restarts>0) {
    parameters$ls50 <- parameters$ls50 * rnorm(1, 1, 0.1)
    parameters$lsdelta <- parameters$lsdelta  * rnorm(1, 1, 0.4)
    parameters$logR0_m_est <- rep(log(1), 11) + rnorm(11, 0, 0.4)
    Recall(data, parameters,  map, Random, control, restarts-1)
  }
  list(opt=opt, obj=obj, rep=rep, sdreport=sdreport, chk=chk)
}

#' Title
#'
#' @param data
#' @param parameters
#' @param map
#' @param Random
#' @param control
#' @param restarts
#'
#' @return
#' @export
#'
Compare_OM_Assess <- function(x, SimMod, Assessment) {

  # OM
  Mod <- SimMod$OM_DF %>% filter(Sim==x)
  nyears <- SimMod$OM_DF$Year %>% unique() %>% length()
  n_recent_months <- SimMod$Data_TS_DF$Month_ind %>% unique() %>% length()
  month_ind <- rev(seq(nyears*12, by=-1, length.out=n_recent_months))
  Mod_data <- SimMod$OM_DF %>% filter(Sim==x, Month_ind%in%month_ind)
  Mod_data$Month_ind <- 1:nrow(Mod_data)
  Mod_data$Model <- 'OM'
  Mod_data$Year_Month <- paste(Mod_data$Year, Mod_data$Month, sep="_")

  # Assessment
  nMonths <- length(Assessment$obj$env$data$CPUE)
  nYears <- nMonths/12
  currentYr <- Assessment$obj$env$data$currentYr
  Years <- (currentYr-nYears+1):currentYr
  Years <- rep(Years, each=12)
  Months <- rep(month.abb, nYears)

  Est_TS_DF <- data.frame(Sim=x,
                          Year=Years,
                          Month=Months,
                          Month_ind=1:nMonths,
                          SB_fished=Assessment$rep$SB_m,
                          N_fished=apply( Assessment$rep$N_m,2,sum),
                          Catch= Assessment$rep$predCB,
                          SPR= Assessment$rep$SPR,
                          Effort= Assessment$rep$StEffort,
                          F_mort= Assessment$rep$F_m,
                          Rec_Devs=exp(Assessment$rep$logRec_Devs),
                          Index= Assessment$rep$stpredIndex,
                          Model='Assessment'
  )
  Est_TS_DF$Year_Month <- paste(Est_TS_DF$Year, Est_TS_DF$Month, sep="_")
  bind_rows(Mod_data, Est_TS_DF)
}

