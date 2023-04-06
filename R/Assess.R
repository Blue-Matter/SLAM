


#' Initialize Parameters for SLAM Estimation
#'
#' @param data An object of class `Data`
#' @param as50 initial value for age at 50% selection
#' @param as95 initial value for age at 95% selection
#' @param Feq_init initial value for initial equilibrium fishing mortality
#' @param F_ts initial values for monthly fishing mortality. Must be length `length(Data$Years)`
#' @param sigmaR standard deviation of log-normal recruitment deviations
#' @param sigmaF_m standard deviation of random-walk penalty for monthly fishing mortality
#' @param sigmaR0 standard deviation of random-walk penalty for monthly (seasonal) recruitment
#'
#' @return A list of class `Parameters`
#' @export
Initialize_Parameters <- function(data,
                                  as50=4, as95=6,
                                  Feq_init=0.15,
                                  F_ts=0.1,
                                  sigmaR=0.6,
                                  sigmaF_m=0.4,
                                  sigmaR0=0.1) {

  if (!inherits(data, 'Data'))
    stop('First argument must be object of class `Data`')

  parameters <- New_Parameters()
  parameters$ls50 <- log(as50)
  parameters$lsdelta <- log(as95-as50)

  n_ts <- ncol(data$CAW)

  parameters$logF_minit <- log(Feq_init)
  parameters$logF_ts <- rep(log(F_ts), n_ts)

  parameters$log_sigmaF_m <- log(sigmaF_m)

  parameters$logR0_m_est <- rep(0, 11)
  parameters$log_sigmaR0 <- log(sigmaR0) # sd for random walk penalty for monthly recruitment
  parameters$logRec_Devs <- rep(log(1),  n_ts)
  parameters$log_sigmaR  <- log(sigmaR) # monthly rec dev sd (usually fixed)
  class(parameters) <- 'Parameters'
  parameters
}



Initialize_Parameters_OM <- function(Simulation, Data, sim=1) {

  parameters <- list()
  parameters$ls50 <- log(Simulation$Exploitation$SA50)
  parameters$lsdelta <- log(Simulation$Exploitation$SA95-Simulation$Exploitation$SA50)

  OM <- Simulation$Time_Series %>% filter(Sim==sim)
  n_ts <- length(OM$Month_ind)
  n_ts_data <- Data$Effort %>% length()
  ts <- (n_ts-n_ts_data+1):n_ts

  # Calculate initial equilibrium F
  OM_at_Age <- Simulation$At_Age_Time_Series  %>% filter(Sim==sim)
  init_data_ts <- OM_at_Age %>% filter(Month_ind==ts[1])
  init_ts <- OM_at_Age %>% filter(Month_ind==1)
  nage <- 15
  Zs <- -log(init_data_ts$N_fished[2:nage]/init_ts$N_unfished_eq[1:(nage-1)])
  Zs <- c(max(Simulation$LifeHistory$M_at_Age), Zs)
  Fs <- Zs - Simulation$LifeHistory$M_at_Age
  Fs <- Fs[is.finite(Fs)]

  parameters$logF_minit <- log(mean(Fs))
  parameters$logF_ts <- log(OM$F_mort[ts])

  parameters$log_sigmaF_m <- log(5)
  parameters$logR0_m_est <- log(Simulation$LifeHistory$R0_m[2:12]/mean(Simulation$LifeHistory$R0_m[2:12]))
  parameters$log_sigmaR0 <- log(5) # sd for random walk penalty for monthly recruitment
  parameters$logRec_Devs <- log(OM$Rec_Devs[ts])
  parameters$log_sigmaR  <- log(Simulation$LifeHistory$sigmaR) # monthly rec dev sd (usually fixed)
  parameters
}



set_data_types <- function(data, Data_types) {
  Data_types <- strsplit(Data_types, "\\+")[[1]]
  if (!'CAW' %in% Data_types) data$Fit_CAW <-0
  if (!'CAA' %in% Data_types) data$Fit_CAA <-0
  if (!'Index' %in% Data_types) data$Fit_Index <-0
  if (!'Effort' %in% Data_types) data$Fit_Effort <-0
  data
}



Construct_Data_OM <- function(sim=1,
                              SimMod,
                              CAW_Monthly_ESS=100,
                              CAA_Monthly_ESS=100,
                              Effort_CV=0.2,
                              CPUE_CV=0.2,
                              Fit_Effort=1,
                              Fit_CPUE=1,
                              Fit_CAW=1,
                              Fit_CAA=1,
                              use_Frwpen=1,
                              use_R0rwpen=1,
                              Data_types=NULL
                              ) {

  data <- list()
  # Assumed life-history parameters
  data$Weight_Age <- SimMod$LifeHistory$Weight_Age_Mean
  data$Weight_Age_SD <- SimMod$LifeHistory$Weight_Age_SD
  data$Mat_at_Age <- SimMod$LifeHistory$Maturity_at_Age
  data$M_at_Age <- SimMod$LifeHistory$M_at_Age
  data$PSM_at_Age <- SimMod$LifeHistory$Post_Spawning_Mortality
  data$h <- SimMod$LifeHistory$steepness
  nAges <- length(data$Weight_Age)

  # CAW Data
  data$WghtBins <- SimMod$Data$Weight_Bins
  data$WghtMids <- SimMod$Data$Weight_Mids

  CAW_DF <- SimMod$Data_CAW_DF %>% filter(Sim==sim) %>% select(Month_ind, Weight, Count)
  nBins <- length(unique(CAW_DF$Weight))
  nMonths <- length(unique(CAW_DF$Month_ind))

  CAW <- matrix(CAW_DF$Count, nrow=nBins, nMonths)
  data$CAW <- CAW
  data$CAW_ESS <- rep(CAW_Monthly_ESS, nMonths)

  # CAA Data
  CAA_DF <- SimMod$Data_CAA_DF %>% filter(Sim==sim) %>% select(Month_ind, Age, Count)
  CAA <- matrix(CAA_DF$Count, nrow=nAges, nMonths)
  data$CAA <- CAA
  data$CAA_ESS <- rep(CAA_Monthly_ESS, nMonths)

  # Effort Index
  Effort_DF <- SimMod$Data_TS_DF %>% filter(Sim==sim) %>% select(Effort)
  data$Effort <- Effort_DF$Effort
  data$Effort_SD <- rep(Effort_CV, nMonths)

  # CPUE Index
  CPUE_DF <- SimMod$Data_TS_DF %>% filter(Sim==sim) %>% select(CPUE)
  data$CPUE <- CPUE_DF$CPUE
  data$CPUE_SD <- rep(CPUE_CV, nMonths)

  data$n_years <- floor(nMonths/12)

  # Options
  data$Fit_Effort <- Fit_Effort
  data$Fit_CPUE <- Fit_CPUE
  data$Fit_CAW <- Fit_CAW
  data$Fit_CAA <- Fit_CAA

  # Penalties
  data$use_Frwpen <- use_Frwpen
  data$use_R0rwpen <- use_R0rwpen

  data$model <- 'SLAM'
  data$currentYr <- max(SimMod$OM_DF$Year)
  data$Month_ind <- 1:length(data$Effort)

  if (!is.null(Data_types)) {
    data <- set_data_types(data, Data_types=Data_types)
  }
  data
}



#' Conduct an Assessment
#'
#' @param Data A `Data` object
#' @param Parameters A `Parameters` object
#' @param Assumed_h Assumed value of steepness for the stock-recruitment relationship
#' @param max_ESS Maximum effective sample size for the catch-at-weight data
#' @param Est_Rec_Devs Logical. Estimate recruitment deviations (TRUE) or assume constant recruitment (FALSE)
#' @param control Optional controls for optimizer
#' @param ... Additional arguments pass to TMB
#'
#' @return A list of class `Assess`
#' @export
#'
Assess <- function(Data, Parameters=NULL,
                   Assumed_h=0.7,
                   max_ESS=200,
                   Est_Rec_Devs=ifelse(length(Data$Year)>=24, TRUE, FALSE),
                   control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20),
                   ...) {

  if (!inherits(Data, 'Data'))
    stop('First argument must be object of class `Data`')

  if (is.null(Parameters)) {
    Parameters <- Initialize_Parameters(Data)
  }

  Check_Parameters(Parameters, Data)

  dots <- list(...)
  if (!is.null(dots$map)) {
    map <- dots$map
  } else {
    map <- list(log_sigmaF_m=factor(NA),
                log_sigmaR=factor(NA),
                log_sigmaR0=factor(NA))
  }

  if (!Est_Rec_Devs) {
    map$logRec_Devs <- rep(factor(NA), length(Parameters$logRec_Devs))
    Parameters$log_sigmaR <- log(0.001)
  }

  if (!is.null(map$log_sigmaR)) {
    Random <- NULL
  } else {
    Random <- 'logRec_Devs'
  }

  Data$h <- Assumed_h
  message('Assuming a BH-SRR steepness of ', Assumed_h)

  Data$CAW_ESS[Data$CAW_ESS>max_ESS] <- max_ESS
  outData <- Data

  Data <- Update(Data)
  # drop
  Data$Year <- Data$Month <- Data$Metadata <- NULL

  do_opt <- opt_TMB_model(Data, Parameters, map, Random, control, restarts=10)
  do_opt$map <- map
  do_opt$Data <- outData
  do_opt$Parameters <- Parameters
  class(do_opt) <- 'Assess'
  do_opt

}



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

