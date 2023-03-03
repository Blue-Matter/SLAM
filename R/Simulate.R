
#' Simulate the fishery dynamics for a short-lived semelparous species
#'
#'
#' @param LifeHistory A named `list` object with the life-history parameters,
#' or the file path to a correctly
#' formatted CSV file containing the `LifeHistory` (and optionally `Exploitation`) parameters.
#' The `LifeHistory` object can optionally include the `Exploitation` parameters as well.
#' @param Exploitation An named `list` object with the exploitation parameters or
#'  the file path to a correctly formatted `Exploitation` CSV file.
#'  Not needed if `LifeHistory` already includes  the `Exploitation` parameters.
#' @param nsim The number of stochastic simulations. Default is 100.
#' @param seed The seed for the random number generator
#' @param currentYr The current (most recent) year for the simulations
#' @param silent Hide messages?
#'
#' @details See the example `SLAM::LifeHistory` and `SLAM::Exploitation` objects and help
#' documentation (`?LifeHistory` and `?Exploitation`) for examples.
#'
#' @return A named list of class `Simulation` containing:
#' \itemize{
#'   \item LifeHistory: The `LifeHistory` Object used for the simulations
#'   \item Exploitation: The `Exploitation` Object used for the simulations, with the addition of `Sel_at_Age` (selectivity schedule)
#'   \item At_Age_Time_Series: A data frame with the following columns
#'     \itemize{
#'     \item Sim: The simulation number
#'     \item Age: The age class
#'     \item Year: The year
#'     \item Month: Abbreviated calendar month
#'     \item Month_ind: Numerical index for the monthly timestep
#'     \item N_fished: The number of fished animals in each age class for each timestep
#'     \item N_unfished: The number of unfished animals in each age class for each timestep
#'     \item N_unfished_eq: The equilibrium number of unfished animals in each age class for each timestep
#'     \item B_fished: The biomass of fished animals in each age class for each timestep
#'     \item B_unfished: The biomass of unfished animals in each age class for each timestep
#'     \item SB_fished: The spawning biomass of fished animals in each age class for each timestep
#'     \item SB_unfished: The spawning biomass of unfished animals in each age class for each timestep
#'     \item Catch: The catch biomass in each age class for each timestep
#'     \item Catch_n: The catch numbers in each age class for each timestep
#'     }
#'   \item Time_Series: A data frame with the following columns
#'     \itemize{
#'     \item Sim: The simulation number
#'     \item Year: The year
#'     \item Month: Abbreviated calendar month
#'     \item Month_ind: Numerical index for the monthly timestep
#'     \item N_fished: The total number of fished animals in each timestep
#'     \item N_unfished: The total number of unfished animals in each timestep
#'     \item B_fished: The biomass of fished animals in each timestep
#'     \item B_unfished: The biomass of unfished animals in each timestep
#'     \item SB_fished: The spawning biomass of fished animals in each timestep
#'     \item SB_unfished: The spawning biomass of unfished animals in each timestep
#'     \item Catch: The catch biomass in each timestep
#'     \item Recruits: The number of recruits in each timestep
#'     \item SPR: The equilibrium spawning potential ratio in each timestep
#'     \item Effort: The fishing effort in each timestep
#'     \item F_mort: The fishing mortality in each timestep
#'     \item Rec_Devs: The recruitment deviation in each timestep
#'     }
#' }
#' @export
Simulate <- function(LifeHistory=NULL,
                     Exploitation=NULL,
                     nsim=100,
                     seed=101,
                     currentYr=format(Sys.Date(), "%Y"),
                     silent=FALSE) {

  set.seed(seed)

  # Import the parameters
  if (is.null(LifeHistory)) {
    if (!silent) message('Using Example Life History Parameters from `SLAM::LifeHistory`')
    LifeHistory <- SLAM::LifeHistory
  }

  # Read in CSVs
  if (inherits(LifeHistory, 'character')) {
    if (!silent) message('Reading Life History Parameters from: ', LifeHistory)
    Parameters <- Import_Parameters(LifeHistory)
    LifeHistory <- Parameters$LifeHistory
    Exploitation <- Parameters$Exploitation
    if (!is.null(Exploitation))
      if (!silent) message('Exploitation parameters also detetected and imported')
  }

  if (!is.null(Exploitation) & inherits(Exploitation, 'character')) {
    message('Reading Exploitation Parameters from: ', Exploitation)
    Exploitation <- Parameters$Exploitation
  }

  if (is.null(Exploitation)) {
    if (!silent) message('Using Example Exploitation Parameters from `SLAM::Exploitation`')
    Exploitation <- SLAM::Exploitation
  }

  # Check Objects
  import_nms <- names(LifeHistory)
  missing <- names(SLAM::LifeHistory)[!names(SLAM::LifeHistory) %in% import_nms]
  if (length(missing)>0) {
    stop('`LifeHistory` object is missing: ', paste(missing, collapse=', '))
  }
  import_nms <- names(Exploitation)
  missing <- names(SLAM::Exploitation)[!names(SLAM::Exploitation) %in% import_nms]
  if (length(missing)>0) {
    stop('`Exploitation` object is missing: ', paste(missing, collapse=', '))
  }


  # Import Life-History Parameters
  Stock_Name <- LifeHistory$Stock_Name
  Species <- LifeHistory$Species
  Ages <- LifeHistory$Ages
  maxage <- max(Ages)
  nAge <- length(Ages)
  Weight_Age_Mean <- LifeHistory$Weight_Age_Mean
  Weight_Age_SD <- LifeHistory$Weight_Age_SD
  M_at_Age <- LifeHistory$M_at_Age
  Maturity_at_Age <- LifeHistory$Maturity_at_Age
  Post_Spawning_Mortality <- LifeHistory$Post_Spawning_Mortality

  R0_bar <- LifeHistory$R0_bar
  R0_m <- LifeHistory$R0_m
  R0_m <- R0_m/sum(R0_m)
  sigmaR <- LifeHistory$sigmaR
  h <- LifeHistory$steepness
  if (is.null(h)) h <- 0.999
  if (is.na(h)) h <- 0.999

  # Import Exploitation Parameters
  Fleet_Name <- Exploitation$Fleet_Name
  SA50 <- Exploitation$SA50
  SA95 <- Exploitation$SA95
  nmonths <- Exploitation$nmonths
  Effort_pattern <- Exploitation$Effort_pattern
  Effort_month <- Exploitation$Effort_month
  Effort_current <- Exploitation$Effort_current
  Effort_cv <- Exploitation$Effort_cv
  q <- Exploitation$q
  q_cv <- Exploitation$q_cv
  HARA_power <- Exploitation$HARA_power
  nts <- nmonths

  # Parameter Checks
  #TODO

  # Model Fishery Dynamics

  # Generate Monthly Recruitment Deviations
  Rec_Devs <- exp(rnorm((nts+maxage)*nsim, -0.5*sigmaR^2, sigmaR))
  Rec_Devs <- matrix(Rec_Devs, nrow=nts+maxage, ncol=nsim)
  Rec_Devs <- t(Rec_Devs)

  # Generate Fishing Effort and Mortality
  qs_dev <- rlnorm(nmonths*nsim, -0.5*q_cv^2, q_cv)
  qs <- q * matrix(qs_dev, nrow=nsim, ncol=nmonths)

  if (is.null(Effort_month)) {
    Effort_month <- generate_Effort(Exploitation, LifeHistory, nsim)
  }
  F_Month <- Effort_month * qs
  if (any(F_Month[,1]>0)) {
    if (!silent)
      message('`Exploitation` parameters had positive fishing effort in initial month. \nChanging so first month is unfished (F=0)')
    F_Month[,1] <- 0
  }

  # Set up Arrays
  M_at_Age <- replicate(nsim, LifeHistory$M_at_Age)
  M_at_Age <- replicate(nts, M_at_Age)
  M_at_Age <- aperm(M_at_Age, c(2,1,3))
  F_at_Age <- Z_at_Age <- array(0, dim=dim(M_at_Age))

  N_Age_fished <- array(0, dim=c(nsim, nAge, nts))
  N_Age_unfished <- N_Age_unfished_eq <- SPR <- SB_Age_fished <- SB_Age_unfished <-
    SB_Age_unfished_eq <- Catch_Age <- N_Age_fished
  B_Age_fished <- B_Age_unfished <- B_Age_unfished_eq <- array(0, dim=c(nsim, nAge, nts))
  SB_fished <- SB_unfished <- SB_unfished_eq <- matrix(0, nrow=nsim, ncol=nts)

  Sel_at_Age <-  1/(1+exp(-log(19)*((Ages-SA50)/(SA95-SA50))))
  Sel_at_Age <- replicate(nsim, Sel_at_Age)
  Sel_at_Age <- replicate(nts, Sel_at_Age)
  Sel_at_Age <- aperm(Sel_at_Age, c(2,1,3))

  ind <- as.matrix(expand.grid(1:nsim, 1:nAge, 1:nts))
  F_at_Age[ind] <- F_Month[ind[,c(1,3)]] * Sel_at_Age[ind]
  Z_at_Age <- M_at_Age + F_at_Age

  # Calculate Spawning Potential Ratio
  theta_0 <- array(1, dim=c(nsim, nAge))
  theta_F <- array(0, dim=c(nsim, nAge,nts))
  theta_F[,1,] <- 1
  for (a in 2:nAge) {
    theta_0[,a] <- theta_0[,a-1]*exp(-M_at_Age[a-1])*(1-Post_Spawning_Mortality[a-1])
  }
  for (t in 1:nts) {
    for (a in 2:nAge) {
      theta_F[,a,t] <- theta_F[,a-1,t]*exp(-(Z_at_Age[,a-1,t]))*(1-Post_Spawning_Mortality[a-1])
    }
  }

  Maturity_at_Age_array <- t(replicate(nsim, Maturity_at_Age))
  Weight_Age_Mean_array <- t(replicate(nsim, Weight_Age_Mean))
  Post_Spawning_Mortality_array <- t(replicate(nsim, Post_Spawning_Mortality))
  E0 <- apply(theta_0*Maturity_at_Age_array*Weight_Age_Mean_array, 1, sum)
  E0 <- replicate(nts, E0)

  Maturity_at_Age_array <- replicate(nts, Maturity_at_Age_array)
  Weight_Age_Mean_array <- replicate(nts, Weight_Age_Mean_array)
  Post_Spawning_Mortality_array <- replicate(nts, Post_Spawning_Mortality_array)
  Ef <- apply(theta_F*Maturity_at_Age_array*Weight_Age_Mean_array,c(1,3), sum)
  SPR <- Ef/E0

  # Initialize unfished population
  N_Age_unfished[,1,1] <-  R0_bar * R0_m[1] * Rec_Devs[,1]
  N_Age_unfished_eq[,1,1] <-  R0_bar * R0_m[1]
  for (a in 2:nAge) {
    age <- a -1
    if (age %% 12==0) {
      m <- -(age%%12)+1
    } else {
      m <- -(age %% 12)+13
    }

    cumM <- apply(M_at_Age[,1:(a-1),1, drop=FALSE], 1, sum)
    cumPSM <- prod(1-Post_Spawning_Mortality[1:(a-1)])

    N_Age_unfished[,a,1] <-  R0_bar * R0_m[m]* exp(-cumM) * cumPSM  * Rec_Devs[,a]
    N_Age_unfished_eq[,a,1] <-  R0_bar * R0_m[m]* exp(-cumM) * cumPSM
  }

  N_Age_fished <- N_Age_unfished
  SB_Age_fished[,,1] <- N_Age_fished[,,1] * Weight_Age_Mean_array[,,1] * Maturity_at_Age_array[,,1] * exp(-F_Month[,1]/2)
  SB_Age_unfished[,,1] <- N_Age_unfished_eq[,,1] * Weight_Age_Mean_array[,,1] * Maturity_at_Age_array[,,1]
  SB_Age_unfished_eq[,,1] <- N_Age_unfished_eq[,,1] * Weight_Age_Mean_array[,,1] * Maturity_at_Age_array[,,1]
  SB_fished[,1] <- apply(SB_Age_fished[,,1], 1, sum)
  SB_unfished[,1] <- apply(SB_Age_unfished[,,1], 1, sum)
  SB_unfished_eq[,1] <- apply(SB_Age_unfished_eq[,,1], 1, sum)
  Catch_Age[,,1] <- N_Age_fished[,,1]*((1-Maturity_at_Age_array[,,1])*exp(-M_at_Age[,,1]/2)+Maturity_at_Age_array[,,1]*exp(-Post_Spawning_Mortality_array[,,1]/2))*(1-exp(-F_at_Age[,,1]))

  B_Age_fished[,,1] <- N_Age_fished[,,1] * Weight_Age_Mean_array[,,1]
  B_Age_unfished[,,1]  <- N_Age_unfished[,,1] * Weight_Age_Mean_array[,,1]
  B_Age_unfished_eq[,,1]  <- N_Age_unfished_eq[,,1] * Weight_Age_Mean_array[,,1]

  # Fished population - Loop over remaining time-steps
  for (t in 2:nts) {
    month <- t%%12
    if (month==0) month <-12
    for (a in 1:maxage) {
      N_Age_fished[,a+1, t] <- N_Age_fished[,a,t-1]*exp(-Z_at_Age[,a,t-1])*(1-Post_Spawning_Mortality_array[,a,t-1])
      N_Age_unfished[,a+1, t] <- N_Age_unfished[,a,t-1]*exp(-M_at_Age[,a,t-1])*(1-Post_Spawning_Mortality_array[,a,t-1])
      N_Age_unfished_eq[,a+1, t] <- N_Age_unfished_eq[,a,t-1]*exp(-M_at_Age[,a,t-1])*(1-Post_Spawning_Mortality_array[,a,t-1])
    }

    B_Age_fished[,,t] <- N_Age_fished[,,t] * Weight_Age_Mean_array[,,t]
    B_Age_unfished[,,t]  <- N_Age_unfished[,,t] * Weight_Age_Mean_array[,,t]
    B_Age_unfished_eq[,,t]  <- N_Age_unfished_eq[,,t] * Weight_Age_Mean_array[,,t]

    # Spawning Biomass mid time-step
    SB_Age_fished[,,t] <- N_Age_fished[,,t] * Weight_Age_Mean_array[,,t] * Maturity_at_Age_array[,,t] * exp(-F_Month[,t]/2)
    SB_Age_unfished[,,t] <- N_Age_unfished[,,t] * Weight_Age_Mean_array[,,t] * Maturity_at_Age_array[,,t]
    SB_Age_unfished_eq[,,t] <- N_Age_unfished_eq[,,t] * Weight_Age_Mean_array[,,t] * Maturity_at_Age_array[,,t]
    SB_fished[,t] <- apply(SB_Age_fished[,,t], 1, sum)
    SB_unfished[,t] <- apply(SB_Age_unfished[,,t], 1, sum)
    SB_unfished_eq[,t] <- apply(SB_Age_unfished_eq[,,t], 1, sum)

    # Catch during time-step
    Catch_Age[,,t] <- N_Age_fished[,,t]*((1-Maturity_at_Age_array[,,t])*exp(-M_at_Age[,,t]/2)+Maturity_at_Age_array[,,t]*exp(-Post_Spawning_Mortality_array[,,t]/2))*(1-exp(-F_at_Age[,,t]))

    # Recruitment to age-0 this time-step
    N_Age_fished[,1, t] <- sapply(1:nsim, function(x)
      BH_SRR(R0=R0_bar*R0_m[month], h=h, SB=SB_fished[x,t], SBpR=E0[x,t]) * Rec_Devs[x,t])

    N_Age_unfished[,1, t] <- sapply(1:nsim, function(x)
      BH_SRR(R0=R0_bar*R0_m[month], h=h, SB=SB_unfished[x,t], SBpR=E0[x,t]) * Rec_Devs[x,t])

    N_Age_unfished_eq[,1, t] <- sapply(1:nsim, function(x)
      BH_SRR(R0=R0_bar*R0_m[month], h=h, SB=SB_unfished_eq[x,t], SBpR=E0[x,t]))
  }

  # Calculate Monthly Time-series
  B_unfished_Month <- apply(N_Age_unfished * Weight_Age_Mean_array, c(1,3), sum)
  B_unfished_eq_Month <- apply(N_Age_unfished_eq * Weight_Age_Mean_array, c(1,3), sum)
  SB_unfished_Month <- apply(N_Age_unfished * Weight_Age_Mean_array * Maturity_at_Age_array, c(1,3), sum)

  B_Month <- apply(N_Age_fished * Weight_Age_Mean_array, c(1,3), sum)
  Catch_N_Month <- apply(Catch_Age, c(1,3), sum)
  Catch_B_Month <- apply(Catch_Age * Weight_Age_Mean_array, c(1,3), sum)

  N_unfished <- apply(N_Age_unfished, c(1,3), sum)
  N_fished <- apply(N_Age_fished, c(1,3), sum)


  # Return info
  currentYr <- as.numeric(currentYr)
  nyears <- nmonths/12
  Years <- (currentYr-nyears+1):currentYr
  Years <- rep(Years, each=12)
  Months <- rep(month.abb, nyears)

  Exploitation$Sel_at_Age <- Sel_at_Age

  At_Age_Time_Series <- data.frame(Sim=1:nsim,
                                   Age=rep(Ages, each=nsim),
                                   Year=rep(Years, each=nsim*nAge),
                                   Month=rep(Months, each=nsim*nAge),
                                   Month_ind=rep(1:nts, each=nsim*nAge),
                                   N_fished=as.vector(N_Age_fished),
                                   N_unfished=as.vector(N_Age_unfished),
                                   N_unfished_eq=as.vector(N_Age_unfished),
                                   B_fished=as.vector(B_Age_fished),
                                   B_unfished=as.vector(N_Age_unfished * Weight_Age_Mean_array),
                                   SB_fished=as.vector(SB_Age_fished),
                                   SB_unfished=as.vector(SB_Age_unfished),
                                   Catch=as.vector(Catch_Age* Weight_Age_Mean_array),
                                   Catch_n=as.vector(Catch_Age)
                                   )

  Time_Series <- data.frame(Sim=1:nsim,
                            Year=rep(Years, each=nsim),
                            Month=rep(Months, each=nsim),
                            Month_ind=rep(1:nts, each=nsim),
                            N_fished=as.vector(N_fished),
                            N_unfished=as.vector(N_unfished),
                            B_fished=as.vector(B_Month),
                            B_unfished=as.vector(B_unfished_Month),
                            B_unfished_eq=as.vector(B_unfished_eq_Month),
                            SB_fished=as.vector(SB_fished),
                            SB_unfished=as.vector(SB_unfished_Month),
                            SB_unfished_eq=as.vector(SB_unfished_eq),
                            Catch=as.vector(Catch_B_Month),
                            Recruits=as.vector(N_Age_fished[,1,]),
                            SPR=as.vector(SPR),
                            Effort=as.vector(Effort_month),
                            F_mort=as.vector(F_Month),
                            Rec_Devs=as.vector(Rec_Devs[,(maxage+1):ncol(Rec_Devs)]))


  out <- list()
  out$LifeHistory <- LifeHistory
  out$Exploitation <- Exploitation
  out$At_Age_Time_Series <- At_Age_Time_Series
  out$Time_Series <- Time_Series
  out$currentYr <- currentYr
  class(out) <- 'Simulation'
  out
}

generate_Effort <- function(Exploitation, LifeHistory, nsim) {

  if (!Exploitation$Effort_pattern %in% c('Stable', 'Increasing', 'Decreasing')) {
    stop("Effort not valid. Must be one of: c('Stable', 'Increasing', 'Decreasing')")
  }

  # Calculate optimal monthly seasonal fishing pattern
  opt_F_pattern <- calculate_optimal_fishing(LifeHistory, Exploitation,
                                             opt_type=1, utilpow=Exploitation$HARA_power)
  month_opt_Eff <- opt_F_pattern$F_m/mean(opt_F_pattern$F_m)
  nmonths <- Exploitation$nmonths
  Effort_pattern <- Exploitation$Effort_pattern
  Effort_current <- Exploitation$Effort_current
  Effort_cv <- Exploitation$Effort_cv
  Effort_dev <- rlnorm(nmonths*nsim, -0.5*Effort_cv^2, Effort_cv)
  Effort_dev <- matrix(Effort_dev, nrow=nsim, ncol=nmonths)

  # Stable
  if (Effort_pattern=='Stable') {
    Effort_month <- rep(0, nmonths)
    ramp_up <- 1/6 * nmonths
    Effort_month[1:ramp_up] <- seq(0, Effort_current, length.out=ramp_up)
    Effort_month[(ramp_up+1):nmonths] <- Effort_current
    Effort_month <- Effort_month * month_opt_Eff
  }

  # Increasing
  if (Effort_pattern=='Increasing') {
    Effort_month <- rep(0, nmonths)
    ramp_up <- 1/6 * nmonths
    Effort_month[1:ramp_up] <- seq(0, Effort_current*0.4, length.out=ramp_up)
    Effort_month[(ramp_up+1):nmonths] <- seq(Effort_current*0.4, Effort_current, length.out=nmonths-ramp_up)
    Effort_month <- Effort_month * month_opt_Eff
  }
  # Decreasing
  if (Effort_pattern=='Decreasing') {
    Effort_month <- rep(0, nmonths)
    ramp_up <- 1/6 * nmonths
    Effort_month[1:ramp_up] <- seq(0, Effort_current*1.4, length.out=ramp_up)
    Effort_month[(ramp_up+1):nmonths] <- seq(Effort_current*1.4, Effort_current, length.out=nmonths-ramp_up)
    Effort_month <- Effort_month * month_opt_Eff
  }

  Effort_month <- t(replicate(nsim, Effort_month))
  Effort_month * Effort_dev
}

#' Generate Sampled Data from a Simulated Fishery
#'
#' @param Simulation An  object of class `Simulation` generated by `Simulate()`
#' @param Sampling A `list` object with the parameters for generating data
#' @param seed The seed for the random number generator
#' @param silent Hide messages?
#' @export
Generate_Data <- function(Simulation=NULL, Sampling=NULL, seed=101, silent=FALSE) {

  set.seed(seed)

  if (!inherits(Simulation, 'Simulation'))
    stop('`Simulation` must be class `Simulation`. Use `Simulate()`')

  if (is.null(Sampling)) {
    if (!silent) message('Using Example Sampling Parameters from `SLAM::Sampling`')
    Sampling <- SLAM::Sampling
  }
  if (inherits(Sampling, 'character')) {
    Sampling <-  Import_Sampling(Sampling)
  }

  # Import Data Sampling Parameters
  n_recent_months <- Sampling$n_recent_months
  Rel_Sample_Month <- Sampling$Rel_Sample_Month
  CPUE_CV <- Sampling$CPUE_CV
  Catch_CV <- Sampling$Catch_CV
  Effort_CV <- Sampling$Effort_CV

  Weight_Bin_Max <- Sampling$Weight_Bin_Max
  Weight_Bin_Width <- Sampling$Weight_Bin_Width
  Weight_Bins <- seq(0, Weight_Bin_Max, Weight_Bin_Width)
  nBins <- length(Weight_Bins)-1
  Weight_Mids <- seq(0.5*Weight_Bin_Width, by=Weight_Bin_Width, length.out=nBins)

  Sampling$Weight_Bins <- Weight_Bins
  Sampling$Weight_Mids <- Weight_Mids

  CAW_Annual_Sample_Size <- Sampling$CAW_Annual_Sample_Size
  CAW_Annual_ESS <- Sampling$CAW_Annual_ESS

  if (is.null(Sampling$CAA_Annual_Sample_Size))
    Sampling$CAA_Annual_Sample_Size <- Sampling$CAW_Annual_Sample_Size
  if (is.null(Sampling$CAA_Annual_ESS))
    Sampling$CAA_Annual_ESS <- Sampling$CAW_Annual_ESS
  CAA_Annual_Sample_Size <- Sampling$CAA_Annual_Sample_Size
  CAA_Annual_ESS <- Sampling$CAA_Annual_ESS

  # Import Life-History Parameters
  Stock_Name <- Simulation$LifeHistory$Stock_Name
  Species <- Simulation$LifeHistory$Species
  Ages <- Simulation$LifeHistory$Ages
  maxage <- max(Ages)
  nAge <- length(Ages)
  Weight_Age_Mean <- Simulation$LifeHistory$Weight_Age_Mean
  Weight_Age_SD <- Simulation$LifeHistory$Weight_Age_SD
  M_at_Age <- Simulation$LifeHistory$M_at_Age
  Maturity_at_Age <- Simulation$LifeHistory$Maturity_at_Age
  Post_Spawning_Mortality <- Simulation$LifeHistory$Post_Spawning_Mortality

  nsim <- length(unique(Simulation$At_Age_Time_Series$Sim))
  nts  <- length(unique(Simulation$At_Age_Time_Series$Month_ind))

  # Generate Data (for all years)
  # Calculate Age-Weight Key (assume log-normal variability in weight-at-age)
  AWK <- matrix(0, nrow=nAge, ncol=nBins)
  mu <- log(Weight_Age_Mean) -0.5*Weight_Age_SD^2
  AWK[,1] <- plnorm(Weight_Bins[2], mu, Weight_Age_SD) # probability of weight-at-age
  for (i in 2:(nBins-1)) {
    AWK[,i] <- plnorm(Weight_Bins[i+1], mu, Weight_Age_SD) -
      plnorm(Weight_Bins[i], mu, Weight_Age_SD)
  }
  AWK[,nBins] <- 1 - plnorm(Weight_Bins[nBins], mu, Weight_Age_SD)

  # Catch-at-Weight - expected
  Catch_Age <- array(NA, dim=c(nsim, nAge, nts))
  Catch_Age[] <- Simulation$At_Age_Time_Series$Catch_n
  CAW_exp <- array(0, dim=c(nsim, nBins, nts))
  for (ts in 1:nts) {
    CAW_exp[,,ts] <- t(sapply(1:nsim, function(x)
      apply(AWK * Catch_Age[x,,ts], 2, sum)))
  }

  # Catch (absolute)
  Catch_B_Month <- array(NA, dim=c(nsim, nts))
  Catch_B_Month[] <- Simulation$Time_Series$Catch
  Catch_Sample <- Catch_B_Month * exp(rnorm(nts*nsim, -0.5*Catch_CV^2, Catch_CV))

  # CPUE
  B_Month <- array(NA, dim=c(nsim, nts))
  B_Month[] <- Simulation$Time_Series$B_fished
  CPUE_Sample <- B_Month/apply(B_Month, 1, mean) * exp(rnorm(nts*nsim, -0.5*CPUE_CV^2, CPUE_CV))

  # Effort (Index)
  Effort_Month <- array(NA, dim=c(nsim, nts))
  Effort_Month[] <- Simulation$Time_Series$Effort
  Effort_Sample <- Effort_Month/apply(Effort_Month, 1, mean) * exp(rnorm(nts*nsim, -0.5*Effort_CV^2, Effort_CV))

  # Generate Monthly Samples for n_recent_months
  sample_ts <- (nts-n_recent_months+1):nts
  n_sample_ts <- length(sample_ts)
  Catch_Sample <- Catch_Sample[,sample_ts, drop=FALSE]
  CPUE_Sample <- CPUE_Sample[,sample_ts, drop=FALSE]
  CPUE_Sample <- CPUE_Sample/apply(CPUE_Sample, 1, mean)
  Effort_Sample <- Effort_Sample[,sample_ts, drop=FALSE]
  Effort_Sample <- Effort_Sample/apply(Effort_Sample, 1, mean)

  CAW_Sample <- array(0, dim=c(nsim, nBins, n_sample_ts))
  Rel_Sample_Month <- Rel_Sample_Month/sum(Rel_Sample_Month)

  for (i in seq_along(sample_ts)) {
    ts <- sample_ts[i]
    month <- ts %%12
    if (month==0) month <- 12
    ts_sample_size <- CAW_Annual_Sample_Size * Rel_Sample_Month[month]
    ts_ess <- CAW_Annual_ESS * Rel_Sample_Month[month]

    CAW_Sample[,,i] <- t(sapply(1:nsim, function(x)
      if (sum(CAW_exp[x,,ts])>0) {
        val <-  ts_sample_size * (rmultinom(1, size=ts_ess, prob=CAW_exp[x,,ts]))/ts_ess
      } else {
        val <- rep(0, nBins)
      }
    ))
  }

  CAA_Sample <- array(0, dim=c(nsim, nAge, n_sample_ts))
  for (i in seq_along(sample_ts)) {
    ts <- sample_ts[i]
    month <- ts %%12
    if (month==0) month <- 12
    ts_sample_size <- CAA_Annual_Sample_Size * Rel_Sample_Month[month]
    ts_ess <- CAA_Annual_ESS * Rel_Sample_Month[month]

    CAA_Sample[,,i] <- t(sapply(1:nsim, function(x)
      if (sum(Catch_Age[x,,ts])>0) {
        val <-  ts_sample_size * (rmultinom(1, size=ts_ess, prob=Catch_Age[x,,ts]))/ts_ess
      } else {
        val <- rep(0, nAge)
      }
    ))
  }

  # Data
  month_ind <- sample_ts %%12
  month_ind[month_ind==0] <- 12

  Months <- month.abb[month_ind]
  n_recent_years <- max(table(Months))
  currentYr <- Simulation$currentYr
  Years <- (currentYr-n_recent_years+1):currentYr

  Data_TS_DF <- data.frame(Sim=1:nsim,
                           Year=rep(Years, each=nsim),
                           Month=rep(Months, each=nsim),
                           Month_ind=rep(sample_ts, each=nsim),
                           Catch=as.vector(Catch_Sample),
                           CPUE=as.vector(CPUE_Sample),
                           Effort=as.vector(Effort_Sample))

  Data_CAW_DF <- data.frame(Sim=1:nsim,
                            Weight=rep(Weight_Mids, each=nsim),
                            Year=rep(Years, each=nsim*nBins),
                            Month=rep(Months, each=nsim*nBins),
                            Month_ind=rep(sample_ts, each=nsim*nBins),
                            Count=as.vector(CAW_Sample))

  Data_CAA_DF <- data.frame(Sim=1:nsim,
                            Age=rep(Ages, each=nsim),
                            Year=rep(Years, each=nsim*nAge),
                            Month=rep(Months, each=nsim*nAge),
                            Month_ind=rep(sample_ts, each=nsim*nAge),
                            Count=as.vector(CAA_Sample))



  out <- list()
  out$Simulation <- Simulation
  out$Sampling <- Sampling
  out$Data <- list(TS=Data_TS_DF,
                   CAW=Data_CAW_DF,
                   CAA=Data_CAA_DF,
                   Weight_Bins=Weight_Bins,
                   Weight_Mids=Weight_Mids)

  class(out) <- 'Simulated'
  out
}






