Simulation <- Simulate(LifeHistory, Exploitation, nsim=2)

Write_Data2CSV.Simulation <- function(Simulation, csvfile, n_years=3, sim=1, Sampling=NULL, dev=T) {

  if (dev) {
    dir <- 'inst'
  } else {
    dir <- system.file(package='SLAM')
  }


  # Fishery Info
  df_meta <- data.frame(Field=c('Fishery', 'Location', 'Species', 'Common_Name', 'Author', 'Date'),
                   Value=c('Example Octopus Fishery', 'Indonesia', 'Octopus cyanea', 'Day octopus',
                           'T. Author (t.author@email.com)',
                           'March 2023'))


  # At-Age Schedules
  n_age <- length(Simulation$LifeHistory$Weight_Age_Mean)

  Ages <- Simulation$LifeHistory$Ages
  Weight_Age_Mean <- Simulation$LifeHistory$Weight_Age_Mean
  Weight_Age_SD <- Simulation$LifeHistory$Weight_Age_SD
  Maturity_at_Age <- Simulation$LifeHistory$Maturity_at_Age
  M_at_Age <- Simulation$LifeHistory$M_at_Age
  Post_Spawning_Mortality <- Simulation$LifeHistory$Post_Spawning_Mortality

  df_age <- data.frame(Ages, Weight_Age_Mean, Weight_Age_SD, Maturity_at_Age, M_at_Age, Post_Spawning_Mortality)

  # CAW Binned Data
  if (is.null(Sampling)) {
    Sampling <- SLAM::Sampling
    Sampling$n_recent_months <- 12*n_years
  }
  Data <- Generate_Data(Simulation, Sampling)

  Sampled_Data <- Import(Data)
  CAW <- Sampled_Data$CAW
  dd <- dim(CAW)
  nBins <- dd[1]
  mat <- matrix(0, length(Sampled_Data$Year), nBins)
  df_out <- data.frame(Year=Sampled_Data$Year, Month=Sampled_Data$Month, CAW)
  colnames(df_out) <- c('Year', 'Month', Sampled_Data$Weight_Mids)
  df_CAW <- df_out
  df_CAW$Month <- match(df_CAW$Month, month.abb)

  # Effort Data
  df <- Data$Data$TS %>% filter(Sim==sim) %>% select(Year, Month,
                                                     Effort_Mean=Effort)
  df$Effort_SD <- Data$Sampling$Effort_CV
  df_Effort <- df %>% select(Year, Month, Effort_Mean, Effort_SD)
  df_Effort$Effort_Mean <- round(df_Effort$Effort_Mean,2) *100
  df_Effort$Month <- match(df_Effort$Month, month.abb)

  # Index of Abundance
  df_Abund <- Data$Data$TS %>% filter(Sim==sim) %>% select(Year, Month,
                                                        Mean=CPUE)
  df_Abund$SD <- Data$Sampling$CPUE_CV
  df_Abund <- df_Abund %>% select(Year, Month, Mean, SD)
  df_Abund$Mean <- round(df_Abund$Mean,2) *100
  df_Abund$Month <- match(df_Abund$Month, month.abb)

  # write csv

  # Meta-data
  cat('Field, Value\n', file=csvfile)
  for (i in 1:nrow(df_meta)) {
    cat(paste0(paste(df_meta[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  # At-Age Schedules
  cat('\n', file=csvfile, append = TRUE)
  cat('Age Schedules\n', file=csvfile,  append = TRUE)
  cat(paste(names(df_age), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_age)) {
    cat(paste0(paste(df_age[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  # CAW Data
  cat('\n', file=csvfile, append = TRUE)
  cat('CAW Binned Data\n', file=csvfile,  append = TRUE)
  cat(paste(names(df_CAW), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_CAW)) {
    cat(paste0(paste(df_CAW[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  # Index Data

  # effort
  cat('\n', file=csvfile, append = TRUE)
  cat('Index Data\n', file=csvfile,  append = TRUE)
  names(df_Effort)[3:4] <- c('Mean', 'SD')
  df_Effort$Type <- 'Effort'
  cat(paste(names(df_Effort), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_Effort)) {
    cat(paste0(paste(df_Effort[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  # abundance
  df_Abund$Type <- 'Biomass'
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_Abund)) {
    cat(paste0(paste(df_Abund[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }


  # Parameters
  sigmaF <- 0.4
  sigmaR0 <- 0.3

  cat('\n', file=csvfile, append = TRUE)
  cat('Parameters, Value, Description\n', file=csvfile,  append = TRUE)
  cat(paste('Steepness', Simulation$LifeHistory$steepness, 'Assumed steepness of Beverton-Holt Stock-Recruit Relationship', sep=', '), '\n', file=csvfile,  append = TRUE)
  cat(paste('sigmaR', Simulation$LifeHistory$sigmaR, 'Assumed standard deviation of log-normal recruitment deviations', sep=', '), '\n', file=csvfile,  append = TRUE)
  cat(paste('sigmaF', sigmaF, 'Standard deviation for penalty for random walk in monthly F (log-space)', sep=', '), '\n', file=csvfile,  append = TRUE)
  cat(paste('sigmaR0', sigmaR0, 'Standard deviation for penalty for random walk in monthly recruitment pattern (log-space)', sep=', '), '\n', file=csvfile,  append = TRUE)

}


Write_Data2CSV.Simulation(Simulation, csvfile='inst/Data_Example.csv', dev=T)


