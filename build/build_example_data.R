Simulation <- Simulate(LifeHistory, Exploitation, nsim=2)

Write_Data2CSV <- function(Simulation, csvfile, n_years=3, sim=1, Sampling=NULL, dev=F) {

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
  df_out <- data.frame(Year=Sampled_Data$Year, Month=Sampled_Data$Month, t(CAW))
  colnames(df_out) <- c('Year', 'Month', Sampled_Data$Weight_Mids)
  df_CAW <- df_out

  # Effort Data
  df <- Data$Data$TS %>% filter(Sim==sim) %>% select(Year, Month,
                                                     Effort_Mean=Effort)
  df$Effort_SD <- Data$Sampling$Effort_CV
  df_Effort <- df %>% select(Year, Month, Effort_Mean, Effort_SD)
  df_Effort$Effort_Mean <- round(df_Effort$Effort_Mean,2) *100

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

  # Effort Data
  cat('\n', file=csvfile, append = TRUE)
  cat('Effort Data\n', file=csvfile,  append = TRUE)
  cat(paste(names(df_Effort), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_Effort)) {
    cat(paste0(paste(df_Effort[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

}


Write_Data2CSV(Simulation, csvfile='inst/Data_Example.csv', dev=T)


