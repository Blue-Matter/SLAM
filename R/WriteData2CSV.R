Write_Data2CSV <- function(Data, csvfile, dir=NULL, overwrite=FALSE) {

  if (is.null(dir)) dir <- getwd()
  if (file.exists(file.path(dir, csvfile)) & !overwrite) {
    stop(file.path(dir, csvfile), ' already exists. Use `overwrite=TRUE` to overwrite')
  } else {
    message('Creating: ', file.path(dir, csvfile))
    file.create(file.path(dir, csvfile))
  }

  csvfile <- file.path(dir, csvfile)

  df_meta <- Data$Metadata

  # Meta-data
  cat('Field, Value\n', file=csvfile)
  for (i in 1:nrow(df_meta)) {
    cat(paste0(paste(df_meta[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  # At-Age Schedules
  n_age <- length(Data$Ages)

  Ages <- Data$Ages
  Weight_Age_Mean <- Data$Weight_Age_Mean
  Weight_Age_SD <- Data$Weight_Age_SD
  Maturity_at_Age <- Data$Maturity_at_Age
  M_at_Age <- Data$M_at_Age
  Post_Spawning_Mortality <- Data$Post_Spawning_Mortality

  df_age <- data.frame(Ages, Weight_Age_Mean, Weight_Age_SD, Maturity_at_Age, M_at_Age, Post_Spawning_Mortality)

  cat('\n', file=csvfile, append = TRUE)
  cat('Age Schedules\n', file=csvfile,  append = TRUE)
  cat(paste(names(df_age), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_age)) {
    cat(paste0(paste(df_age[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  CAW <- Data$CAW
  dd <- dim(CAW)
  Catch_n <- matrix(Data$Catch_n, nrow=nrow(Data$CAW), ncol=ncol(Data$CAW))
  df_out <- data.frame(Year=Data$Year, Month=Data$Month, Data$CAW*Catch_n )
  colnames(df_out) <- c('Year', 'Month', Data$Weight_Mids)
  df_CAW <- df_out

  cat('\n', file=csvfile, append = TRUE)
  cat('CAW Binned Data\n', file=csvfile,  append = TRUE)
  cat(paste(names(df_CAW), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df_CAW)) {
    cat(paste0(paste(df_CAW[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }

  # Index Data


  cat('\n', file=csvfile, append = TRUE)
  cat('Index Data\n', file=csvfile,  append = TRUE)

  df_temp <- data.frame(Year=Data$Year, Month=Data$Month, Mean=NA, SD=NA, Type=NA)
  df_effort <- df_temp
  df_catch <- df_temp
  df_index <- df_temp

  # effort
  if (!all(is.na(Data$Effort_Mean))) {
    df_effort$Type <- 'Effort'
    df_effort$Mean <- Data$Effort_Mean
    df_effort$SD <- Data$Effort_SD
  }

  # catch
  if (!all(is.na(Data$Catch_Mean))) {

    df_catch$Type <- 'Catch'
    df_catch$Mean <- Data$Catch_Mean
    df_catch$SD <- Data$Catch_SD
  }

  # biomass
  if (!all(is.na(Data$Index_Mean))) {
    df_index$Type <- 'Biomass'
    df_index$Mean <- Data$Index_Mean
    df_index$SD <- Data$Index_SD
  }

  df <- bind_rows(df_effort, df_catch, df_index) %>% filter(is.na(Type)==FALSE)


  cat(paste(names(df), collapse=', '), file=csvfile,  append=TRUE)
  cat('\n', file=csvfile, append = TRUE)
  for (i in 1:nrow(df)) {
    cat(paste0(paste(df[i,], collapse=', '), '\n'), file=csvfile, append=TRUE)
  }



  # Parameters
  sigmaR <- 0.4
  sigmaF <- 0.4
  sigmaR0 <- 0.3
  Steepness <- 0.85
  if (length(Data$sigmaR)>0) sigmaR <- Data$sigmaR
  if (length(Data$sigmaF)>0) sigmaF <- Data$sigmaF
  if (length(Data$sigmaR0)>0) sigmaR0 <- Data$sigmaR0
  if (length(Data$Steepness)>0) Steepness <- Data$Steepness


  cat('\n', file=csvfile, append = TRUE)
  cat('Parameters, Value, Description\n', file=csvfile,  append = TRUE)
  cat(paste('Steepness', Steepness, 'Assumed steepness of Beverton-Holt Stock-Recruit Relationship', sep=', '), '\n', file=csvfile,  append = TRUE)
  cat(paste('sigmaR', sigmaR, 'Assumed standard deviation of log-normal recruitment deviations', sep=', '), '\n', file=csvfile,  append = TRUE)
  cat(paste('sigmaF', sigmaF, 'Standard deviation for penalty for random walk in monthly F (log-space)', sep=', '), '\n', file=csvfile,  append = TRUE)
  cat(paste('sigmaR0', sigmaR0, 'Standard deviation for penalty for random walk in monthly recruitment pattern (log-space)', sep=', '), '\n', file=csvfile,  append = TRUE)

}
