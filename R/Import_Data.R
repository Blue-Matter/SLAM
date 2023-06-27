
#' Print the file path for example data sets included in SLAM
#'
#' @return The file paths for the example data
#' @export
#'
#' @examples
#' Example_Data()
Example_Data <- function() {
  dir <- system.file(package='SLAM')
  fls <- file.path(dir, list.files(dir, pattern='csv'))
  ind <- which(basename(fls) == "Simulation_Parameters.csv")
  fls <- fls[-ind]
  fls
}

#' Import and organize data for assessment model
#'
#' @param csvfile The full file path to an CSV file structured for SLAM data
#' @param silent Logical. Hide the messages?
#' @param ... Additional parameters (not currently used)
#'
#' @return A `Data` object for the `Assess` function
#' @export
Import <- function(csvfile,
                   silent=FALSE,
                   ...) {
  if (inherits(csvfile, 'Simulated')) {
    UseMethod('Import')
  } else {

    tempin <- strsplit(readLines(csvfile), ",")
    somechar <- function(x) sum(nchar(x) > 0)
    Ncol <- max(unlist(lapply(tempin, somechar)))
    col.names <- paste0("V", 1:Ncol)

    datasheet <- read.csv(csvfile, header = FALSE,
                          colClasses = "character", col.names=col.names,
                          stringsAsFactors = FALSE, row.names = NULL)



    data <- Make_Data(TRUE)
    # meta-data
    for (i in 1:nrow(data$Metadata)) {
      ind <- which(data$Metadata$Field[i]==datasheet$V1)
      if (length(ind)<1)
        warning(data$Metadata$Field[i], ' not found in ', basename(csvfile))
      data$Metadata$Value[i] <- datasheet[ind,2]
    }

    # age
    ind <- which(datasheet$V1 =='Ages')
    age_rows <- suppressWarnings(as.numeric(datasheet$V1[(ind+1):length(datasheet$V1)]))
    age_rows <- 1:(min(which(is.na(age_rows)))-1)
    row_ind <- (ind+age_rows)
    age_df <- datasheet[row_ind,]
    df <- data.frame(age_df)
    colnames(df) <- trimws(datasheet[ind,])
    df <- df[,which(nchar(df[1,])>0)]

    for (i in 1:ncol(df)) {
      data[[colnames(df)[i]]] <-as.numeric(df[[colnames(df)[i]]])

    }

    # CAW
    CAW_ind <- which(datasheet$V1 =='Year')[1]
    year_rows <- suppressWarnings(as.numeric(datasheet$V1[(CAW_ind+1):length(datasheet$V1)]))
    year_rows <- 1:(min(which(is.na(year_rows)))-1)
    row_ind <- (CAW_ind+year_rows)
    CAW_df <- datasheet[row_ind,3:ncol(datasheet)]
    CAW_df <- apply(CAW_df, 2, as.numeric)
    n_months <- nrow(CAW_df)

    data$Weight_Mids <- Weight_Mids <- as.numeric(datasheet[CAW_ind,3:ncol(datasheet)])
    by <-Weight_Mids[2]- Weight_Mids[1]
    Weight_Bins <- seq(Weight_Mids[1]-0.5*by, by=by, length.out=length(Weight_Mids)+1)
    data$Weight_Bins <- Weight_Bins

    CAW_n <- apply(CAW_df, 1, sum)
    data$CAW <- CAW_df/CAW_n
    data$CAW_ESS <- CAW_n
    data$Year <- as.numeric(datasheet$V1[row_ind])
    data$Month <- as.numeric(datasheet$V2[row_ind])

    # effort
    Effort_ind <- which(datasheet$V1 =='Year')[2]
    year_rows <- (Effort_ind+1):length(datasheet$V1)

    effort_year <- as.numeric(datasheet$V1[year_rows])
    effort_month <- as.numeric(datasheet$V2[year_rows])
    effort_mean <- as.numeric(datasheet$V3[year_rows])
    effort_sd <- as.numeric(datasheet$V4[year_rows])

    Year_Month <- paste(data$Year, data$Month, sep='_')
    Effort_Year_Month <- paste(effort_year,effort_month, sep='_')
    data$Effort_Mean <- rep(NA, n_months)
    data$Effort_SD <- rep(NA, n_months)
    ind <- match(Effort_Year_Month, Year_Month)
    data$Effort_Mean[ind] <- effort_mean
    data$Effort_SD[ind] <- effort_sd

    # index - not currently used
    data$Index_Mean <- rep(NA, n_months)
    data$Index_SD <- rep(NA, n_months)
    return(data)
  }
}
#
# import_at_age <- function(XLData, data) {
#   At_Age_Data <- XLData$`At-Age-Schedules`
#   if (is.null(At_Age_Data)) {
#     stop('Could not import At-Age-Schedules from ', xlfile, ' . Is the worksheet named `At-Age-Schedules`?')
#   }
#   Names <- At_Age_Data$At.Age.Schedules %>% as.vector()
#   Ncol <- ncol(At_Age_Data)
#   n_age <- Ncol-1
#
#   for (i in seq_along(Names)) {
#     dd <- At_Age_Data[i,2:(Ncol)]
#     dd <- as.matrix(dd)
#     colnames(dd) <- NULL
#     data[[Names[i]]] <- as.vector(dd)
#   }
#
#   data
# }
#
# import_caw_data <- function(XLData, data, BinWidth=NULL, BinMax=NULL) {
#   CAW_Data <- XLData$`CAW-Data`
#   if (is.null(CAW_Data)) {
#     stop('Could not import CAW-Data from ', xlfile, ' . Is the worksheet named `CAW-Data`?')
#   }
#   cnames <- colnames(CAW_Data)
#   if (length(cnames)>3) {
#     # binned data
#     BinMids <- as.numeric(cnames[3:length(cnames)])
#     by <- BinMids[2]-BinMids[1]
#     Bins <- seq(BinMids[1]-0.5*by, by=by, length.out=length(BinMids)+1)
#     nMonths <- nrow(CAW_Data)
#     nBins <- length(BinMids)
#     data$CAW <- as.matrix(CAW_Data[,3:length(cnames)])
#     data$Weight_Bins <- Bins
#     data$Weight_Mids <- BinMids
#     data$CAW <- t(data$CAW)
#     data$CAW_ESS <- apply(data$CAW, 2, sum)
#     data$Year <- XLData$`CAW-Data`$Year
#     data$Month <- XLData$`CAW-Data`$Month
#
#   } else {
#     message('Raw CAW data detected')
#     if (is.null(BinWidth)) {
#       message('Argument `BinWidth` not set. Using default value of 0.1')
#       BinWidth <- 0.1
#     }
#
#     if (!is.null(BinMax)) {
#       maxBin <- BinMax
#       XLData$`CAW-Data`$Weight[XLData$`CAW-Data`$Weight>maxBin] <- maxBin
#       message('Setting maximum bin width to: ', maxBin)
#     } else {
#       maxBin <- ceiling(max(XLData$`CAW-Data`$Weight))
#       message('Setting maximum bin width to maximum observed weight: ', maxBin)
#     }
#     BinMids <- seq(0.5*BinWidth, by=BinWidth, maxBin)
#     nBins <- length(BinMids)
#     Bins <- seq(BinMids[1]-0.5*BinWidth, by=BinWidth, length.out=length(BinMids)+1)
#     Years <- unique(XLData$`CAW-Data`$Year)
#     CAW <- XLData$`CAW-Data`
#     Year_Month <- CAW %>% distinct(Year, Month)
#     nMonths <- Year_Month %>% nrow()
#     data$CAW <- matrix(0, nBins, nMonths)
#     for (i in 1:nMonths) {
#       ts <- Year_Month[i,]
#       caw_yr <- CAW %>% filter(Year %in% ts$Year, Month %in% ts$Month)
#       data$CAW[, i] <- as.numeric(table(cut(caw_yr$Weight, Bins)))
#     }
#     data$CAW_ESS <- apply(data$CAW, 2, sum)
#     Year_Month <- data.frame(Year=XLData$`CAW-Data`$Year, Month=XLData$`CAW-Data`$Month) %>%
#       distinct(Year, Month)
#     data$Year <- Year_Month$Year
#     data$Month <- Year_Month$Month
#   }
#   data
# }
#
# import_ts_data <- function(XLData, data) {
#   TS_data <- XLData$`Time-Series-Data`
#   TS_Year_Month <- data.frame(Year=TS_data$Year, Month=TS_data$Month)
#   CAW_Year_Month <- data.frame(Year=data$Year, Month=data$Month)
#
#   # match Year and Months
#   match <- FALSE
#   if (all(dim(TS_Year_Month) ==dim(CAW_Year_Month))) {
#     match <- all(TS_Year_Month == CAW_Year_Month)
#   }
#   if (!match)
#     stop('`CAW-Data` and `Time-Series-Data` must both contain all Years and Months. Use NA for missing values.')
#
#   data$Effort_Mean <- TS_data$Effort_Mean
#   data$Effort_SD <- TS_data$Effort_SD
#   data$Index_Mean <- TS_data$Index_Mean
#   data$Index_SD <- TS_data$Index_SD
#   data
# }









#' Import Simulated Data
#'
#' @param Sampled_Data An object of class `Simulated` generated by `Generate_Data`
#' @param sim Simulation number
#' @param Fit_Effort Logical. Fit effort data?
#' @param Fit_CPUE Logical. Fit CPUE data?
#' @param Fit_CAW Logical. Fit catch-at-weight data?
#' @param Fit_CAA Logical. Fit catch-at-age data?
#' @param use_Frwpen Logical. Use the random walk penalty for F?
#' @param use_R0rwpen Logical. Use the random walk penalty for seasonal recruitment?
#' @param Data_types Optional character vector specifying data types to fit to.
#'
#' @return A `Data` object
#' @export
Import.Simulated <- function(Sampled_Data=NULL,
                             sim=1,
                             Fit_Effort=1,
                             Fit_CPUE=1,
                             Fit_CAW=1,
                             Fit_CAA=1,
                             use_Frwpen=1,
                             use_R0rwpen=1,
                             CAW_ESS=200,
                             Effort_CV=0.2,
                             Index_CV=0.3,
                             Data_types=NULL) {

  data <- list()
  # Assumed life-history parameters
  data$Weight_Age_Mean <- Sampled_Data$Simulation$LifeHistory$Weight_Age_Mean
  data$Weight_Age_SD <- Sampled_Data$Simulation$LifeHistory$Weight_Age_SD
  data$Maturity_at_Age <- Sampled_Data$Simulation$LifeHistory$Maturity_at_Age
  data$M_at_Age <- Sampled_Data$Simulation$LifeHistory$M_at_Age
  data$Post_Spawning_Mortality <- Sampled_Data$Simulation$LifeHistory$Post_Spawning_Mortality
  data$h <- Sampled_Data$Simulation$LifeHistory$steepness
  nAges <- length(data$Weight_Age_Mean)

  # CAW Data
  data$Weight_Bins <- Sampled_Data$Data$Weight_Bins
  data$Weight_Mids <- Sampled_Data$Data$Weight_Mids

  CAW_DF <- Sampled_Data$Data$CAW %>% filter(Sim==sim) %>% select(Month_ind, Weight, Count)
  nBins <- length(unique(CAW_DF$Weight))
  nMonths <- length(unique(CAW_DF$Month_ind))

  CAW <- matrix(CAW_DF$Count, nrow=nBins, nMonths)
  data$CAW <- CAW
  CAW_Monthly_ESS <- Sampled_Data$Sampling$CAW_Annual_ESS/12
  CAW_Monthly_ESS[CAW_Monthly_ESS>CAW_ESS] <- CAW_ESS
  data$CAW_ESS <- rep(CAW_Monthly_ESS, nMonths)

  # CAA Data
  CAA_DF <- Sampled_Data$Data$CAA %>% filter(Sim==sim) %>% select(Month_ind, Age, Count)
  CAA <- matrix(CAA_DF$Count, nrow=nAges, nMonths)
  data$CAA <- CAA
  CAA_Monthly_ESS <- Sampled_Data$Sampling$CAA_Annual_ESS/12
  data$CAA_ESS <- rep(CAA_Monthly_ESS, nMonths)

  # Effort Index
  Effort_DF <- Sampled_Data$Data$TS %>% filter(Sim==sim) %>% select(Effort)
  data$Effort_Mean <- Effort_DF$Effort
  data$Effort_SD <- rep(Effort_CV, nMonths)

  # CPUE Index
  CPUE_DF <- Sampled_Data$Data$TS %>% filter(Sim==sim) %>% select(CPUE)
  data$Index_Mean <- CPUE_DF$CPUE
  data$Index_SD <- rep(Index_CV, nMonths)

  # Year and Month info
  tt <- Sampled_Data$Data$TS %>% filter(Sim==sim)
  data$Month_ind <- tt$Month_ind
  data$n_month <- length(data$Month_ind)
  data$Year <- tt$Year
  data$Month <- tt$Month

  # Options
  data$Fit_Effort <- Fit_Effort
  data$Fit_Index <- Fit_CPUE
  data$Fit_CAW <- Fit_CAW
  data$Fit_CAA <- Fit_CAA

  # Penalties
  data$use_Frwpen <- use_Frwpen
  data$use_R0rwpen <- use_R0rwpen

  data$model <- 'SLAM'

  if (!is.null(Data_types)) {
    data <- set_data_types(data, Data_types=Data_types)
  }
  class(data) <- 'Data'
  data
}


