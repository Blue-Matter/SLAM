
#' Process Case Study Data
#'
#' @param Sites Names of the village
#' @param Case_Study_Data A data.frame
#'
#' @return A list
#' @export
Process_Data <- function(Sites, Case_Study_Data,
                         maxWeightBin=4, weight_bin=0.25) {

  agg_villages <- strsplit(Sites, ',')[[1]] %>% trimws

  village_data <- Case_Study_Data %>% filter(Village %in% agg_villages)

  # ---- Format Dates ----
  village_data$Year <- village_data$year
  village_data$Date <-  lubridate::date(village_data$Fishing.Date)
  village_data$Date2 <-  lubridate::my(paste(village_data$Month, village_data$Year, sep="-"))
  village_data$Month_Name <- month.abb[village_data$Month]

  # ---- Weight Composition Data -----
  weight_data <- village_data %>% select(Date2,
                                         Sex,
                                         Length=Mantle.length..Cm.,
                                         Weight=Octopus.individual.weight..kg.) %>%
    mutate(Weight=as.numeric(Weight)) %>%
    filter(is.na(Length)==FALSE, is.na(Weight)==FALSE,
           Weight>0, Length>0)


  by <- weight_bin
  Weight_Bins <- seq(0, to=maxWeightBin, by=by)
  Weight_Mids <- seq(0.5*by, by=by, length.out=length(Weight_Bins)-1)

  weight_data$Weight[weight_data$Weight>=maxWeightBin] <- maxWeightBin
  WeightData_Binned <- weight_data %>% group_by(Date2) %>%
    reframe(bin = cut(Weight, breaks = Weight_Bins))

  WeightData_Binned <- WeightData_Binned %>% group_by(Date=Date2) %>%
    reframe(Count=as.numeric(table(bin)), Bin=Weight_Mids)

  WeightData_Binned$Year <- lubridate::year(WeightData_Binned$Date)
  WeightData_Binned$Month <- lubridate::month(WeightData_Binned$Date)


  # Continuous Dates
  Years <- range(WeightData_Binned$Year)
  Years <- seq(Years[1], Years[2], by=1)

  Date_DF <- data.frame(Year=Years, Month=rep(1:12, each=length(Years))) %>%
    arrange(Year) %>%
    mutate(Date=lubridate::my(paste(Month, Year, sep="-")))

  WeightData_Binned <- left_join(Date_DF, WeightData_Binned, by = join_by(Year, Month, Date))
  WeightData_Binned$Sites <- Sites

  keep <- range(which(!is.na(WeightData_Binned$Count)))
  keep_rows <- keep[1]:keep[2]
  WeightData_Binned <- WeightData_Binned[keep_rows,]


  # --- Effort Data ----
  village_data$Duration <- suppressWarnings(as.numeric(village_data$Fishing.Duration..hour.))

  if (!all(is.na(unique(village_data$Duration )))) {

    Effort_DF <- village_data %>% group_by(Date=Date2) %>%
      summarize(n_fisher=length(unique(Fisher)),
                mean_effort=mean(Duration, na.rm=TRUE),
                sd_effort=sd(log(Duration), na.rm=TRUE),
                effort=mean_effort*n_fisher) %>%
      ungroup() %>%
      mutate(rel_effort=effort/mean(effort, na.rm=TRUE))

    Effort_DF$sd_effort[Effort_DF$sd_effort<0.1] <- mean(Effort_DF$sd_effort[!Effort_DF$sd_effort==0], na.rm=T)

    Effort_DF$Year <- lubridate::year(Effort_DF$Date)
    Effort_DF$Month <- lubridate::month(Effort_DF$Date)

    Effort_DF <- left_join(Date_DF, Effort_DF, by = join_by(Year, Month, Date))
    keep <- range(which(!is.na(Effort_DF$rel_effort)))
    keep_rows <- keep[1]:keep[2]
    Effort_DF <- Effort_DF[keep_rows,]
    do_CPUE <- TRUE
    if (nrow(Effort_DF)<2)
      do_CPUE <- FALSE
  } else {
    Effort_DF <- NA
    do_CPUE <- FALSE
  }



  # ---- CPUE Data ----
  # Mean CPUE by Fisher and Month
  CPUE_DF <- village_data %>% group_by(Date=Date2, Fisher) %>%
    filter(is.na(Duration)==FALSE,
           is.na(Octopus.individual.weight..kg.)==FALSE) %>%
    summarize(Catch=sum(Octopus.individual.weight..kg., na.rm=TRUE),
              Effort=sum(Duration, na.rm=TRUE)) %>%
    mutate(CPUE=Catch/Effort, lCPUE=log(CPUE)) %>%
    ungroup()

  if(nrow(CPUE_DF)<1 | !do_CPUE) {
    CPUE_DF <- NA
  } else {
    CPUE_DF$Year <- lubridate::year(CPUE_DF$Date)
    CPUE_DF$Month <- lubridate::month(CPUE_DF$Date)
    CPUE_DF$Year_Month <- paste(CPUE_DF$Month, CPUE_DF$Year, sep="_")

    mod <- lm(lCPUE~Fisher+Year_Month, data=CPUE_DF)

    Estimated_CPUE <- exp(predict(mod)+sd(summary(mod)$residuals)^2/2)
    Est_CPUE975 <- exp(predict(mod)+predict(mod, se.fit=T)$se*1.96+sd(summary(mod)$residuals)^2/2)
    Est_CPUE025 <- exp(predict(mod)-predict(mod, se.fit=T)$se*1.96+sd(summary(mod)$residuals)^2/2)

    Est_CPUE <- aggregate(Estimated_CPUE~ CPUE_DF$Year_Month, FUN=mean, na.rm=T)
    Est_CPUE025 <- aggregate(Est_CPUE025~ CPUE_DF$Year_Month, FUN=mean, na.rm=T)
    Est_CPUE975 <- aggregate(Est_CPUE975~ CPUE_DF$Year_Month, FUN=mean, na.rm=T)

    SD <- predict(mod, se.fit=T)
    logSD <- aggregate(SD$se.fit~ CPUE_DF$Year_Month, FUN=mean, na.rm=T)

    CPUE_DF <- CPUE_DF %>% distinct(Date) %>%
      mutate(CPUE=Est_CPUE$Estimated_CPUE/mean(Est_CPUE$Estimated_CPUE),
             logSD=logSD$`SD$se.fit`)

    CPUE_DF$Year <- lubridate::year(CPUE_DF$Date)
    CPUE_DF$Month <- lubridate::month(CPUE_DF$Date)

    CPUE_DF <- left_join(Date_DF, CPUE_DF, by = join_by(Year, Month, Date))
    keep <- range(which(!is.na(CPUE_DF$CPUE)))
    keep_rows <- keep[1]:keep[2]
    CPUE_DF <- CPUE_DF[keep_rows,]
  }




  # ---- Return data.frames
  out <- list()
  out$CAW <- WeightData_Binned
  out$Effort <- Effort_DF
  out$CPUE <- CPUE_DF
  out$Sites <- Sites
  out
}

#' Summarize Case Study Data
#'
#' @param data a data.frame returned by `Process_Data`
#'
#' @return A data.frame
#' @export
Summarize_Data <- function(village, Case_Study_Data, Case_Study_Sites) {

  # ---- Combine adjacent villages
  agg_sites <- Case_Study_Sites %>% filter(Village==village) %>% select(Combined)

  if (nrow(agg_sites)<1) {
    agg_villages <- village
    site_village <- village

  } else {
    agg_villages <- Case_Study_Sites %>% filter(Combined == agg_sites$Combined) %>% select(Village)
    agg_villages <- sort(agg_villages$Village)
    site_village <- village
  }

  village_data <- Case_Study_Data %>% filter(Village %in% agg_villages)

  # ---- Format Dates ----
  village_data$Year <- village_data$year
  village_data$Date <-  lubridate::date(village_data$Fishing.Date)
  village_data$Date2 <-  lubridate::my(paste(village_data$Month, village_data$Year, sep="-"))
  village_data$Month_Name <- month.abb[village_data$Month]

  # ---- Catch-Weight Data ----
  weight_data <- village_data %>% select(Date2,
                                         Sex,
                                         Length=Mantle.length..Cm.,
                                         Weight=Octopus.individual.weight..kg.) %>%
    mutate(Weight=as.numeric(Weight)) %>%
    filter(is.na(Length)==FALSE, is.na(Weight)==FALSE,
           Weight>0, Length>0)

  df <- weight_data %>% group_by(Date=Date2) %>% summarize(n_CAW=length(Weight),
                                                           Catch=sum(Weight))

  # Continuous Dates
  df$Year <- lubridate::year(df$Date)
  df$Month <- lubridate::month(df$Date)
  Years <- range(df$Year)
  Years <- seq(Years[1], Years[2], by=1)

  Date_DF <- data.frame(Year=Years, Month=rep(1:12, each=length(Years))) %>%
    arrange(Year) %>%
    mutate(Date=lubridate::my(paste(Month, Year, sep="-")))

  df <- left_join(Date_DF, df, by = join_by(Year, Month, Date))

  df$Name <- village
  df$Sites <- paste(agg_villages, collapse=', ')


  # ---- Effort -----
  village_data$Duration <- suppressWarnings(as.numeric(village_data$Fishing.Duration..hour.))
  effort_df <- village_data %>% group_by(Date=Date2) %>% summarize(Effort=mean(Duration, na.rm=TRUE))

  df <- left_join(df, effort_df, by='Date') %>%
    select(Name, Sites, Date, Year, Month, n_CAW, Catch, Effort)

  keep <- range(which(!is.na(df$n_CAW)))
  keep_rows <- keep[1]:keep[2]
  df <- df[keep_rows,]

  df$n_months <- nrow(df)
  df$missing_months <- sum(is.na(df$n_CAW))
  df$n_sites <- length(agg_villages)
  df
}


Make_Data_Objects <- function(i, data_list, Case_Study_Sites, save.dir=NULL) {

  obs_data <- data_list[[i]]
  Sites <- obs_data$Sites

  sites <- strsplit(Sites, ', ')[[1]] %>% trimws()

  data <- Make_Data(TRUE)

  # Metadata
  Location <- Case_Study_Sites %>% filter(Village %in% sites) %>% select(Province.Provinsi)

  if (nrow(Location)<1) {
    Location <- 'Unknown'
  } else {
    Location <- unique(Location$Province.Provinsi)
  }


  data$Metadata$Value <- c(Sites,
                           Location,
                           'Octopus cyanea',
                           'Day octopus',
                           'A. Hordyk (adrian@bluematterscience.com',
                           'June 2023')

  # Life History
  data$Ages <-LifeHistory$Ages
  data$Weight_Age_Mean <- LifeHistory$Weight_Age_Mean
  data$Weight_Age_SD <- LifeHistory$Weight_Age_SD
  data$Maturity_at_Age <- LifeHistory$Maturity_at_Age
  data$M_at_Age <- LifeHistory$M_at_Age
  data$Post_Spawning_Mortality <- LifeHistory$Post_Spawning_Mortality


  # CAW
  Weight_Mids <- unique(obs_data$CAW$Bin)
  Weight_Mids <- Weight_Mids[!is.na(Weight_Mids)]
  by <- Weight_Mids[2] - Weight_Mids[1]
  Weight_Bins <-  seq(Weight_Mids[1]-0.5*by, by=by, length.out=length(Weight_Mids)+1)
  nbins <- length(Weight_Mids)

  Years <- unique(obs_data$CAW$Year)
  Time_DF <- expand.grid(Year=Years, Month=1:12) %>% arrange(Year)
  n_months <- nrow(Time_DF)

  data$Weight_Bins <- Weight_Bins
  data$Weight_Mids <- Weight_Mids

  data$CAW <- matrix(NA, nrow=n_months, ncol=nbins)
  data$CAW_ESS <- rep(NA, n_months)

  nsamp <- rep(NA, n_months)
  for (y in 1:nrow(Time_DF)) {
    temp <- obs_data$CAW %>% filter(Year==Time_DF$Year[y], Month==Time_DF$Month[y])
    if (nrow(temp)>0) {
      data$CAW_ESS[y] <- sum(temp$Count)
      data$CAW[y,] <- temp$Count/ data$CAW_ESS[y]
    }
  }

  # Effort
  data$Effort_Mean <- rep(NA, n_months)
  data$Effort_SD <- rep(NA, n_months)

  for (y in 1:nrow(Time_DF)) {
    if(!class(obs_data$Effort)=='logical') {
      if (nrow(obs_data$Effort)>1) {
        temp2 <- obs_data$Effort %>% filter(Year==Time_DF$Year[y], Month==Time_DF$Month[y])
        if(nrow(temp2)>0) {
          data$Effort_Mean[y] <- temp2$rel_effort
          data$Effort_SD[y] <- temp2$sd_effort
        }
      }
    }
  }



  # Index
  data$Index_Mean <- rep(NA, n_months)
  data$Index_SD <- rep(NA, n_months)

  for (y in 1:nrow(Time_DF)) {
    if(!class(obs_data$CPUE)=='logical') {
      if (nrow(obs_data$CPUE)>1) {
        temp2 <- obs_data$CPUE %>% filter(Year==Time_DF$Year[y], Month==Time_DF$Month[y])
        if(nrow(temp2)>0) {
          data$Index_Mean[y] <- temp2$CPUE
          data$Index_SD[y] <- temp2$logSD
        }
      }
    }
  }

  # Drop first and last NA years
  keep <- range(which(!is.na(data$CAW_ESS)))
  keep_rows <- keep[1]:keep[2]
  data$CAW <- data$CAW[keep_rows,]
  data$CAW_ESS <- data$CAW_ESS[keep_rows]
  data$Year <- Time_DF$Year[keep_rows]
  data$Month <- Time_DF$Month[keep_rows]

  data$Effort_Mean <- data$Effort_Mean[keep_rows]
  data$Effort_SD <- data$Effort_SD[keep_rows]
  data$Index_Mean <-   data$Index_Mean[keep_rows]
  data$Index_SD <-   data$Index_SD[keep_rows]

  if (!is.null(save.dir)) {
    nm <- paste(sites, collapse="_")
    nm <- paste0(nm, '.rdata')
    saveRDS(data, file.path(save.dir, nm))
  }
  invisible(data)
}