#' Process Case Study Data
#'
#' @param village Name of the village
#' @param Case_Study_Data A data.frame
#'
#' @return A list
#' @export
Process_Data <- function(village, Case_Study_Data) {
  village_data <- Case_Study_Data %>% filter(Village==village)

  # ---- Format Dates ----
  village_data$Year <- village_data$year
  village_data$Date <-  lubridate::date(village_data$Fishing.Date)
  village_data$Date2 <-  lubridate::my(paste(village_data$Month, village_data$Year, sep="-"))
  village_data$Month_Name <- month.abb[village_data$Month]

  # --- Effort Data ----
  village_data$Duration <- as.numeric(village_data$Fishing.Duration..hour.)

  effort_data <- village_data %>% filter(is.na(Duration) ==FALSE)
  if (nrow(effort_data)<1) {
    Effort_Month_Year <- NA
  } else {
    # only keep duration < 10 hours
    # effort_data <- effort_data %>% filter(Duration <=10)

    Effort_Month_Year <- effort_data %>%
      distinct(Date, Year, Month, Month_Name, Fisher, Duration, Date2, Only.Fishing=Only.fishing.) %>%
      group_by(Date2, Year, Month, Month_Name) %>%
      summarise(Effort=sum(Duration, na.rm=TRUE), .groups='keep')

    Effort_Month_Year$Village <- village
    Effort_Month_Year$Month_Name <- factor(    Effort_Month_Year$Month_Name, levels=month.abb, ordered = TRUE)
  }



  # ---- Weight Composition Data -----
  weight_data <- village_data %>% select(Year, Month, Month_Name,
                                         Sex,
                                         Length=Mantle.length..Cm.,
                                         Weight=Octopus.individual.weight..kg.) %>%
    mutate(Length=as.numeric(Length), Weight=as.numeric(Weight)) %>%
    filter(is.na(Length)==FALSE, is.na(Weight)==FALSE,
           Weight>0, Length>0)

  maxWeightBin <- 4

  by <- 0.25
  Weight_Bins <- seq(0, to=maxWeightBin, by=by)
  Weight_Mids <- seq(0.5*by, by=by, length.out=length(Weight_Bins)-1)

  weight_data$Weight[weight_data$Weight>=maxWeightBin] <- maxWeightBin
  WeightData_Binned <- weight_data %>% group_by(Year, Month, Month_Name) %>%
    reframe(bin = cut(Weight, breaks = Weight_Bins))

  WeightData_Binned <- WeightData_Binned %>% group_by(Year, Month, Month_Name) %>%
    reframe(Count=as.numeric(table(bin)), Bin=Weight_Mids)

  WeightData_Binned$Village <- village

  WeightData_Binned$Month_Name <- factor(WeightData_Binned$Month_Name, levels=month.abb, ordered = TRUE)

  # ---- Return data.frames
  out <- list()
  out$Effort <- Effort_Month_Year
  out$CAW <- WeightData_Binned
  out$Village <- village

  out
}

#' Summarize Case Study Data
#'
#' @param data a data.frame returned by `Process_Data`
#'
#' @return A data.frame
#' @export
Make_Summary_Table <- function(data) {

  Site_Name <- data$Village

  # Summarize CAW
  CAW_df <- data$CAW %>% group_by(Site=Village, Year, Month_Name ) %>% summarize(CAW_n=sum(Count))

  # Summarize Effort
  if (is.null(dim(data$Effort))) {
    df <- CAW_df
    df$Effort <- NA
  } else {
    Effort_df <- data$Effort %>% ungroup() %>% select(Site=Village, Year, Month_Name, Effort=Effort )
    df <- left_join(CAW_df, Effort_df, by = join_by(Site, Year, Month_Name))
  }

  years <- df$Year %>% range()
  years <- unique(years)
  Time_DF <- expand.grid(Year=years, Month_Name=month.abb) %>% arrange(Year)
  Time_DF$Month_Name <- factor(Time_DF$Month_Name, levels=month.abb, ordered = TRUE)
  Time_DF$Site <- Site_Name
  df2 <- left_join(Time_DF, df, by = join_by(Site, Year, Month_Name))

  df2 %>% select(Site, Year, Month=Month_Name, CAW_n, Effort)

}

