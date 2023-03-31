


#' Title
#'
#' @param Simulation
#' @param Est
#' @param sim
#'
#' @return
#' @export
#'
#' @examples
plot_F <- function(Simulation, Est, sim=1) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% Data$Month_ind)
  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Simulated=OM$F_mort,
                   Predicted=Est$F_m)
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))

  df2 <- df %>% tail(12)
  print(median((df2$Predicted-df2$Simulated)/df2$Simulated))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    expand_limits(y=0) +
    theme_bw()
}

#' Title
#'
#' @param Simulation
#' @param Est
#' @param sim
#'
#' @return
#' @export
#'
#' @examples
plot_SB <- function(Simulation, Est, sim=1) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% Data$Month_ind)


  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Simulated=OM$SB_fished/mean(OM$SB_unfished_eq),
                   Predicted=Est$SB_m/mean(Est$SB0_m))
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    expand_limits(y=0) +
    theme_bw()
}

#' Title
#'
#' @param Simulation
#' @param Est
#' @param sim
#'
#' @return
#' @export
#'
#' @examples
plot_N <- function(Simulation, Est, sim=1) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% Data$Month_ind)


  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Simulated=OM$N_fished,
                   Predicted=colSums(Est$N_m))
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    expand_limits(y=0) +
    theme_bw()
}

#' Title
#'
#' @param Simulation
#' @param Est
#' @param sim
#'
#' @return
#' @export
#'
#' @examples
plot_SPR <- function(Simulation, Est, sim=1) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% Data$Month_ind)


  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Simulated=OM$SPR,
                   Predicted=Est$SPR)
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    expand_limits(y=0) +
    theme_bw()
}


#' Title
#'
#' @param Simulation
#' @param Est
#' @param sim
#'
#' @return
#' @export
#'
#' @examples
plot_B <- function(Simulation, Est, sim=1) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% Data$Month_ind)


  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Simulated=OM$B_fished/OM$B_unfished,
                   Predicted=Est$B_m/Est$SB0_m)
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    theme_bw()
}






