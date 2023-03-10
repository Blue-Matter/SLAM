
plot_CPUE_fit <- function(Data, Est) {

  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Observed=Data$CPUE,
                   Predicted=Est$stpredIndex)
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Observed), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    expand_limits(y=0) +
    theme_bw()


}

plot_Effort_fit <- function(Data, Est) {
  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Observed=Data$Effort,
                   Predicted=Est$StEffort)

  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))

  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Observed), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    expand_limits(y=0) +
    theme_bw()

}

plot_CAA_fit <- function(Data, Est, months=1:12) {
  nmonths <- length(months)
  nage <- length(Data$M_at_Age)
  Ages <- 0:(nage-1)
  obs <- Data$CAA[,months]/apply(Data$CAA[,months], 2, sum)
  pred <- Est$predCAA[,months]
  Year <- sort(Data$Year)[months]
  Month <- Data$Month[months]
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=rep(Year, each=nage),
                   Month=rep(Month, each=nage),
                   Year_Month=rep(Year_Month,each=nage),
                   Age=rep(Ages, nmonths),
                   Observed=as.vector(obs),
                   Predicted=as.vector(pred))
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  df$Year <- factor(df$Year, ordered = TRUE)
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)

  df <- df %>% tidyr::pivot_longer(., cols=c(Observed, Predicted))

  ggplot(df, aes(x=Age)) +
    geom_line(aes(y=value, color=name)) +
    facet_grid(Year~Month) +
    expand_limits(y=0) +
    theme_bw()
}

plot_CAW_fit <- function(Data, Est, months=1:12) {
  nmonths <- length(months)
  WghtMids <- Data$WghtMids
  nbins <- length(WghtMids)

  obs <- Data$CAW[,months]/apply(Data$CAW[,months], 2, sum)
  pred <- Est$predCAW[,months]
  Year <- sort(Data$Year)[months]
  Month <- Data$Month[months]
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=rep(Year, each=nbins),
                   Month=rep(Month, each=nbins),
                   Year_Month=rep(Year_Month,each=nbins),
                   Weight=rep(WghtMids, nmonths),
                   Observed=as.vector(obs),
                   Predicted=as.vector(pred))

  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))

  df$Year <- factor(df$Year, ordered = TRUE)
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)

  df <- df %>% tidyr::pivot_longer(., cols=c(Observed, Predicted))

  ggplot(df, aes(x=Weight)) +
    geom_line(aes(y=value, color=name)) +
    facet_grid(Year~Month) +
    expand_limits(y=0) +
    theme_bw()
}


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






