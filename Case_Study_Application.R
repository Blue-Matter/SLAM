
# To do:

# - F_minit has massive SD -
# check selectivy parameters gradients
# update Methods for selectity-at-weight (if this is fixed above)

# -
BaseCase$sdreport



library(SLAM)
library(dplyr)
library(ggplot2)
library(ggthemes)

source('R/functions.R')
source('R/PopDynamics.R')

# Application to the Case Study ----

## Base Case
data <- SLAM::casestudydata

## Lower M
data_lowM <- data
data_lowM$M_at_Age <- rep(0.1, length(data_lowM$M_at_Age))

## Higher M
data_highM <- data
data_highM$M_at_Age <- rep(0.2, length(data_lowM$M_at_Age))

## Lower h
data_lowh <- data
data_lowh$h <- 0.6

## Higher h
data_highh <- data
data_highh$h <- 0.8

## Smaller size
data_smaller <- data
data_smaller$Weight_Age <- data$Weight_Age * 0.9

## Larger size
data_larger <- data
data_larger$Weight_Age <- data$Weight_Age * 1.1

BaseCase <- Assess(data)
LowerM <- Assess(data_lowM)
HigherM <- Assess(data_highM)
Lowerh <- Assess(data_lowh)
Higherh <- Assess(data_highh)
SmallerW <- Assess(data_smaller)
LargerW <- Assess(data_larger)





Process_Assess <- function(assess, Name='Base Case', firstyr=2017) {
  nts <- length(assess$rep$SPR)
  nyr <- ceiling(nts/12)
  years <- firstyr:(firstyr+nyr-1)
  Years <- rep(years, each=12)
  Months <- 1:12
  Dates_DF <- data.frame(Year=Years, Month=month.abb[Months])
  Dates_DF$Year_Month <- paste(Dates_DF$Month, Dates_DF$Year, sep="-")
  Dates_DF$Year_Month <- factor(Dates_DF$Year_Month,
                                levels=unique(Dates_DF$Year_Month),
                                ordered=TRUE)
  Dates_DF$Date <- lubridate::my(paste(Dates_DF$Month, Dates_DF$Year, sep="-"))
  Dates_DF <- Dates_DF[1:nts,]
  Dates_DF$t <- 1:nts

  # Effort
  DF_eff  <- data.frame(Observed=assess$obj$env$data$Effort,
                        Estimated=assess$rep$StEffort,
                        t=1:nts)
  DF_eff$Estimated[is.na(DF_eff$Observed)] <- NA
  # DF_eff <- tidyr::pivot_longer(DF_eff, cols=1:2, values_to="values")
  DF_eff$Type <- "Effort"
  DF_eff <- left_join(DF_eff, Dates_DF, by='t')
  DF_eff$Name <- Name
  DF_eff$Year_Month  <- factor(DF_eff$Year_Month, levels=unique(DF_eff$Year_Month),
                                ordered=TRUE)


  # CPUE
  DF_cpue  <- data.frame(Observed=assess$obj$env$data$CPUE,
                         Estimated=assess$rep$stpredCPUE,
                         t=1:nts)
  DF_cpue$Estimated[is.na(DF_cpue$Observed)] <- NA
  # DF_cpue <- tidyr::pivot_longer(DF_cpue, cols=1:2, values_to="values")
  DF_cpue$Type <- "CPUE"
  DF_cpue <- left_join(DF_cpue, Dates_DF, by='t')
  DF_cpue$Name <- Name
  DF_cpue$Year_Month  <- factor(DF_cpue$Year_Month, levels=unique(DF_cpue$Year_Month),
                               ordered=TRUE)

  # Catch-at-Weight
  Wmids <- assess$obj$env$data$WghtMids
  CAW_obs <- assess$obj$env$data$CAW
  CAW_pred<- assess$rep$predCAW
  CAW_DF <- data.frame(WMids=Wmids,
                       CAW_obs=as.numeric(CAW_obs),
                       CAW_pred=as.numeric(CAW_pred),
                       t=rep(1:nts, each=length(Wmids)))
  CAW_DF <- left_join(CAW_DF, Dates_DF, by='t')
  CAW_DF <- CAW_DF %>% group_by(Date) %>%
    mutate(n=sum(CAW_obs)) %>%
    filter(n>0)

  CAW_DF$Year_Month <- paste(CAW_DF$Month, CAW_DF$Year, sep="-")
  CAW_DF$Year_Month  <- factor(CAW_DF$Year_Month, levels=unique(CAW_DF$Year_Month),
                               ordered=TRUE)
  CAW_DF$Name <- Name

  # Predicted seasonal recruitment
  df_rec <- data.frame(Month=month.abb, Rec=assess$rep$R0_m)
  df_rec$Month <- factor(df_rec$Month, levels=month.abb, ordered=TRUE)
  df_rec$Name <- Name

  # is data present
  data_na_df <- data.frame(t=1:nts,
                        cpue=is.na(assess$obj$env$data$CPUE),
                        effort=is.na(assess$obj$env$data$Effort),
                        caw=apply(assess$obj$env$data$CAW, 2, sum)==0)

  data_na_df$no_data <- !data_na_df %>% apply(1, prod) == 0

  # Predicted fishing mortality
  df_predF <- Dates_DF
  df_predF$predF <-assess$rep$F_m
  df_predF$Name <- Name
  df_predF <- left_join(df_predF, data_na_df, by='t')

  # Predicted SPR
  df_predSPR <- Dates_DF
  df_predSPR$SPR <-assess$rep$SPR
  df_predSPR$Name <- Name
  df_predSPR <- left_join(df_predSPR, data_na_df, by='t')

  # Predicted selectivity-at-age
  selA <- assess$rep$selA
  ages <- seq(0, by=1, length.out=length(selA))
  df_predSatA <- data.frame(Age=ages, selA=selA)
  df_predSatA$Name <- Name

  # Optimize for F
  utilpow <- 0.4
  opt <- Optimize(assess$obj$env$data, assess$rep$R0_m,
                  assess$rep$selA, assumed_h=assess$obj$env$data$h,
                  utilpow = utilpow)

  opt_df <- data.frame(Month=month.abb[1:12],
                       F=opt$F_m,
                       C=opt$predCB,
                       SPR=opt$SPR,
                       utilpow=utilpow,
                       Name=Name)
  out <- list()
  out$Effort <- DF_eff
  out$CPUE <- DF_cpue
  out$CAW <- CAW_DF
  out$predR0m <- df_rec
  out$predF <- df_predF
  out$predSPR <- df_predSPR
  out$predSatA <- df_predSatA
  out$opt <- opt_df
  out
}


# --- Plot Fits ----

BaseCase_Results <- Process_Assess(BaseCase)
LowerM_Results <- Process_Assess(LowerM, 'Lower M')
HigherM_Results <- Process_Assess(HigherM, 'Higher M')
Lowerh_Results <- Process_Assess(Lowerh, 'Lower h')
Higherh_Results <-Process_Assess(Higherh, 'Higher h')
SmallerW_Results <- Process_Assess(SmallerW, 'Smaller size')
LargerW_Results <- Process_Assess(LargerW, 'Larger size')

assessList <- list(BaseCase_Results,
                   LowerM_Results, HigherM_Results,
                   Lowerh_Results, Higherh_Results,
                   SmallerW_Results, LargerW_Results)


plot_CPUE <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$CPUE$Name <- factor(res$CPUE$Name, levels=unique(res$CPUE$Name), ordered = TRUE)

  ggplot(res$CPUE) +
    geom_line(aes(x=Year_Month, y=Observed, group=1), size=1.1) +
    geom_point(aes(x=Year_Month, y=Observed, group=1), size=2) +
    geom_line(aes(x=Year_Month, y=Estimated, group=Name, linetype=Name,
                  color=Name)) +
    labs(x='Date', y='CPUE', color='Run Name', linetype='Run Name') +
    expand_limits(y=0) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))
}

plot_Effort <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$Effort$Name <- factor(res$Effort$Name, levels=unique(res$Effort$Name), ordered = TRUE)

  ggplot(res$Effort) +
    geom_line(aes(x=Year_Month, y=Observed, group=1), size=1.1) +
    geom_point(aes(x=Year_Month, y=Observed, group=1), size=2) +
    geom_line(aes(x=Year_Month, y=Estimated, group=Name, linetype=Name,
                  color=Name)) +
    expand_limits(y=0) +
    labs(x='Date', y='Relative Effort', color='Run Name', linetype='Run Name') +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))

}

plot_CAW <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$CAW$Name <- factor(res$CAW$Name, levels=unique(res$CAW$Name), ordered = TRUE)
  res$CAW <- res$CAW %>% group_by(Name, Date) %>%
    mutate(CAW_pred2=CAW_pred*n)

  ggplot(res$CAW, aes(x=WMids)) +
    facet_wrap(~Year_Month, ncol=7, scales='free_y') +
    geom_bar(aes(y=CAW_obs),  stat='identity', data=res$CAW %>% filter(Name=='Base Case')) +
    geom_line(aes(y=CAW_pred2, color=Name, linetype=Name, group=Name)) +
    labs(y="Frequency", x='Weight Bins (kg)', color='Run Name', linetype='Run Name') +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))

}

# --- Plot Estimates ----
plot_predR0m <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$predR0m$Name <- factor(res$predR0m$Name, levels=unique(res$predR0m$Name), ordered = TRUE)

  ggplot(res$predR0m, aes(x=Month, y=Rec, linetype=Name, color=Name, group=Name)) +
    geom_line() +
    labs(y="Estimated Relative Recruitment", x='Month', color='Run Name', linetype='Run Name') +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))
}


plot_predF <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$predF$Name <- factor(res$predF$Name, levels=unique(res$predF$Name), ordered = TRUE)

  # remove estimates when no data
  res$predF$predF[res$predF$no_data] <- NA

  ggplot(res$predF, aes(x=Year_Month, y=predF, linetype=Name, color=Name, group=Name)) +
    geom_line() +
    labs(y="Estimated Fishing Mortality", x='Date',
         color='Run Name', linetype='Run Name') +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))
}

plot_predSPR <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$predSPR$Name <- factor(res$predSPR$Name, levels=unique(res$predSPR$Name), ordered = TRUE)

  # remove estimates when no data
  res$predSPR$SPR[res$predSPR$no_data] <- NA

  ggplot(res$predSPR, aes(x=Year_Month, y=SPR, linetype=Name, color=Name, group=Name)) +
    geom_line() +
    labs(y="Estimated SPR", x='Date',
         color='Run Name', linetype='Run Name') +
    expand_limits(y=c(0,1)) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))
}

plot_predSelect <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$predSatA$Name <- factor(res$predSatA$Name, levels=unique(res$predSatA$Name), ordered = TRUE)
  ggplot(res$predSatA, aes(x=Age, y=selA, linetype=Name, color=Name, group=Name)) +
    geom_line() +
    labs(y="Estimated Selectivity", x='Age (month)',
         color='Run Name', linetype='Run Name') +
    expand_limits(y=c(0,1)) +
    theme_clean() +
    theme(plot.background=element_rect(color='white'),
          legend.background = element_rect(color='white'))

}

plot_predOpt <- function(assessList) {
  res <- purrr::pmap(assessList, bind_rows)
  res$opt$Name <- factor(res$opt$Name, levels=unique(res$opt$Name), ordered = TRUE)
  res$opt <- res$opt %>% tidyr::pivot_longer(cols=2:4)
  res$opt$Month <- factor(res$opt$Month, levels=month.abb, ordered=TRUE)

  ggplot(res$opt, aes(x=Month, y=value, color=Name, linetype=Name,
                      group=Name)) +
    facet_wrap(~name, scales='free_y') +
    expand_limits(y=0) +
    geom_line()

}



out$predF <- df_predF
out$predSPR <- df_predSPR
out$predSatA <- df_predSatA



df_Sel <- data.frame(Age=0:14, Selectivity=BaseCase$rep$selA)
ggplot(df_Sel, aes(x=Age, y=Selectivity)) +
  geom_line(size=1.2) +
  expand_limits(y=0) +
  theme_clean()+
  theme(plot.background=element_rect(color='white')) +
  labs(y="Estimated Selectivity")


estsDF <- Dates_DF
estsDF$SPR <- BaseCase$rep$SPR
estsDF$F <- BaseCase$rep$F_m
ind <- which(is.na(BaseCase$obj$env$data$Effort))
estsDF$SPR[ind] <- NA
estsDF$F[ind] <- NA

estsDF <- estsDF %>% tidyr::pivot_longer(cols=6:7)
ggplot(estsDF, aes(x=Date, y=value, color=name, linetype=name)) +
  geom_line(size=1.2) +
  labs(x="", y='Estimated Value',
       color='', linetype='') +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  expand_limits(y=c(0,1)) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_rect(color='white'),
        legend.background = element_rect(color='white'))



# Calculate mean monthly Fishing mortality
ind <- which(is.na(BaseCase$obj$env$data$CPUE))
estF <- BaseCase$rep$F_m
estF[ind] <- NA
estSPR <- BaseCase$rep$SPR
estSPR[ind] <- NA
estCB <- BaseCase$rep$predCB
estCB[ind] <- NA


df <- data.frame(estF=estF, estSPR=estSPR, estCB=estCB)
df$t <- 1:nrow(df)
df$M <- df$t %%12
df$M[df$M==0] <- 12
df$Month <- month.abb[df$M]
df$Month <- factor(df$Month, levels=month.abb, ordered=TRUE)
df2 <- df %>% group_by(Month) %>% summarize(F=mean(estF, na.rm=TRUE),
                                            SPR=mean(estSPR, na.rm=TRUE),
                                            CB=mean(estCB, na.rm=TRUE))

# optimizing for utility
opt <- Optimize(data, BaseCase$rep$R0_m, BaseCase$rep$selA, assumed_h=data$h,
                utilpow = 0.4)

df2a<- df2 %>% tidyr::pivot_longer(cols=2:4, values_to='Estimated')

df2b <- df2
df2b$F <- opt$F_m
df2b$SPR <- opt$SPR
df2b$CB <- opt$predCB
df2b<- df2b %>% tidyr::pivot_longer(cols=2:4, values_to='Predicted Optimal')

df3 <- left_join(df2a, df2b, by=c("Month", "name"))
df4 <- df3 %>% tidyr::pivot_longer(cols=3:4, names_to='Name')

ggplot(df4, aes(x=Month, group=Name, color=Name)) +
  facet_wrap(~name, scales='free') +
  expand_limits(y=0) +
  geom_line(aes(y=value), size=1.2) +
  theme_clean()+
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_rect(color='white'),
        legend.background = element_rect(color='white'))



