

library(SLAM)
library(dplyr)
library(ggplot2)
library(ggthemes)

source('R/functions.R')
source('R/PopDynamics.R')

# Application to the Case Study ----

## Base Case ----
data <- SLAM::casestudydata

do_assess <- Assess(data)

# Plot Fits
Years <- rep(2017:2022, each=12)
Months <- 1:12
Dates_DF <- data.frame(Year=Years, Month=Months)
Dates_DF$Year_Month <- paste(Dates_DF$Year, Dates_DF$Month, sep="_")
Dates_DF$Date <- lubridate::my(paste(Dates_DF$Month, Dates_DF$Year, sep="-"))

nts <- length(do_assess$rep$F_m)
Dates_DF <- Dates_DF[1:nts,]
Dates_DF$t <- 1:nts

# Plot fits to CAW data
Wmids <- do_assess$obj$env$data$WghtMids
CAW_obs <- do_assess$obj$env$data$CAW
CAW_pred<- do_assess$rep$predCAW
CAW_DF <- data.frame(WMids=Wmids,
                     CAW_obs=as.numeric(CAW_obs),
                     CAW_pred=as.numeric(CAW_pred),
                     t=rep(1:nts, each=length(Wmids)))

CAW_DF <- left_join(CAW_DF, Dates_DF)
CAW_DF$Month <- month.abb[CAW_DF$Month]

# drop months with no observations
CAW_DF <- CAW_DF %>% group_by(Date) %>%
  mutate(n=sum(CAW_obs)) %>%
  filter(n>0)

CAW_DF$Year_Month <- paste(CAW_DF$Month, CAW_DF$Year, sep="-")
CAW_DF$Year_Month  <- factor(CAW_DF$Year_Month, levels=unique(CAW_DF$Year_Month),
                             ordered=TRUE)
ggplot(CAW_DF, aes(x=WMids)) +
  facet_wrap(~Year_Month, ncol=7, scales='free_y') +
  geom_bar(aes(y=CAW_obs),  stat='identity') +
  geom_line(aes(y=CAW_pred*n)) +
  theme_clean()+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        plot.background=element_rect(color='white')) +
  labs(y="Frequency", y='Weight Bins (kg)')

# Plot fits to CPUE and Effort data
DF_cpue  <- data.frame(Observed=do_assess$obj$env$data$CPUE,
                       Estimated=do_assess$rep$stpredCPUE,
                       t=1:nts)
DF_cpue$Estimated[is.na(DF_cpue$Observed)] <- NA
DF_cpue <- tidyr::pivot_longer(DF_cpue, cols=1:2, values_to="values")
DF_cpue$Type <- "CPUE"

DF_eff  <- data.frame(Observed=do_assess$obj$env$data$Effort,
                       Estimated=do_assess$rep$StEffort,
                       t=1:nts)
DF_eff$Estimated[is.na(DF_eff$Observed)] <- NA
DF_eff <- tidyr::pivot_longer(DF_eff, cols=1:2, values_to="values")
DF_eff$Type <- "Effort"

DF_cpue_eff <- bind_rows(DF_eff, DF_cpue)
DF_cpue_eff <- left_join(DF_cpue_eff, Dates_DF)


ggplot(DF_cpue_eff, aes(x=Date, y=values, color=name, linetype=name)) +
  facet_wrap(~Type, scales="free_y") +
  geom_line(size=1.2) +
  theme_clean()+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        plot.background=element_rect(color='white'),
        legend.background = element_rect(color='white')) +
  scale_x_date(date_breaks = "6 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  labs(color='', x='', y='Value', linetype='')


DF_cpue_eff$Month <- month.abb[DF_cpue_eff$Month]


DF_cpue_eff2 <- DF_cpue_eff %>% tidyr::pivot_longer(cols=c(1,2))
DF_cpue_eff2 <- DF_cpue_eff2 %>% tidyr::pivot_longer(cols=c(1,2),
                                                     values_to='Effortv',
                                                     names_to = 'Effort')


plot(do_assess$obj$env$data$CPUE, type='l')
lines(do_assess$rep$stpredCPUE, col='blue')

ind <- which(!is.na(do_assess$obj$env$data$CPUE))
mean(do_assess$rep$stpredCPUE[ind])

ind <- which(!is.na(do_assess$obj$env$data$Effort))
mean(do_assess$rep$StEffort[ind])

mean(do_assess$rep$StEffort)


plot(do_assess$obj$env$data$Effort, type='l')
lines(do_assess$rep$StEffort, col='blue')






# Plot estimates

df_rec <- data.frame(Month=month.abb, Rec=do_assess$rep$R0_m)
df_rec$Month <- factor(df_rec$Month, levels=month.abb, ordered=TRUE)

ggplot(df_rec, aes(x=Month, y=Rec, group=1)) +
  geom_line(size=1.2) +
  expand_limits(y=0) +
  theme_clean()+
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_rect(color='white'),
        legend.background = element_rect(color='white')) +
  labs(y="Estimated Relative Recruitment")

df_Sel <- data.frame(Age=0:14, Selectivity=do_assess$rep$selA)
ggplot(df_Sel, aes(x=Age, y=Selectivity)) +
  geom_line(size=1.2) +
  expand_limits(y=0) +
  theme_clean()+
  theme(plot.background=element_rect(color='white')) +
  labs(y="Estimated Selectivity")


estsDF <- Dates_DF
estsDF$SPR <- do_assess$rep$SPR
estsDF$F <- do_assess$rep$F_m
ind <- which(is.na(do_assess$obj$env$data$Effort))
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
ind <- which(is.na(do_assess$obj$env$data$CPUE))
estF <- do_assess$rep$F_m
estF[ind] <- NA
estSPR <- do_assess$rep$SPR
estSPR[ind] <- NA
estCB <- do_assess$rep$predCB
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
opt <- Optimize(data, do_assess$rep$R0_m, do_assess$rep$selA, assumed_h=data$h,
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

## Sensitivity Tests ----

# Different M
# Different h
# Different Maturity
# Different Growth


