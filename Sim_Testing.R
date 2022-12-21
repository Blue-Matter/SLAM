

library(SLAM)
library(dplyr)
library(ggplot2)
library(ggthemes)

source('R/functions.R')

# Simulation Testing -----

data <- SLAM::casestudydata


## Seasonal Recruitment Scenarios ----

mulist <- list(1, 6.5, 6.5, c(3,9))
sdlist <- list(100, 1, 3, c(1.5,1.5))
reclist <- list()
for (i in 1:length(mulist)) {
  reclist[[i]] <- GenMonthlyRec(mulist[[i]], sdlist[[i]])
}

rec_scen_names <- c('Constant', 'Pulse', 'Diffuse', 'Bi-Modal')
rec_scen_names <- factor(rec_scen_names, ordered=TRUE, levels=rec_scen_names)
rec_pattern_df <- data.frame(Scenario=rep(rec_scen_names, each=12), M=1:12, Rec=unlist(reclist))
rec_pattern_df$Month <- month.abb[rec_pattern_df$M]
rec_pattern_df$Month <- factor(rec_pattern_df$Month, ordered = TRUE, levels=month.abb)

ggplot(rec_pattern_df, aes(x=Month, y=Rec, group=Scenario)) +
  facet_wrap(~Scenario) +
  geom_line() +
  labs(x='Calendar Month', y="Relative Recruitment") +
  theme_clean() +
  theme(plot.background=element_rect(colour = "white"),
        axis.text.x = element_text(angle = 60, hjust=1))

ggsave('Figures/SimTest/Rec_Scenarios.png', width=4, height=4)




## --- Exploitation History ---
nyears <- 10
nts <- nyears * 12
nsim <- 100
set.seed(1001)

calcMonthlyMean <- function(dat, offset=0) {

  nts <- length(dat)
  df <- data.frame(t=1:nts, n=dat)
  df$n[df$n==0] <- NA
  df$M <- (df$t+ offset) %% 12
  df$M[df$M ==0] <- 12
  df$Month <- month.abb[df$M ]
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)
  df %>% group_by(Month) %>%
    summarize(mu=mean(n, na.rm=TRUE), sd=sd(log(n), na.rm=T))
}

mu_sd <- calcMonthlyMean(data$Effort)
mu_sd$sd <- 0.3

Effort <- matrix(rlnorm(nts*nsim, log(mu_sd$mu)-0.5*mu_sd$sd^2, mu_sd$sd), nrow=nts, ncol=nsim)

Effort_DF <- data.frame(t=1:nts,
                        Effort=as.numeric(Effort),
                        Sim=rep(1:nsim, each=nts),
                        Year=rep(1:nyears, each=12))
Effort_DF$M <- Effort_DF$t %% 12
Effort_DF$M[Effort_DF$M==0] <- 12
Effort_DF$Month <- month.abb[Effort_DF$M]
Effort_DF$Year <- Effort_DF$Year + 2011
Effort_DF$Date <- lubridate::my(paste(Effort_DF$Month, Effort_DF$Year, sep="-"))

ggplot(Effort_DF, aes(x=Date, y=Effort, group=Sim)) +
  geom_line(alpha=0.1) +
  geom_line(data=Effort_DF %>% filter(Sim==1), size=1, color='blue') +
  labs(x="", y='Relative Effort') +
  scale_x_date(date_breaks = "6 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_rect(color='white'))

ggsave('Figures/SimTest/Effort.png', width=6, height=4)


# ---- Make Parameters List ----
Scens <- list()

for (s in 1:length(rec_scen_names)) {
  Pars <- list()
  Pars$Name <- rec_scen_names[s]
  Pars$maxage <- 14
  Pars$M <- 0.15
  Pars$Weight_Age <- data$Weight_Age
  Pars$Weight_Age_SD <- data$Weight_Age_SD
  Pars$Mat_at_Age <- data$Mat_at_Age
  Pars$PSM_at_Age <- data$PSM_at_Age
  Pars$h <- data$h
  Pars$WghtBins <-  data$WghtBins
  Pars$WghtMids <-  data$WghtMids

  # recruitment pattern
  ind <- switch(as.character(rec_scen_names[s]), Constant=1, Pulse=2, Diffuse=3, 'Bi-Modal'=4)
  Pars$rec_mu <- reclist[[ind]]
  Pars$Rbar <- 1E5

  # recruitment deviations
  rec_sd <- 0.9
  rec_devs <- exp(rnorm(nts*nsim, -0.5*rec_sd^2, rec_sd))
  Pars$rec_devs <- matrix(rec_devs, nrow=nts, ncol=nsim)

  # selectivity
  Pars$sA50 <- 5
  Pars$sA95 <- 6

  # effort & fishing mortality
  Pars$Effort <- Effort_DF
  Pars$q <- 0.2

  # Observation process
  Pars$nyears <- 5 # number of years of data
  Pars$sigmaI <- 0.2
  Pars$sigmaE <- 0.2

  Pars$CAW_nsamp <- 1000
  Pars$CAW_ESS <- 200

  Scens[[s]] <- Pars
}

# ---- Explore Simulation Properties ----



# ---- Simulate & Estimate ----

for (s in 1:length(Scens)) {
  message('Scenario: ', s,'/', length(Scens))
  # loop over scenarios
  Pars <- Scens[[s]]
  Name <- Pars$Name

  # loop over simulations
  scen_list <- list()
  for (i in 1:nsim) {
    message('Sim: ', i,'/', nsim)
    # simulate
    SimPop <- Simulate(Pars, sim=i)

    # assess
    sim_data <- list()
    nts <- length(SimPop$Biomass)
    ndata_ts <- length(SimPop$Eff_ind)
    data_ts <- (nts-ndata_ts+1):nts
    sim_data$Weight_Age <- Pars$Weight_Age
    sim_data$Weight_Age_SD <- Pars$Weight_Age_SD
    sim_data$M_at_Age <- SimPop$M_at_Age
    sim_data$Mat_at_Age <- Pars$Mat_at_Age
    sim_data$PSM_at_Age <- Pars$PSM_at_Age
    sim_data$h  <- Pars$h
    sim_data$Effort <- SimPop$Eff_ind
    sim_data$Effort_SD <- rep(0.3, ndata_ts)
    sim_data$CPUE <- SimPop$Index
    sim_data$CPUE_SD <- rep(0.3, ndata_ts)
    sim_data$WghtBins <- Pars$WghtBins
    sim_data$WghtMids <- Pars$WghtMids
    sim_data$CAW <- SimPop$CAW_samp
    sim_data$CAW_ESS <- rep(Pars$CAW_ESS, ndata_ts)
    sim_data$model <- 'SLAM'
    sim_data$Fit_Effort <- 1
    sim_data$Fit_CPUE <- 1
    sim_data$use_Frwpen <- 1
    sim_data$use_R0rwpen <- 1
    sim_data$use_Fmeanprior <- 0
    sim_data$F_meanprior <- 1

    do_assess <- try(Assess(sim_data), silent=TRUE)

    scen_list[[i]] <- list(Sim=SimPop, Assess=do_assess, sim_data=sim_data)
  }
  fl <- file.path(paste0('Results/Sim_Test/', Name, '.rda'))
  saveRDS(scen_list, fl)

}

# ---- Summarise Performance ----
scenlist <- list()
for (scen in 1:length(rec_scen_names)) {
  message('Scenario:', scen)

  Name <- rec_scen_names[scen]
  fl <-  file.path(paste0('Results/Sim_Test/', Name, '.rda'))
  obj <- readRDS(fl)

  Rec_scen <- as.character(Name)

  dflist <- list()
  for (sim in 1:nsim) {
    message('Sim:', sim)
    SimPop <- obj[[sim]]$Sim
    Assess_obj <- obj[[sim]]$Assess

    nts <- length(SimPop$Biomass)
    ndata_ts <- length(SimPop$Eff_ind)
    data_ts <- (nts-ndata_ts+1):nts

    dflist[[sim]] <- data.frame(R0_act=SimPop$Rec_Pattern,
                                R0_pred= Assess_obj$rep$R0_m,
                                SPR_act=SimPop$SPR[data_ts],
                                SPR_pred=Assess_obj$rep$SPR,
                                F_act=SimPop$F_m[data_ts],
                                F_pred=Assess_obj$rep$F_m,
                                Sim=sim,
                                Rec_scen=Rec_scen,
                                t=1:ndata_ts)
  }
  scenlist[[scen]] <- do.call('rbind', dflist)
}


DF <- do.call('rbind', scenlist)
DF$m <- 1:12
DF$Month <- month.abb[DF$m]
DF$Month <- factor(DF$Month, levels=month.abb, ordered=TRUE)
DF$Rec_scen <- factor(DF$Rec_scen, levels=unique(rec_scen_names), ordered = TRUE)

# Drop Sims where all SPR estimates >0.95 # failed to converge - would be detected in real applications
# look at fix with SPR or F prior/penalty

tt <- DF %>% group_by(Rec_scen, Sim) %>%
  mutate(highSPR=SPR_pred>0.95) %>%
  group_by(Rec_scen, Sim) %>%
  summarise(highSPR=prod(highSPR)) %>%
  filter(highSPR==TRUE)

tt$Rec_scen_Sim <- paste0(tt$Rec_scen, tt$Sim)
DF$Rec_scen_Sim <- paste0(DF$Rec_scen, DF$Sim)

DF <- DF %>% filter(!Rec_scen_Sim %in%tt$Rec_scen_Sim)

DF_true <- DF %>% filter(Sim==1)

ggplot(DF, aes(x=Month, y=R0_pred, group=Sim)) +
  facet_wrap(~Rec_scen) +
  geom_point(stat='summary', fun=mean, color='white') +
  stat_summary(fun=mean, geom="line", alpha=0.3) +
  geom_line(aes(x=Month, y=R0_act, group=Sim), color='blue', size=1,
            data=DF_true) +
  labs(y="Relative Monthly Recruitment") +
  expand_limits(y=0) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
        plot.background=element_rect(color='white'))

ggsave(paste0('Figures/SimTest/R0_ests.png'), width=6, height=4)



DF <- DF %>% group_by(Rec_scen, Sim) %>%
  mutate(SPR_RE=(SPR_pred-SPR_act)/SPR_act,
         F_RE=(F_pred-F_act)/F_act)


DF %>% group_by(Rec_scen) %>% summarise(F=median(F_RE),
                                        SPR=median(SPR_RE))


DF_SPR <- DF %>% group_by(Rec_scen, Sim) %>%
  summarize(SPR_MRE=median(SPR_RE))

ggplot(DF_SPR, aes(x=Rec_scen, y=SPR_MRE)) +
  geom_hline(yintercept = 0, linetype=2) +
  expand_limits(y=c(-1,1.5)) +
  geom_boxplot(fill='lightgray') +
  theme_clean() +
  labs(x='Recruitment Scenario', y='Median Relative Error SPR') +
  theme(plot.background=element_rect(color='white'),
        legend.position = 'bottom')

ggsave('Figures/SimTest/SPR_MRE.png', width=6, height=4)


DF_F <- DF %>% group_by(Rec_scen, Sim) %>%
  summarize(F_MRE=median((F_pred-F_act)/F_act))

tt <- DF %>% filter(Sim==1)

plot(tt$t, tt$F_act, type='l')
lines(tt$t, tt$F_pred, col='blue')

median((tt$F_pred-tt$F_act)/tt$F_act)


ggplot(DF_F, aes(x=Rec_scen, y=F_MRE)) +
  geom_hline(yintercept = 0, linetype=2) +
  expand_limits(y=c(-1,1)) +
  geom_boxplot(fill='lightgray') +
  theme_clean() +
  labs(x='Recruitment Scenario', y='Median Relative Error Fishing Mortality (F)') +
  theme(plot.background=element_rect(color='white'),
        legend.position = 'bottom')

ggsave('Figures/SimTest/F_MRE.png', width=6, height=4)

DF_F %>% group_by(Rec_scen) %>% summarise(median(F_MRE))
DF_SPR %>% group_by(Rec_scen) %>% summarise(median(SPR_MRE))

ggplot(DF, aes(x=SPR_act, y=SPR_pred)) +
  facet_wrap(~Rec_scen, ncol=2) +
  geom_point(alpha=0.3) +
  expand_limits(y=c(0,1), x=c(0,1)) +
  labs(x='SPR OM', y='SPR Estimated') +
  geom_abline(slope=1, intercept = 0, linetype=2) +
  theme_clean() +
  theme(plot.background=element_rect(color='white'),
        legend.position = 'bottom')
ggsave('Figures/SimTest/SPR_compare.png', width=6, height=6)


axmax <- max(c(DF$F_act, DF$F_pred))

ggplot(DF, aes(x=F_act, y=F_pred)) +
  facet_wrap(~Rec_scen, ncol=2) +
  geom_point(alpha=0.6) +
  expand_limits(y=c(0,axmax), x=c(0,axmax)) +
  labs(x='F OM', y='F Estimated') +
  geom_abline(slope=1, intercept = 0, linetype=2) +
  theme_clean() +
  theme(plot.background=element_rect(color='white'),
        legend.position = 'bottom')
ggsave('Figures/SimTest/F_compare.png', width=6, height=6)

# ---- Calculate optimal F pattern for Recruitment Scenarios ----

scenlist <- list()
for (scen in 1:length(rec_scen_names)) {
  message('Scenario:', scen)
  Name <- as.character(rec_scen_names[scen])
  fl <-  file.path(paste0('Results/Sim_Test/', Name, '.rda'))
  obj <- readRDS(fl)

  sim <- 1
  SimPop <- obj[[sim]]$Sim

  Data_true <- list()
  Data_true$Weight_Age <- SimPop$Pars$Weight_Age
  Data_true$R0 <- SimPop$Pars$Rbar
  nage <- length(Data_true$Weight_Age)
  Data_true$Mat_at_Age <- SimPop$Pars$Mat_at_Age
  Data_true$M_at_Age <- rep(SimPop$Pars$M, nage)
  Data_true$PSM_at_Age <- SimPop$Pars$PSM_at_Age

  utilpow_vec <- seq(0.4, 1, by=0.2)
  util_list <- list()
  for (i in seq_along(utilpow_vec)) {
    opt1 <- Optimize(Data_true, SimPop$Rec_Pattern, SimPop$Sel_at_Age, opt_type=1,
                     utilpow=utilpow_vec[i],
                     assumed_h=SimPop$Pars$h)
    util_list[[i]] <- data.frame(M=1:12, Month=month.abb[1:12],
                                 'Fishing Mortality'=opt1$F_m,
                                 SPR=opt1$SPR,
                                 Catch=opt1$predCB,
                                 utilpow=utilpow_vec[i])
  }
  df <- do.call('rbind', util_list)

  df$Rec_Scen <- Name
  scenlist[[scen]] <- df
}

scen_df <- do.call('rbind', scenlist)
scen_df$Month <- factor(scen_df$Month, levels=month.abb, ordered = TRUE)
scen_df$utilpow <- factor(scen_df$utilpow)
scen_df$Rec_Scen <- factor(scen_df$Rec_Scen, levels=rec_scen_names, ordered = T)

rec_pattern_df2 <- rec_pattern_df
rec_pattern_df2$utilpow <- 1
rec_pattern_df2$Rec_Scen <- rec_pattern_df2$Scenario

ggplot(scen_df, aes(x=Month, group=1)) +
  facet_grid(Rec_Scen~utilpow, scales='free_y') +
  expand_limits(y=c(0, 0.2)) +
  geom_line(aes(y=Fishing.Mortality)) +
  geom_line(data=rec_pattern_df2, aes(y=Rec, group=1), linetype=2) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
        plot.background=element_rect(color='white'),
        legend.background = element_blank(),
        legend.position = 'bottom') +
  labs(y='Fishing Mortality')
ggsave('Figures/SimTest/Fopt_patterns.png', width=8, height=6)


ggplot(scen_df, aes(x=Month, group=1)) +
  facet_grid(Rec_Scen~utilpow, scales='free_y') +
  expand_limits(y=c(0, 0.2)) +
  geom_line(aes(y=SPR)) +
  geom_line(data=rec_pattern_df2, aes(y=Rec, group=1), linetype=2) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
        plot.background=element_rect(color='white'),
        legend.background = element_blank(),
        legend.position = 'bottom') +
  labs(y='SPR')
ggsave('Figures/SimTest/SPRopt_patterns.png', width=8, height=6)

scen_df2 <- scen_df %>% group_by(utilpow, Rec_Scen) %>%
  mutate(ctot=sum(Catch), relCatch=Catch/ctot) %>%
  group_by(Rec_Scen) %>%
  mutate(relC=round(ctot/ctot[utilpow==1],2))

ggplot(scen_df2, aes(x=Month, y=relCatch, group=1)) +
  facet_grid(Rec_Scen~utilpow, scales='free_y') +
  expand_limits(y=c(0, 0.2)) +
  geom_line() +
  geom_text(data=scen_df2 %>% filter(M==1),
            aes(x=1, y=Inf, label=relC), vjust="inward",hjust="inward") +
  geom_line(data=rec_pattern_df2, aes(y=Rec, group=1), linetype=2) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
        plot.background=element_rect(color='white'),
        legend.background = element_blank(),
        legend.position = 'bottom') +
  labs(y='Relative Catch')
ggsave('Figures/SimTest/Catchopt_patterns.png', width=8, height=6)

rec_pattern_df2$utilpow <- NULL
left_join(scen_df, rec_pattern_df2) %>% group_by(Month, Rec_Scen, utilpow) %>%
  summarize(SPR=round(mean(SPR),2),
            F=round(mean(Fishing.Mortality),2),
            C=round(mean(Catch),3),
            Rec=round(mean(Rec),2)) %>%
  filter(Month=='Aug', utilpow==0.4)


left_join(scen_df, rec_pattern_df2) %>% group_by(Rec_Scen, utilpow) %>%
  summarize(SPR=round(mean(SPR),2),
            F=round(mean(Fishing.Mortality),2),
            C=round(mean(Catch),3),
            Rec=round(mean(Rec),2)) %>%
  filter(utilpow==0.8)



# Sensitivity Tests ----

## Incorrect M in Assessment ----
# Diffuse recruitment
# High recsd

# # Simulate with M = 0.1 and estimate with M=0.2
# Pars <- Scens[[7]]
# scen_list <- list()
#
# # loop over simulations
# for (i in 1:nsim) {
#   message('Sim: ', i,'/', nsim)
#   # simulate
#   SimPop <- Simulate(Pars, sim=i)
#
#   # assess
#   sim_data <- list()
#   nts <- length(SimPop$Biomass)
#   ndata_ts <- length(SimPop$Eff_ind)
#   data_ts <- (nts-ndata_ts+1):nts
#   nage <- length(Pars$Weight_Age)
#   sim_data$Weight_Age <- Pars$Weight_Age
#   sim_data$Weight_Age_SD <- Pars$Weight_Age_SD
#   sim_data$M_at_Age <- rep(0.2, nage)
#   sim_data$Mat_at_Age <- Pars$Mat_at_Age
#   sim_data$PSM_at_Age <- Pars$PSM_at_Age
#   sim_data$h  <- Pars$h
#   sim_data$Effort <- SimPop$Eff_ind
#   sim_data$Effort_SD <- rep(0.2, ndata_ts)
#   sim_data$CPUE <- SimPop$Index
#   sim_data$CPUE_SD <- rep(0.2, ndata_ts)
#   sim_data$WghtBins <- Pars$WghtBins
#   sim_data$WghtMids <- Pars$WghtMids
#   sim_data$CAW <- SimPop$CAW_samp
#   sim_data$CAW_ESS <- rep(Pars$CAW_ESS, ndata_ts)
#   sim_data$model <- 'SLAM'
#   sim_data$Fit_Effort <- 1
#   sim_data$Fit_CPUE <- 1
#   sim_data$use_Frwpen <- 1
#   sim_data$use_R0rwpen <- 1
#   sim_data$use_Fmeanprior <- 0
#   sim_data$F_meanprior <- 1
#
#   do_assess <- try(Assess(sim_data), silent=TRUE)
#
#   if (inherits(do_assess, 'try-error')) {
#     # re-run up to 10 times
#     x <- 0
#     while ((inherits(do_assess, 'try-error')) && x < 11) {
#       x <- x+1
#       message(x)
#       do_assess <- try(Assess(sim_data), silent=TRUE)
#     }
#   }
#
#   plot(SimPop$SPR[61:120],type='l')
#   lines(do_assess$rep$SPR, col='blue')
#
#   scen_list[[i]] <- list(Sim=SimPop, Assess=do_assess, sim_data=sim_data)
# }
#
# # Estimate with M = 0.2
#
# # and reverse
#
# grid



# Sensitivity tests:
# - incorrect M in assessment
# - incorrect growth in assessment
# - incorrect growth_sd in assessment
# - senstivity to penalties
# - incorrect maturity in assessment
# - sensitity to incorrect h in optimization

# Case study application

