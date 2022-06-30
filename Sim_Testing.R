

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


## Recruitment Process Error Scenarios ----
recsd_scen_names <- 'High'

## Natural Mortality Scenarios ----
M_scen_names <- 0.15 # c('0.1','0.2')

## Build Scenario Grid ----
grid <- expand.grid(rec_scen_names=rec_scen_names,
                    recsd_scen_names=recsd_scen_names,
                    M_scen_names=M_scen_names)


## --- Exploitation History ---
nyears <- 10
nts <- nyears * 12
nsim <- 100
set.seed(1001)

# re-do this - last 5 years matching data
# assumptions for earlier years
# - another run with missing data to match


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
  labs(x="", y='Relative Effort') +
  scale_x_date(date_breaks = "6 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_blank())
ggsave('Figures/SimTest/Effort.png', width=6, height=4)


# ---- Make Parameters List ----
Scens <- list()

for (s in 1:nrow(grid)) {
  Pars <- list()
  Pars$Name <- paste0(apply(grid[s,], 2, as.character), collapse="_")
  Pars$maxage <- 14
  Pars$M <- grid$M_scen_names[s] %>% as.character() %>% as.numeric()
  Pars$Weight_Age <- data$Weight_Age
  Pars$Weight_Age_SD <- data$Weight_Age_SD
  Pars$Mat_at_Age <- data$Mat_at_Age
  Pars$PSM_at_Age <- data$PSM_at_Age
  Pars$h <- data$h
  Pars$WghtBins <-  data$WghtBins
  Pars$WghtMids <-  data$WghtMids

  # recruitment pattern
  ind <- switch(as.character(grid$rec_scen_names[s]), Constant=1, Pulse=2, Diffuse=3, 'Bi-Modal'=4)
  Pars$rec_mu <- reclist[[ind]]
  Pars$Rbar <- 1E5

  # recruitment deviations
  rec_sd <- switch(as.character(grid$recsd_scen_names[s]), Medium=0.6, High=0.9)
  rec_devs <- exp(rnorm(nts*nsim, -0.5*rec_sd^2, rec_sd))
  Pars$rec_devs <- matrix(rec_devs, nrow=nts, ncol=nsim)

  # selectivity
  Pars$sA50 <- 5.5
  Pars$sA95 <- 6.0

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



# ---- Calculate optimal F pattern for Recruitment Scenarios ----

subgrid <- grid %>% distinct(rec_scen_names, M_scen_names)
subgrid$RecSD <- 'Medium'
subgrid <- subgrid %>% select(rec_scen_names, RecSD, M_scen_names)
scenlist <- list()
for (scen in 1:nrow(subgrid)) {
  message('Scenario:', scen)
  sub_grid <- subgrid[scen,]
  Name <- paste0(apply(sub_grid, 2, as.character), collapse="_")
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

  utilpow_vec <- seq(0.2, 1, by=0.2)
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

  df$Rec_Scen <- sub_grid$rec_scen_names
  df$M <- sub_grid$M_scen_names
  scenlist[[scen]] <- df
}
scen_df <- do.call('rbind', scenlist)
scen_df$Month <- factor(scen_df$Month, levels=month.abb, ordered = TRUE)
scen_df$utilpow <- factor(scen_df$utilpow)

scen_df2 <- scen_df %>% filter(M==0.1)
scen_df2 <- scen_df2 %>% tidyr::pivot_longer(cols=3:5)

rec_pattern_df2 <- rec_pattern_df
rec_pattern_df2$name <- 'Catch'
rec_pattern_df2$value <- rec_pattern_df2$Rec
rec_pattern_df2$utilpow <- NA
rec_pattern_df2$Rec_Scen <- rec_pattern_df2$Scenario
ggplot(scen_df2, aes(x=Month, y=value, color=utilpow,
                     group=utilpow)) +
  facet_grid(name~Rec_Scen, scales='free_y') +
  geom_line() +
  geom_line(data=rec_pattern_df2, linetype=2) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_rect(color='white'),
        legend.background = element_blank(),
        legend.position = 'bottom') +
  labs(y='Value', color='Utility Exponent')

ggsave('Figures/SimTest/Fopt_patterns.png', width=8, height=6)


c_df <- scen_df2 %>% filter(name=='Catch') %>% group_by(Rec_Scen, utilpow ) %>%
  summarise(C=sum(value)) %>%
  mutate(Crel=C/C[utilpow==1])

cols <- 1:4

ggplot(c_df, aes(x=utilpow, y=Crel, color=Rec_Scen,
                 group=Rec_Scen)) +
  geom_line(size=1.2) +
  expand_limits(y=0)+
  scale_color_manual(values=cols) +
  theme_clean() +
  theme(plot.background=element_rect(color='white'),
        legend.background = element_blank()) +
  labs(x='Utility Exponent', y='Relative Catch',
       linetype="Recruitment Pattern", color="Recruitment Pattern")

ggsave('Figures/SimTest/RelCatch.png', width=6, height=4)

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
    sim_data$Effort_SD <- rep(0.2, ndata_ts)
    sim_data$CPUE <- SimPop$Index
    sim_data$CPUE_SD <- rep(0.2, ndata_ts)
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

    if (inherits(do_assess, 'try-error')) {
      # re-run up to 10 times
      x <- 0
      while ((inherits(do_assess, 'try-error')) && x < 11) {
        x <- x+1
        message(x)
        do_assess <- try(Assess(sim_data), silent=TRUE)
      }
    }

    scen_list[[i]] <- list(Sim=SimPop, Assess=do_assess, sim_data=sim_data)
  }
  fl <- file.path(paste0('Results/Sim_Test/', Name, '.rda'))
  saveRDS(scen_list, fl)

}

# ---- Summarise Performance ----
scenlist <- list()
for (scen in 1:nrow(grid)) {
  message('Scenario:', scen)
  sub_grid <- grid[scen,]
  Name <- paste0(apply(sub_grid, 2, as.character), collapse="_")
  fl <-  file.path(paste0('Results/Sim_Test/', Name, '.rda'))
  obj <- readRDS(fl)

  Rec_scen <- as.character(sub_grid[1,1])
  RecSD  <- as.character(sub_grid[1,2])
  M  <- as.numeric(as.character(sub_grid[1,3]))


  dflist <- list()
  for (sim in 1:nsim) {
    message('Sim:', sim)
    SimPop <- obj[[sim]]$Sim
    Assess_obj <- obj[[sim]]$Assess

    nts <- length(SimPop$Biomass)
    ndata_ts <- length(SimPop$Eff_ind)
    data_ts <- (nts-ndata_ts+1):nts
    last_year <- (ndata_ts-11):ndata_ts

    dflist[[sim]] <- data.frame(R0_act=SimPop$Rec_Pattern,
                                R0_pred= Assess_obj$rep$R0_m,
                                SPR_act=SimPop$SPR[data_ts][last_year],
                                SPR_pred=Assess_obj$rep$SPR[last_year],
                                F_act=SimPop$F_m[data_ts][last_year],
                                F_pred=Assess_obj$rep$F_m[last_year],
                                Sim=sim,
                                Rec_scen=Rec_scen,
                                RecSD=RecSD,
                                M=M)
  }
  scenlist[[scen]] <- do.call('rbind', dflist)
}

DF <- do.call('rbind', scenlist)
DF$m <- 1:12
DF$Month <- month.abb[DF$m]
DF$Month <- factor(DF$Month, levels=month.abb, ordered=TRUE)
DF$RecSD <- factor(DF$RecSD,levels=c("Medium", 'High'), ordered = TRUE)
DF$Rec_scen <- factor(DF$Rec_scen, levels=unique(grid$rec_scen_names), ordered = TRUE)

DF_true <- DF %>% filter(Sim==1)

for (rec_scen in unique(grid$rec_scen_names)) {
  ggplot(DF %>% filter(Rec_scen==rec_scen),
         aes(x=Month, y=R0_pred, group=Sim)) +
    facet_grid(M~RecSD) +
    geom_point(stat='summary', fun=mean, color='white') +
    stat_summary(fun=mean, geom="line", alpha=0.3) +
    geom_line(aes(x=Month, y=R0_act, group=Sim), color='blue', size=1,
              data=DF_true %>% filter(Rec_scen==rec_scen)) +
    labs(y="Relative Monthly Recruitment") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          plot.background=element_blank())

  ggsave(paste0('Figures/SimTest/R0_', rec_scen, '.png'), width=6, height=4)
}

DF_SPR <- DF %>% group_by(Rec_scen, RecSD, M, Sim) %>%
  summarize(SPR_MRE=median((SPR_pred-SPR_act)/SPR_act))

ggplot(DF_SPR, aes(x=Rec_scen, y=SPR_MRE, fill=as.factor(M))) +
  facet_grid(~RecSD) +
  geom_hline(yintercept = 0, linetype=2) +
  expand_limits(y=c(-1,2)) +
  geom_boxplot() +
  theme_clean() +
  labs(x='Recruitment Scenario', y='Median Relative Error SPR',
       fill='Natural mortality (M)') +
  theme(plot.background=element_blank(),
        legend.background = element_blank(),
        legend.position = 'bottom')

ggsave('Figures/SimTest/SPR_MRE.png', width=6, height=4)


DF_F <- DF %>% group_by(Rec_scen, RecSD, M, Sim) %>%
  summarize(F_MRE=median((F_pred-F_act)/F_act))

ggplot(DF_F, aes(x=Rec_scen, y=F_MRE, fill=as.factor(M))) +
  facet_grid(~RecSD) +
  geom_hline(yintercept = 0, linetype=2) +
  expand_limits(y=c(-1,2)) +
  geom_boxplot() +
  theme_clean() +
  labs(x='Recruitment Scenario', y='Median Relative Error Fishing Mortality (F)',
       fill='Natural mortality (M)') +
  theme(plot.background=element_blank(),
        legend.background = element_blank(),
        legend.position = 'bottom')

ggsave('Figures/SimTest/F_MRE.png', width=6, height=4)





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

