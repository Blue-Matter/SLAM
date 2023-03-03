library(pkgdown)
pkgdown::build_site()

remotes::install_github('blue-matter/SLAM')

# remotes::install_github("mlysy/TMBtools") # this is needed to compile the TMB code??

# To do:
# - make standard data file structure
# - test and confirm for 12 months of data with Pulse and Constant Recruitment
# - test under perfect conditions with 12 month CAW and CAW & Effort

# add plot(data) function


library(SLAM)

# testing under perfect conditions
LifeHistory_No_Error <- SLAM::LifeHistory
LifeHistory_No_Error$sigmaR <- 0.01

Exploitation_No_Error <- SLAM::Exploitation
Exploitation_No_Error$q_cv <- 0.0001

Simulation <- Simulate(LifeHistory=LifeHistory_No_Error,
                       Exploitation=Exploitation_No_Error,
                       nsim=2)

Perfect_Sampling <- SLAM::Sampling
Perfect_Sampling$CPUE_CV <- 0.0001
Perfect_Sampling$Catch_CV <- 0.0001
Perfect_Sampling$Effort_CV <- 0.0001
Perfect_Sampling$CAW_Annual_ESS <- 1E6
Perfect_Sampling$CAW_Annual_Sample_Size <- 1E6
Perfect_Sampling$CAA_Annual_ESS <- 1E6
Perfect_Sampling$CAA_Annual_Sample_Size <- 1E6
Perfect_Sampling$n_recent_months <- 360

Sampled_Data <- Generate_Data(Simulation, Perfect_Sampling)


Data <- Import_Data(Sampled_Data)

# Parameters <- Initialize_Parameters(Data)

Parameters <- Initialize_Parameters_OM(Simulation)

# ----- Assessment ----
control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20)

map=list(log_sigmaF_m=factor(NA),
         log_sigmaR=factor(NA),
         log_sigmaR0=factor(NA))


nmonths <- length(Data$Month_ind)
if (nmonths<24) {
  map$logRec_Devs <- rep(factor(NA), length(Parameters$logRec_Devs))
}

if (!is.null(map$log_sigmaR)) {
  Random <- NULL
} else {
  Random <- 'logRec_Devs'
}

Year <- Data$Year
Month <- Data$Month

TMB_Data <- Data
TMB_Data$Year <- TMB_Data$Month <- NULL

TMB_Data$Fit_CAW <- 0
TMB_Data$Fit_CPUE <- 0

TMB_Data$use_Frwpen <- 0
TMB_Data$use_R0rwpen <- 0


obj <- TMB::MakeADFun(data=TMB_Data, parameters=Parameters, DLL="SLAM_TMBExports",
                      silent=TRUE, hessian=FALSE, map=map, random=Random)

starts <- obj$par
opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control)),silent=TRUE)

rep <- obj$report()

Est <- rep
Est$nll_joint

plot_Effort_fit(Data, Est)

plot_CPUE_fit(Data, Est)
plot_CAA_fit(Data, Est)
plot_CAW_fit(Data, Est)

plot_F(Simulation, Est)
plot_SB(Simulation, Est)
plot_SPR(Simulation, Est)





do_assess <- Assess(Data, Parameters)
do_assess$opt

do_assess$sdreport

plot(OM_DF$Month_ind, OM_DF$F_mort, type='l', ylim=c(0,max(c(OM_DF$F_mort, Assess_DF$F_mort))))
lines(Assess_DF$Month_ind, Assess_DF$F_mort, col='blue')

median((Assess_DF$F_mort-OM_DF$F_mort)/OM_DF$F_mort)
median((Assess_DF$SPR-OM_DF$SPR)/OM_DF$SPR)

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
    theme_bw()
}


plot_F <- function(Simulation, Est) {
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
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    theme_bw()
}

plot_SB <- function(Simulation, Est) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% Data$Month_ind)


  Year <- sort(Data$Year)
  Month <- Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   Simulated=OM$SB_fished/OM$SB_unfished_eq,
                   Predicted=Est$SB_m/Est$SB0_m)
  df$Year_Month <- factor(df$Year_Month, ordered = TRUE, levels=unique(df$Year_Month))
  ggplot(df, aes(x=Year_Month, group=1)) +
    geom_line(aes(y=Simulated), color='blue') +
    geom_line(aes(y=Predicted), linetype=2) +
    theme_bw()
}

plot_SPR <- function(Simulation, Est) {
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
    theme_bw()
}

































library(SLAM)
library(dplyr)

nsim <- 100

# Don't include CAA data - not available in these fisheries
Scenarios <- Scenario_Grid %>% filter(grepl('CAA', Data_types) == FALSE,
                                      Monthly_Recruitment_Pattern =='Constant',
                                      Conditions=='Idealized')

# Simulate data
SimScenario <- Scenarios %>% filter(Data_n_months==max(Data_n_months)) %>%
  group_by(Conditions) %>%
  filter(Data_n_months==max(Data_n_months), Data_types=='CAW')

Parameters <- get(SimScenario$Name)
LifeHistory <- Parameters$LifeHistory
Exploitation <- Parameters$Exploitation
Data <- Parameters$Data

SimMod <- Simulate(LifeHistory, Exploitation, Data, nsim = nsim)

# Modify simulated data time-series based on scenario
Scen_res <- list()
for (i in 1:nrow(Scenarios)) {
  message('Scenario: ', i, '/', nrow(Scenarios))
  Scenario <- Scenarios[i,]
  Data_types <- Scenario$Data_types
  Data_n_months <- Scenario$Data_n_months

  modify_data <- function(SimMod, Data_n_months) {
    Month_ind <- SimMod$Data_TS_DF$Month_ind %>% unique()
    Month_ind_keep <- Month_ind[(length(Month_ind)-Data_n_months+1):length(Month_ind)]

    SimMod$Data_TS_DF <- SimMod$Data_TS_DF %>% filter(Month_ind %in% Month_ind_keep)

    SimMod$Data_TS_DF <- SimMod$Data_TS_DF %>% group_by(Sim) %>%
      mutate(CPUE=CPUE/mean(CPUE), Effort=Effort/mean(Effort)) %>% ungroup()

    SimMod$Data_CAW_DF <- SimMod$Data_CAW_DF %>% filter(Month_ind %in% Month_ind_keep)
    SimMod$Data_CAA_DF <- SimMod$Data_CAA_DF %>% filter(Month_ind %in% Month_ind_keep)
    SimMod
  }

  SimMod_Scen <- modify_data(SimMod, Data_n_months)

  assess_list <- list()
  for (x in 1:nsim) {
    message(x, '/', nsim)
    data <- Construct_Data_OM(x, SimMod_Scen, Data_types=Data_types)
    Assessment <- Do_Assess(data)
    assess_list[[x]] <- Compare_OM_Assess(x, SimMod_Scen, Assessment)
  }

  assess_DF <- do.call('rbind', assess_list)

  # Last 12 months
  n_months <- assess_DF$Month_ind %>% max()
  sel_months <- (n_months-11):n_months
  assess_DF <- assess_DF %>% filter(Month_ind %in% sel_months)
  assess_DF <- bind_cols(assess_DF, Scenario)

  saveRDS(assess_DF, file.path('Results/Sim_Testing', paste0(Scenario$Name, '.rda')))

  Scen_res[[i]] <- assess_DF
}
RE_Results <- do.call('rbind', Scen_res)

library(ggplot2)

DF$RE_F <- NA
DF$RE_SPR <- NA
DF <- RE_Results %>%
  group_by(Sim, Data_types , Data_n_months, Monthly_Recruitment_Pattern, Conditions) %>%
  summarise(RE_F=(F_mort[Model=='Assessment']-F_mort[Model=='OM'])/F_mort[Model=='OM'],
         RE_SPR=(SPR[Model=='Assessment']-SPR[Model=='OM'])/SPR[Model=='OM'])

DF$Data_n_months <- factor(DF$Data_n_months)
DF$SPR <- DF$RE_SPR
DF$F <- DF$RE_F

DF_long <- DF %>% tidyr::pivot_longer(., cols=c('SPR', 'F'))

ggplot(DF_long %>% filter(Monthly_Recruitment_Pattern=='Constant'), aes(x=Data_n_months, y=value)) +
  facet_grid(name~Data_types, scales='free') +
  geom_boxplot(aes(fill=Conditions)) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  labs(y='Relative Error', x="n Months")


ggplot(DF_long %>% filter(Monthly_Recruitment_Pattern=='Pulse'), aes(x=Data_n_months, y=value)) +
  facet_grid(name~Data_types, scales='free') +
  geom_boxplot(aes(fill=Conditions)) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  labs(y='Relative Error', x="n Months")


OM_DF <- OM_Assess_DF %>% filter(Model=='OM')
Assess_DF <- OM_Assess_DF %>% filter(Model=='Assessment')







plot(OM_DF$Month_ind, OM_DF$F_mort, type='l', ylim=c(0,max(c(OM_DF$F_mort, Assess_DF$F_mort))))
lines(Assess_DF$Month_ind, Assess_DF$F_mort, col='blue')

median((Assess_DF$F_mort-OM_DF$F_mort)/OM_DF$F_mort)
median((Assess_DF$SPR-OM_DF$SPR)/OM_DF$SPR)



rep <- Assessment$rep
rep$F_minit
rep$F_m
rep$logRec_Devs %>% exp()
OM_DF$F_mort

par(mfrow=c(3,4))
dd <- dim(data$CAA)
ms <- ((dd[2] - 12 + 1):dd[2])
for (i in ms) {
  plot(0:14, data$CAA[,i]/sum(data$CAA[,i]), type='l')
  lines(0:14, rep$predCAA[,i], col='blue')
}

par(mfrow=c(2,2))
plot(data$CPUE, type='l')
lines(rep$stpredIndex, col='blue')

plot(data$Effort, type='l')
lines(rep$StEffort, col='blue')

Mod <- SimMod$OM_DF %>% filter(Sim==x)
month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
Mod_data <- SimMod$OM_DF %>% filter(Sim==x, Month_ind%in%month_ind)

plot(Mod_data$F_mort, type='l', ylim=c(0, max(c(Mod_data$F_mort, rep$F_m))))
lines(rep$F_m, col='blue')
nonzero <- which(Mod_data$F_mort>0)
median((rep$F_m[nonzero]-Mod_data$F_mort[nonzero])/Mod_data$F_mort[nonzero])

rep$F_minit
rep$F_m

par(mfrow=c(3,4))
dd <- dim(data$CAA)
ms <- ((dd[2] - 12 + 1):dd[2])
for (i in ms) {
  plot(0:14, data$CAA[,i]/sum(data$CAA[,i]), type='l')
  lines(0:14, rep$predCAA[,i], col='blue')
}

rep$nll_joint
rep$CAAnll

# Run test for no obs error
# Run test with process error

# set default parameters
# TODO - deal with missing data !!

library(SLAM)
library(dplyr)
library(ggplot2)

nsim <- 100

Scenarios <- Scenario_Grid %>% filter(Conditions=='Idealized',
                                      Monthly_Recruitment_Pattern=='Pulse')



Scenarios <- Scenario_Grid %>% filter(Conditions=='Idealized',
                                      Monthly_Recruitment_Pattern=='Constant')


Scenarios <- Scenario_Grid %>% filter(Conditions=="Process+Observation Error",
                                      Monthly_Recruitment_Pattern=='Pulse')

Scenarios <- Scenario_Grid

set_data_types <- function(data, Data_types) {
  Data_types <- strsplit(Data_types, "\\+")[[1]]
  if (!'CAW' %in% Data_types) data$Fit_CAW <-0
  if (!'Index' %in% Data_types) data$Fit_CPUE <-0
  if (!'Effort' %in% Data_types) data$Fit_Effort <-0
  data
}

Scen_res <- list()
for (i in 1:nrow(Scenarios)) {
# for (i in 10:nrow(Scenarios)) {
  message('Scenario ' , i, '/', nrow(Scenarios))
  Scenario <- Scenarios[i,]
  Parameters <- get(Scenario$Name)
  LifeHistory <- Parameters$LifeHistory
  Exploitation <- Parameters$Exploitation
  Data <- Parameters$Data

  SimMod <- Simulate(LifeHistory, Exploitation, Data, nsim = nsim)

  assess_list <- list()
  for (x in 1:nsim) {
    message(x, '/', nsim)
    data <- Construct_Data_OM(x, SimMod)
    data <- set_data_types(data, Data_types=Scenario$Data_types)

    data$use_F_seasonrwpen <- 0
    data$use_Frwpen <- 1
    data$use_R0rwpen <- 1

    parameters <- Initialize_Parameters(data,
                                        as50=4, as95=6,
                                        Feq_init=0.05,
                                        sigmaR=0.5,
                                        sigmaF_m=0.8,
                                        sigmaF_season=0.9,
                                        sigmaR0=0.5)

    map=list(log_sigmaR=factor(NA),
             log_sigmaR0=factor(NA),
             log_sigmaF_m=factor(NA),
             log_sigmaF_season=factor(NA))

    # Random effects
    Random <- NULL

    if (is.null(map$log_sigmaR))
      Random <- c(Random, 'logRec_Devs')

    n_years <- floor(ncol(data$CAW)/12)

    if (n_years==1) {
      map$logRec_Devs <- rep(factor(NA), length(parameters$logRec_Devs))
    }

    control=list(eval.max=2E4, iter.max=2E4, abs.tol=1E-20)
    do_opt <- opt_TMB_model(data, parameters, map, Random, control, restarts=10)

    # ----------------------------------------------------------------------- #
    Assessment <- do_opt
    OM_Assess_DF <- Compare_OM_Assess(x, SimMod, Assessment)

    # Plot OM vs Estimates
    OM_DF <- OM_Assess_DF %>% filter(Model=='OM')
    Assess_DF <- OM_Assess_DF %>% filter(Model=='Assessment')

    plot(OM_DF$Month_ind, OM_DF$F_mort, type='l', ylim=c(0,max(OM_DF$F_mort)))
    lines(Assess_DF$Month_ind, Assess_DF$F_mort, col='blue')

    median((Assess_DF$F_mort-OM_DF$F_mort)/OM_DF$F_mort)


    rep <- Assessment$rep
    par(mfrow=c(2,2))
    plot(data$CPUE, type='l')
    lines(rep$stpredIndex, col='blue')

    plot(data$Effort, type='l')
    lines(rep$StEffort, col='blue')

    Mod <- SimMod$OM_DF %>% filter(Sim==x)
    month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
    Mod_data <- SimMod$OM_DF %>% filter(Sim==x, Month_ind%in%month_ind)

    plot(Mod_data$F_mort, type='l', ylim=c(0, max(Mod_data$F_mort)))
    lines(rep$F_m, col='blue')
    nonzero <- which(Mod_data$F_mort>0)
    median((rep$F_m[nonzero]-Mod_data$F_mort[nonzero])/Mod_data$F_mort[nonzero])

    rep$F_minit
    rep$F_m


    par(mfrow=c(3,4))
    dd <- dim(data$CAW)
    ms <- ((dd[2] - 12 + 1):dd[2])
    for (i in ms) {
      plot(data$WghtMids, data$CAW[,i]/sum(data$CAW[,i]), type='l')
      lines(data$WghtMids, rep$predCAW[,i], col='blue')
    }

    # ----------------------------------------------------------------------- #



    assess_list[[x]] <- Compare_OM_Assess(x, SimMod, do_opt)

    # save to disk
    out <- list()
    out$do_opt <- do_opt
    out$assess_DF <- assess_list[[x]]
    name <- paste0(i, '_', Scenario$Name, '.rda')
    out.file <- file.path('results', 'sim_testing', name)
    saveRDS(out, file=out.file)
  }

  assess_DF <- do.call('rbind', assess_list)



  # Calculate relative error - for last 12 months
  n_months <- assess_DF$Month_ind %>% max()
  sel_months <- (n_months-11):n_months
  RE_DF <- assess_DF %>% filter(Month_ind %in% sel_months) %>%
    group_by(Sim) %>%
    summarize(RE_F=median((F_mort[Model=='Assessment']-F_mort[Model=='OM'])/F_mort[Model=='OM']),
              RE_SPR=median((SPR[Model=='Assessment']-SPR[Model=='OM'])/SPR[Model=='OM']),)
  RE_DF <- bind_cols(RE_DF, Scenario)
  Scen_res[[i]] <- RE_DF
}

# ---- Process Results ----
res.dir <- 'C:/Users/User/Documents/GitHub/SLAM/results/sim_testing'
fls <- list.files(res.dir)

DF.list <- list()
for (i in 1:nrow(Scenario_Grid)) {
  Scenario <- Scenario_Grid[i,]

  fl <- paste0(paste(i, Scenario$Name, sep='_'), '.rda')
  assess_DF <- readRDS(file.path(res.dir, fl))$assess_DF

  RE_DF <- assess_DF %>% group_by(Sim) %>%
    summarize(RE_F=median((F_mort[Model=='Assessment']-F_mort[Model=='OM'])/F_mort[Model=='OM']),
              RE_SPR=median((SPR[Model=='Assessment']-SPR[Model=='OM'])/SPR[Model=='OM']),)
  RE_DF <- bind_cols(RE_DF, Scenario)

  DF.list[[i]] <- RE_DF
}
DF <- do.call('rbind', DF.list)

DF$Data_n_months <- factor(DF$Data_n_months)
DF$SPR <- DF$RE_SPR
DF$F <- DF$RE_F

DF_long <- DF %>% tidyr::pivot_longer(., cols=c('SPR', 'F'))

ggplot(DF_long %>% filter(Monthly_Recruitment_Pattern=='Constant'), aes(x=Data_n_months, y=value)) +
  facet_grid(name~Data_types, scales='free') +
  geom_boxplot(aes(fill=Conditions)) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  labs(y='Relative Error', x="n Months")


ggplot(DF_long %>% filter(Monthly_Recruitment_Pattern=='Pulse'), aes(x=Data_n_months, y=value)) +
  facet_grid(name~Data_types, scales='free') +
  geom_boxplot(aes(fill=Conditions)) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  labs(y='Relative Error', x="n Months")



DF_long %>% filter(Monthly_Recruitment_Pattern=='Pulse', Data_n_months==12)





data$use_Eff_rwpen <- 0
data$use_R0rwpen <- 0
data$Fit_CPUE <- 0
data$Fit_Effort <- 0

do_opt <- opt_TMB_model(data, parameters, map, Random, control, restarts=10)


Assessment <- do_opt
OM_Assess_DF <- Compare_OM_Assess(x, SimMod, Assessment)

# Plot OM vs Estimates
OM_DF <- OM_Assess_DF %>% filter(Model=='OM')
Assess_DF <- OM_Assess_DF %>% filter(Model=='Assessment')

plot(OM_DF$Month_ind, OM_DF$F_mort, type='l', ylim=c(0,max(OM_DF$F_mort)))
lines(Assess_DF$Month_ind, Assess_DF$F_mort, col='blue')

median((Assess_DF$F_mort-OM_DF$F_mort)/OM_DF$F_mort)


rep <- Assessment$rep
par(mfrow=c(2,2))
plot(data$CPUE, type='l')
lines(rep$stpredIndex, col='blue')

plot(data$Effort, type='l')
lines(rep$StEffort, col='blue')

Mod <- SimMod$OM_DF %>% filter(Sim==x)
month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
Mod_data <- SimMod$OM_DF %>% filter(Sim==x, Month_ind%in%month_ind)

plot(Mod_data$F_mort, type='l', ylim=c(0, max(Mod_data$F_mort)))
lines(rep$F_m, col='blue')
nonzero <- which(Mod_data$F_mort>0)
median((rep$F_m[nonzero]-Mod_data$F_mort[nonzero])/Mod_data$F_mort[nonzero])


par(mfrow=c(3,4))
dd <- dim(data$CAW)
ms <- ((dd[2] - 12 + 1):dd[2])
for (i in ms) {
  plot(data$WghtMids, data$CAW[,i]/sum(data$CAW[,i]), type='l')
  lines(data$WghtMids, rep$predCAW[,i], col='blue')
}

do_opt$rep$F_minit
do_opt$rep$F_m

do_opt$rep$Za_init




Run_Sim_Test <- function(Scenario, nsim=50) {

  Parameters <- get(Scenario$Name)
  LifeHistory <- Parameters$LifeHistory
  Exploitation <- Parameters$Exploitation
  Exploitation$Effort_Annual_SD <- rep(0.1, length(Exploitation$Effort_Annual_SD))
  Data <- Parameters$Data
  Data$Weight_Bin_Max <- 4

  SimMod <- Simulate(LifeHistory, Exploitation, Data, nsim = nsim)

  assess_list <- list()
  for (x in 1:nsim) {
    assess_list[[x]] <- Run_Assess(x, SimMod, data_types=Scenario$Data_types)
  }
  assess_DF <- do.call('rbind', assess_list)


  # Calculate relative error
  RE_DF <- assess_DF %>% group_by(Sim) %>%
    summarize(RE_F=median((F_mort[Model=='Assessment']-F_mort[Model=='OM'])/F_mort[Model=='OM']),
              RE_SPR=median((SPR[Model=='Assessment']-SPR[Model=='OM'])/SPR[Model=='OM']),)
  RE_DF <- bind_cols(RE_DF, Scenario)

}

set_data_types <- function(data, Data_types) {
  Data_types <- strsplit(Data_types, "\\+")[[1]]
  if (!'CAW' %in% Data_types) data$Fit_CAW <-0
  if (!'Index' %in% Data_types) data$Fit_CPUE <-0
  if (!'Effort' %in% Data_types) data$Fit_Effort <-0
  data
}





Run_Assess <- function(x, SimMod, data_types, map=list(log_sigmaF=factor(NA),
                                                       log_sigmaR=factor(NA),
                                                       log_sigmaR0=factor(NA))) {

  data <- Construct_Data_OM(SimMod, sim=x)
  data <- set_data_types(data, Data_types=data_types)

  Assessment <- Do_Assess(data, map=map, log_sigmaF=log(0.5))

  OM_Assess_DF <- Compare_OM_Assess(x, SimMod, Assessment)
  OM_Assess_DF
}

Scen_res <- list()
for (i in 1:nrow(Scenarios)) {
  message('Scenario ', i, ' of ', nrow(Scenarios))
  Scenario <- Scenarios[i,]
  Scen_res[[i]] <- Run_Sim_Test(Scenario, nsim)
}

DF <- do.call('rbind', Scen_res)
DF$Data_n_months <- factor(DF$Data_n_months)

DF <- DF %>% filter(Data_types!='Index+Effort')
ggplot(DF, aes(x=Data_n_months, y=RE_F)) +
  facet_grid(~Data_types) +
  geom_boxplot()







# Plot OM vs Estimates
OM_DF <- OM_Assess_DF %>% filter(Model=='OM')
Assess_DF <- OM_Assess_DF %>% filter(Model=='Assessment')

plot(OM_DF$Month_ind, OM_DF$F_mort, type='l', ylim=c(0,max(OM_DF$F_mort)))
lines(Assess_DF$Month_ind, Assess_DF$F_mort, col='blue')

# plot(OM_DF$Month_ind, OM_DF$SPR, type='l', ylim=c(0,1))
# lines(Assess_DF$Month_ind, Assess_DF$SPR, col='blue')

median((Assess_DF$F_mort-OM_DF$F_mort)/OM_DF$F_mort)
# mean((Assess_DF$SPR-OM_DF$SPR)/OM_DF$SPR)

rep <- Assessment$rep

par(mfrow=c(2,2))
plot(data$CPUE, type='l')
lines(rep$stpredIndex, col='blue')

plot(data$Effort, type='l')
lines(rep$StEffort, col='blue')

Mod <- SimMod$OM_DF %>% filter(Sim==x)
month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
Mod_data <- SimMod$OM_DF %>% filter(Sim==x, Month_ind%in%month_ind)

plot(Mod_data$F_mort, type='l', ylim=c(0, max(Mod_data$F_mort)))
lines(rep$F_m, col='blue')
nonzero <- which(Mod_data$F_mort>0)
median((rep$F_m[nonzero]-Mod_data$F_mort[nonzero])/Mod_data$F_mort[nonzero])



OM_DF <- OM_Assess_DF %>% filter(Model=='OM')
Assess_DF <- OM_Assess_DF %>% filter(Model=='Assessment')

plot(OM_DF$Month_ind, OM_DF$F_mort, type='l', ylim=c(0,max(OM_DF$F_mort)))
lines(Assess_DF$Month_ind, Assess_DF$F_mort, col='blue')

rep <- Assessment$rep


par(mfrow=c(3,4))
dd <- dim(data$CAW)
ms <- ((dd[2] - 12 + 1):dd[2])
for (i in ms) {
  plot(data$WghtMids, data$CAW[,i]/sum(data$CAW[,i]), type='l')
  lines(data$WghtMids, rep$predCAW[,i], col='blue')
}




# ---- Only CAW Data - Constant  ----


Scenarios <- Scenario_Grid %>% filter(Conditions=='Process+Observation Error',
                                      Monthly_Recruitment_Pattern=='Constant',
                                      Data_types=='CAW')

for (i in 1:nrow(Scenarios)) {
  message('Scenario ', i, ' of ', nrow(Scenarios))
  Scenario <- Scenarios[i,]

  Parameters <- get(Scenario$Name)
  LifeHistory <- Parameters$LifeHistory
  Exploitation <- Parameters$Exploitation
  Data <- Parameters$Data
  SimMod <- Simulate(LifeHistory, Exploitation, Data, nsim = nsim)



  for (x in 1:nsim) {
    message(x, '/', nsim)
    data <- Construct_Data_OM(SimMod, sim=x)
    data <- set_data_types(data, Data_types)

  }

}












control <- list(eval.max=1e4, iter.max=1e4,
                step.min=1, step.max=1,
                trace=0, abs.tol=1e-20)

starts <- obj$par
opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control)),
           silent=TRUE)

rep <- obj$report(obj$env$last.par.best)

Mod <- SimMod$OM_DF %>% filter(Sim==sim)
month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
Mod_data <- SimMod$OM_DF %>% filter(Sim==sim, Month_ind%in%month_ind)

par(mfrow=c(2,2))
plot(data$CPUE, type='l')
lines(rep$stpredIndex, col='blue')

plot(data$Effort, type='l')
lines(rep$StEffort, col='blue')

plot(Mod_data$F_mort, type='l', ylim=c(0, max(Mod_data$F_mort)))
lines(rep$F_m, col='blue')
nonzero <- which(Mod_data$F_mort>0)
mean(rep$F_m[nonzero]/Mod_data$F_mort[nonzero])


plot(rep$R0_m, type='l', ylim=c(0, max(rep$R0_m)))

plot(Mod_data$SPR, type='l', ylim=c(0,1))
lines(rep$SPR, col='blue')

n_months<- ncol(data$CAW)
par(mfrow=c(3,4))
for (i in 1:12) {
  plot(data$WghtMids, data$CAW[,i]/sum(data$CAW[,i]), type='l')
  lines(data$WghtMids, rep$predCAW[,i], col='blue')
}






tt <- SimMod$Data_CAW_DF %>% filter(Sim==1)
plot(tt$Weight, tt$Count, type='l')


mean(rep$F_m[nonzero]/Mod_data$F_mort[nonzero])


Fests <- data.frame(OM=rep(NA,50))
Fests$Ests <- NA



for (sim in 1:50) {
  message(sim, '/50')
  # Fit Assessment Model
  data <- Construct_Data_OM(SimMod, sim=sim)
  data$CPUE_SD <- rep(0.01, length(data$CPUE_SD))
  data$Effort_SD <- rep(0.01, length(data$CPUE_SD))
  parameters <- Initialize_Parameters_OM(SimMod)


  # Do Assessment
  map=list(log_sigmaF=factor(NA),
           log_sigmaR=factor(NA),
           log_sigmaR0=factor(NA))

  if (!is.null(map$log_sigmaR)) {
    Random <- NULL
  } else {
    Random <- 'logRec_Devs'
  }

  map_all_parameters <- function(parameters) {
    list_out <- parameters
    for (x in seq_along(list_out)) {
      list_out[[x]] <- factor(rep(NA, length(list_out[[x]])))
    }
    list_out
  }

  # map <- map_all_parameters(parameters)
  # map$ls50 <- map$lsdelta <- NULL
  #

  obj <- TMB::MakeADFun(data=data, parameters=parameters, DLL="SLAM_TMBExports",
                        silent=TRUE, map=map, hessian=FALSE, random = Random)

  step.min <- 1
  step.max <- 1
  control <- list(eval.max=1e4, iter.max=1e4,
                  step.min=step.min, step.max=step.max,
                  trace=0, abs.tol=1e-20)

  starts <- obj$par
  opt <- try(suppressWarnings(nlminb(starts, obj$fn, obj$gr, control = control)),
             silent=TRUE)

  rep <- obj$report(obj$env$last.par.best)

  Mod <- SimMod$OM_DF %>% filter(Sim==sim)
  month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
  Mod_data <- SimMod$OM_DF %>% filter(Sim==sim, Month_ind%in%month_ind)

  Fests$OM[sim] <- mean(Mod_data$F_mort)
  Fests$Ests[sim] <- mean(rep$F_m)

}

plot(Fests$OM, type='l')
lines(Fests$Ests, col='blue')
median((Fests$Ests-Fests$OM)/Fests$OM)



rep <- obj$report(obj$env$last.par.best)


Mod <- SimMod$OM_DF %>% filter(Sim==sim)
month_ind <- rev(seq(Exploitation$nyears*12, by=-1, length.out=Data$n_recent_months))
Mod_data <- SimMod$OM_DF %>% filter(Sim==1, Month_ind%in%month_ind)

par(mfrow=c(1,2))
plot(data$CPUE, type='l')
lines(rep$stpredIndex, col='blue')

plot(data$Effort, type='l')
lines(rep$StEffort, col='blue')


plot(Mod_data$F_mort, type='l', ylim=c(0, 0.5))
lines(rep$F_m, col='blue')
mean(rep$F_m/Mod_data$F_mort)









Zs <- LifeHistory$M_at_Age+ rep$F_minit
PSM <- LifeHistory$Post_Spawning_Mortality
ind <- 2:15
ind2 <- ind-1
tt <- c(rep$N_m[1,1], rep$N_unfished[ind2,12]* exp(-Zs[ind2]) * (1-PSM[ind2]))
plot(tt, type='l')
lines(rep$N_m[,1], col='blue')

rep$N_unfished[1,12]* exp(-Zs[1]) * (1-PSM[1])
rep$N_m[2,1]


Zs <- -log(N_Age_fished[1,2:15,241]/N_Age_unfished_eq[1,1:14,240])
Zs <- c(Zs,Inf)
Zs <- Zs/(1-LifeHistory$Post_Spawning_Mortality)
Finit <- Zs -  LifeHistory$M_at_Age
Finit[Finit==Inf] <- 5
logF_minit <- log(Finit)


N_Age_unfished[1,1:14,240] * exp(-Zs)
N_Age_fished[1,2:15,241]

TS_OM <- get_TS_OM(SimMod)
Est_DF <- get_TS_Predict(rep)

plot(TS_OM$N_fished, type='l')
lines(Est_DF$N_fished, col='blue')


plot(TS_OM$SB_fished, type='l')
lines(Est_DF$SB_fished, col='blue')




get_TS_OM <- function(SimMod, sim=1) {
  OM_DF <-   Data_Y_M <- SimMod$Data_TS_DF %>% filter(Sim==1) %>% distinct(Year, Month)
  OM_sub <-SimMod$OM_DF %>% filter(Sim==1, Year%in%Data_Y_M$Year, Month%in%Data_Y_M$Month)
  OM_sub$Model <- 'OM'
  OM_sub$Year_Month <- paste(OM_sub$Year, OM_sub$Month, sep="_")
  OM_sub
}


get_TS_Predict <- function(rep) {
  nMonths <- length(obj$env$data$CPUE)
  nYears <- nMonths/12
  currentYr <- obj$env$data$currentYr
  Years <- (currentYr-nYears+1):currentYr
  Years <- rep(Years, each=12)
  Months <- rep(month.abb, nYears)

  Est_TS_DF <- data.frame(Year=Years,
                          Month=Months,
                          Month_ind=1:nMonths,
                          N_fished=apply(rep$N_m,2,sum),
                          SB_fished=rep$SB_m,
                          Catch=rep$predCB,
                          SPR=rep$SPR,
                          Effort=rep$StEffort,
                          F_mort=rep$F_m,
                          Rec_Devs=exp(rep$logRec_Devs),
                          CPUE=rep$stpredCPUE,
                          Model='Assessment'
  )
  Est_TS_DF$Year_Month <- paste(Est_TS_DF$Year, Est_TS_DF$Month, sep="_")
  Est_TS_DF

}






plot(TS_OM$F_mort, type='l', ylim=c(0, 0.1))
lines(Est_DF$F_mort, col='blue')


par(mfrow=c(1,2))
plot(data$CPUE, type='l')
lines(Est_DF$CPUE, col='blue')


plot(data$Effort, type='l')
lines(Est_DF$Effort, col='blue')

plot(TS_OM$SPR, type='l', ylim=c(0,1))
lines(Est_DF$SPR, col='blue')














plot_fishery_dynamics <- function(SimMod, sim=1) {

  df <- data.frame(Month=month.abb,
                   Rel_Recruitment=SimMod$LifeHistory$R0_m,
                   Rel_Effort=SimMod$Exploitation$Effort_Month_Mean) %>%
    tidyr::pivot_longer(., cols=2:3)
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)

  ggplot(df, aes(x=Month, y=value, color=name, group=name)) +
    geom_line() +
    expand_limits(y=0)

  # plot effort
  OM_DF <- SimMod$OM_DF %>% filter(Sim==sim)
  plot(OM_DF$Effort, type='l')


}






# Evaluate Fits





plot(rep$R0_m, type='l', ylim=c(0, 1/12))

# Evaluate Predictions


DF <- bind_rows(TS_OM, Est_OM)

ggplot(DF, aes(x=Year_Month,
               y=SPR,
               group=Model,
               color=Model)) +
  geom_line()


data$CAW[is.na(data$CAW)] <- 0
data$Effort[data$Effort==0] <- 0.001
data$CPUE[data$CPUE==0] <- 0.001
parameters$logF_m[parameters$logF_m==-Inf] <- log(0.001)








# Make Estimates data.frame





Data_TS <- MyMod$Data_TS_DF %>% filter(Sim==sim)


plot(Data_TS$Month_ind, Data_TS$CPUE, type='l')
lines(Est_TS_DF$Month_ind, Est_TS_DF$CPUE, col='blue')

plot(Data_TS$Month_ind, Data_TS$Effort, type='l')
lines(Est_TS_DF$Month_ind, Est_TS_DF$Effort, col='blue')


plot(Data_TS$Month_ind, mod$SPR[301:360], type='l')
lines(Est_TS_DF$Month_ind, Est_TS_DF$SPR, col='blue')

# Compare Estimates with OM
plot(rep$R0_m, type='l')

head(mod)



plot(rep$selA, type='l')


names(rep)
rep$F_m %>% length()
MyMod$OM_DF %>% filter(Sim==sim)

rep$SPR


sdreport <- TMB::sdreport(obj, obj$env$last.par.best)

sdreport$pdHess

rep$F_m
rep$F_month_NLL











plot(mod$Month_ind, mod$Effort, type='l')

plot(mod$Month_ind, mod$Catch, type='l')
plot(mod$Month_ind, mod$SB_fished/mod$SB_unfished, type='l', ylim=c(0,1))

lines(mod$Month_ind, mod$SPR, type='l', col='blue')

plot(mod$SB_fished/mod$SB_unfished, mod$Recruits, type='l', ylim=c(0,max(mod$Recruits)))
plot( mod$Recruits, type='l', ylim=c(0,max(mod$Recruits)))

