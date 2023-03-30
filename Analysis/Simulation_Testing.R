# library(pkgdown)
# pkgdown::build_site()

remotes::install_github('blue-matter/SLAM')

# remotes::install_github("mlysy/TMBtools") # this is needed to compile the TMB code??

# To do:
# - test with a new version of R - does it install dependencies?
# - make standard data file structure
# - test and confirm for 12 months of data with Pulse and Constant Recruitment
# - test under perfect conditions with 12 month CAW and CAW & Effort

# add plot(data) function

library(SLAM)
Example_Data()

Data <- Import_Data('C:/Users/User/AppData/Local/R/win-library/4.2/SLAM/Data_Example_Binned.xlsx')

Report(Data)

# set max CAW_ESS
# checks ??
MyAssess <- Assess(Data)



xlfile <- 'inst/Data_Example_Binned.xlsx'
xlfile <- 'inst/Data_Example_Raw.xlsx'

# deal with missing values
# missing CAW, Effort, and Index data?

# add weight units?
# add effort units?

devtools::load_all()

data <- Import_Data('inst/Data_Example_Binned.xlsx')
data <- Import_Data('inst/Data_Example_Raw.xlsx', BinWidth = 1, BinMax=5)

x <- data


Data <- Import_Data()

plot(Data)

MyAssess <- Assess(Data)

Report(Assess)

library(ggplot2)
library(purrr)

library(readxl)
library(xlsx)


class(data)

library(cowplot)
library(ggplot2)







mu <- log(Simulation$LifeHistory$Weight_Age_Mean) -0.5*Simulation$LifeHistory$Weight_Age_SD^2
AWK[,1] <- plnorm(Weight_Bins[2], mu, Simulation$LifeHistory$Weight_Age_SD)


tt <- Import_Data(Sampled_Data)



plot(Data)
Report(Data)


plot(Assess)
Report(Assess)




calc_F_RE <- function(sim, assess, Simulation, n_months=12) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% assess$Data$Month_ind)
  Year <- sort(assess$Data$Year)
  Month <- assess$Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   OM=OM$F_mort,
                   Estimate=assess$rep$F_m)

  df <- df %>% tail(n_months)
  # median((df$Estimate-df$OM)/df$OM)
  (median(df$Estimate)-median(df$OM))/median(df$OM)
}

calc_SB_SB0_RE <- function(sim, assess, Simulation, n_months=12) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% assess$Data$Month_ind)
  Year <- sort(assess$Data$Year)
  Month <- assess$Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   OM=OM$SB_fished/OM$SB_unfished_eq,
                   Estimate=assess$rep$SB_m/assess$rep$SB0_m)

  df <- df %>% tail(n_months)
  # median((df$Estimate-df$OM)/df$OM)

  (median(df$Estimate)-median(df$OM))/median(df$OM)
}



calc_SPR_RE <- function(sim, assess, Simulation, n_months=12) {
  OM <- Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% assess$Data$Month_ind)
  Year <- sort(assess$Data$Year)
  Month <- assess$Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   OM=OM$SPR,
                   Estimate=assess$rep$SPR)

  df <- df %>% tail(n_months)
  # median((df$Estimate-df$OM)/df$OM)
  (median(df$Estimate)-median(df$OM))/median(df$OM)
}



Sim_Test <- function(x, grid, Simulation, Sampling, nsim) {

  # Generate Data
  Sampling$n_recent_months <- grid$n_months[x]
  Sampled_Data <- Generate_Data(Simulation, Sampling = Sampling)

  # Do assessment
  outlist <- list()
  for (i in 1:nsim) {
    message(i, '/', nsim)
    Data <- Import_Data(Sampled_Data, sim=i, Data_types = grid$data_types[x])

    Parameters <- Initialize_Parameters(Data)
    assess <- Assess(Data, Parameters)

    df <- Simulation$Time_Series %>% filter(Sim==i)
    df <- df %>% tail(12)

    outlist[[i]] <- data.frame(Sim=i,
                               F=median(df$F_mort),
                               SPR= median(df$SPR),
                               SB_SB0= median(df$SB_fished/df$SB_unfished_eq),
                               RE_F=calc_F_RE(i, assess, Simulation),
                               RE_SPR=calc_SPR_RE(i, assess, Simulation),
                               RE_SB_SB0=calc_SB_SB0_RE(i, assess, Simulation))
  }

  DF <- do.call('rbind', outlist)
  DF$n_months <-  grid$n_months[x]
  DF$Data_types <- grid$data_types[x]
  DF
}




# ---- Continuous Recruitment ----
nsim <- 200

n_months_vector <- c(12, 24, 36, 60, 120)
data_types_vector <- c('CAW', 'CAW+Effort', 'CAW+Effort+Index')
grid <- expand.grid(n_months=n_months_vector, data_types=data_types_vector,
                    stringsAsFactors = FALSE)

Simulation <- Simulate(LifeHistory, Exploitation, nsim=nsim)

Sampling <- SLAM::Sampling
for (x in 1:nrow(grid)) {
  message('Scenario: ', x, '/', nrow(grid))
  DF <- Sim_Test(x, grid=grid,Simulation, Sampling, nsim=nsim)
  # save
  fl <- paste0(paste(grid[x,], collapse = '_'), '.rda')
  saveRDS(DF, file.path('Analysis/sim_testing/results/continuous', fl))
}

outList <- list()
for (x in 1:nrow(grid)) {
  fl <- paste0(paste(grid[x,], collapse = '_'), '.rda')
  obs <- readRDS(file.path('Analysis/sim_testing/results/continuous', fl))
  outList[[x]] <- obs
}
DF <- purrr::list_rbind(outList)
fl <- 'Continuous.rda'
saveRDS(DF, file.path('Analysis/sim_testing/results', fl))


# ---- Pulse Recruitment ----
Simulation <- Simulate(LifeHistory_Pulse, Exploitation, nsim=nsim)

Sampling <- SLAM::Sampling
for (x in 1:nrow(grid)) {
  message('Scenario: ', x, '/', nrow(grid))
  DF <- Sim_Test(x, grid=grid,Simulation, Sampling, nsim=nsim)
  # save
  fl <- paste0(paste(grid[x,], collapse = '_'), '.rda')
  saveRDS(DF, file.path('Analysis/sim_testing/results/pulse', fl))
}

outList <- list()
for (x in 1:nrow(grid)) {
  fl <- paste0(paste(grid[x,], collapse = '_'), '.rda')
  obs <- readRDS(file.path('Analysis/sim_testing/results/pulse', fl))
  outList[[x]] <- obs
}
DF <- purrr::list_rbind(outList)
fl <- 'Pulse.rda'
saveRDS(DF, file.path('Analysis/sim_testing/results', fl))



# ----- Process Results ----

continuous_DF <- readRDS(file.path('Analysis/sim_testing/results', 'Continuous.rda'))
pulse_DF <- readRDS(file.path('Analysis/sim_testing/results', 'Pulse.rda'))

continuous_DF$Scenario <- 'Continuous'
pulse_DF$Scenario <- 'Pulse'

DF <- bind_rows(continuous_DF,pulse_DF)

DF$n_months <- factor(DF$n_months)
DF <- DF %>% tidyr::pivot_longer(., cols=c(RE_F, RE_SPR, RE_SB_SB0))

ggplot(DF, aes(x=n_months, y=value, fill=Scenario)) +
  facet_grid(name~Data_types, scales='free_y') +
  geom_boxplot()

# Add SB/SB0 and SB/SB_opt
tt <- DF %>% filter(Data_types=='CAW+Effort+Index',
                    Scenario=='Continuous', name=='RE_F', n_months==120)
sd(tt$value)

tt <- DF %>% filter(Data_types=='CAW',
                    Scenario=='Continuous', name=='RE_F', n_months==12)
sd(tt$value)

tt <- DF %>% filter(Data_types=='CAW+Effort+Index',
              Scenario=='Continuous', name=='RE_SB_SB0', n_months==24)

tt %>% filter(value==max(value))

df <- readRDS('Analysis/sim_testing/results/continuous/24_CAW+Effort+Index.rda')
df$RE_SB_SB0 %>% which.max()


