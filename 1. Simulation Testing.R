remotes::install_github('blue-matter/SLAM')



# TO DO
# - update methodology
# - write up simulation testing methods and results

library(SLAM)
library(ggplot2)

nsim <- 100

# Simulate historical fishery ----

Simulation <- Simulate(LifeHistory, Exploitation, nsim=nsim)
Simulation_Pulse <- Simulate(LifeHistory_Pulse, Exploitation, nsim=nsim)


# Recruitment Scenario Plot ----
df1 <- data.frame(Month=month.abb, Rec=Simulation$LifeHistory$R0_m, Scenario='Constant')
df2 <- data.frame(Month=month.abb, Rec=Simulation_Pulse$LifeHistory$R0_m, Scenario='Pulse')

df <- bind_rows(df1,df2)
df$Month <- factor(df$Month, levels=month.abb, ordered = TRUE)
ggplot(df, aes(x=Month, y=Rec, color=Scenario, group=Scenario)) +
  geom_line() +
  expand_limits(y=0) +
  theme_bw() +
  labs(y='Relative Recruitment')

ggsave('img/Simulation_Testing/Recruitment_Scenarios.png', width=6, height=4)



# Historical Effort Plot ----
df <- Simulation$Time_Series
df$Month_n <- match(df$Month, month.abb)
df$Date <- as.Date(paste(df$Year, '-', df$Month_n, "-01", sep=""))

ggplot(df %>% filter(Sim %in% 1:4), aes(x=Date, y=F_mort, group=Sim)) +
  facet_wrap(~Sim) +
  expand_limits(y=0) +
  geom_line() +
  scale_x_date(date_breaks = "12 month", date_labels =  "%b %Y") +
  theme_bw()

ggplot(df %>% filter(Sim %in% 1:4), aes(x=Date, y=B_fished/B_unfished_eq, group=Sim)) +
  facet_wrap(~Sim) +
  expand_limits(y=0) +
  geom_line() +
  scale_x_date(date_breaks = "12 month", date_labels =  "%b %Y") +
  theme_bw()




# Simulation Testing Scenarios ----
n_months_vector <- c(12, 24, 36, 48, 60)
data_types_vector <- c('CAW', 'CAW+Effort', 'CAW+Effort+Index' )
grid <- expand.grid(n_months=n_months_vector, data_types=data_types_vector,
                    stringsAsFactors = FALSE)


## Constant Recruitment ----
for (x in 1:nrow(grid)) {
  message('Scenario: ', x, '/', nrow(grid))
  DF <- Sim_Test(x, grid=grid, Simulation, Sampling)
  # save
  fl <- paste0(paste(grid[x,], collapse = '_'), '_continuous.rda')
  saveRDS(DF, file.path('Results/Simulation_Testing/', fl))
}

## Pulse Recruitment ----
for (x in 1:nrow(grid)) {
  message('Scenario: ', x, '/', nrow(grid))
  DF <- Sim_Test(x, grid=grid, Simulation_Pulse, Sampling)
  # save
  fl <- paste0(paste(grid[x,], collapse = '_'), '_pulse.rda')
  saveRDS(DF, file.path('Results/Simulation_Testing/', fl))
}


# ---- Calculate Reference Points ----

calc_Fref_RE <- function(sim, ll, utilpow=0.3) {
  assess <- ll$assess[[sim]]
  Simulation <- ll$Simulation

  n_months_total <- ncol(Simulation$Exploitation$Sel_at_Age[sim,,])
  sel_at_age <- Simulation$Exploitation$Sel_at_Age[sim,,n_months_total]

  OM <- calculate_optimal_fishing(R0_m=Simulation$LifeHistory$R0_m,
                                  steepness=Simulation$LifeHistory$steepness,
                                  Weight_Age_Mean=Simulation$LifeHistory$Weight_Age_Mean,
                                  Maturity_at_Age=Simulation$LifeHistory$Maturity_at_Age,
                                  M_at_Age=Simulation$LifeHistory$M_at_Age,
                                  Post_Spawning_Mortality=Simulation$LifeHistory$Post_Spawning_Mortality,
                                  sel_at_age=sel_at_age,
                                  opt_type=1,
                                  utilpow=utilpow)

  Est <- calculate_optimal_fishing(R0_m=assess$rep$R0_m,
                                   steepness=assess$Data$h,
                                   Weight_Age_Mean=assess$Data$Weight_Age_Mean,
                                   Maturity_at_Age=assess$Data$Maturity_at_Age,
                                   M_at_Age=assess$Data$M_at_Age,
                                   Post_Spawning_Mortality=assess$Data$Post_Spawning_Mortality,
                                   sel_at_age=assess$rep$selA,
                                   opt_type=1, utilpow=utilpow)


  OM_DF <- Simulation$Time_Series %>% filter(Sim==sim) %>% tail(12)

  est <- median(tail(assess$rep$F_m,12)/Est$F_m)
  om <- median(OM_DF$F_mort/OM$F_m)

  data.frame(RE=(est-om)/om,
             n_months=ll$n_months,
             Data_types=ll$Data_types, Var='Fref')

}


# ---- Process Results ----

# Constant Recruitment
sim_results <- list.files("Results/Simulation_Testing")

results_list <- list()
for (i in seq_along(sim_results)) {
  message(i, '/', length(sim_results) )
  fl <- sim_results[i]
  ll <- readRDS(file.path("Results/Simulation_Testing", fl))
  nsim <- ll$assess %>% length()

  Scenario <- 'Continuous'
  if (grepl('pulse', fl)) {
    Scenario <- 'Pulse'
  }

  F_RE <- lapply(1:nsim, calc_F_RE, ll=ll) %>% do.call('rbind', .)
  F_RE$Scenario <- Scenario

  SPR_RE <- lapply(1:nsim, calc_SPR_RE, ll=ll) %>% do.call('rbind', .)
  SPR_RE$Scenario <- Scenario


  relF_RE <- lapply(1:nsim, calc_Fref_RE, ll=ll) %>% do.call('rbind', .)
  relF_RE$Scenario <- Scenario

  results_list[[i]] <- bind_rows(F_RE, SPR_RE, relF_RE)
}

results_df <- do.call('rbind', results_list)
results_df$n_months <- factor(results_df$n_months)

ggplot(results_df, aes(x=n_months, y=RE, fill=Scenario)) +
  facet_grid(Var~Data_types, scales = 'free') +
  geom_hline(yintercept = 0, linetype=2) +
  geom_boxplot() +
  theme_bw()



