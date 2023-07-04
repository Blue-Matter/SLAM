# remotes::install_github('blue-matter/SLAM')

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
df1 <- Simulation$Time_Series
df2 <- Simulation_Pulse$Time_Series
df1$Month_n <- match(df1$Month, month.abb)
df1$Date <- as.Date(paste(df1$Year, '-', df1$Month_n, "-01", sep=""))

df2$Month_n <- match(df2$Month, month.abb)
df2$Date <- as.Date(paste(df2$Year, '-', df2$Month_n, "-01", sep=""))

df1$Scenario <- 'Constant'
df2$Scenario <- 'Pulse'
effort_df <- bind_rows(df1, df2)

ggplot(effort_df %>% filter(Sim %in% 1), aes(x=Date, y=F_mort, group=Sim)) +
  facet_grid(~Scenario) +
  expand_limits(y=0) +
  geom_line() +
  scale_x_date(date_breaks = "12 month", date_labels =  "%Y",
               limits =range(effort_df$Date),
               expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=1)) +
  labs(x='Month', y='Fishing mortality (F)') +
  scale_y_continuous(expand = c(0, 0))

ggsave('img/Simulation_Testing/Fishing_Mortality.png', width=8, height=4)


# ---- Life History Plots -----
Ages <- Simulation$LifeHistory$Ages
Weight_Age_Mean <- Simulation$LifeHistory$Weight_Age_Mean
Weight_Age_SD <- Simulation$LifeHistory$Weight_Age_SD

age_df <- data.frame(Age=Ages, Mean=Weight_Age_Mean)
mu <- log(Weight_Age_Mean) -0.5*Weight_Age_SD^2
age_df$Upper <- exp(qnorm(0.75, mu, Weight_Age_SD ))
age_df$Lower <- exp(qnorm(0.25, mu, Weight_Age_SD ))

ggplot(age_df, aes(x=Age, ymin=Lower, ymax=Upper, y=Mean)) +
  geom_ribbon(fill='grey') +
  geom_line() +
  theme_bw() +
  labs(x='Age (months)', y='Weight (kg)')

ggsave('img/Simulation_Testing/Weight_Age.png', width=4, height=4)


Ages <- Simulation$LifeHistory$Ages
Maturity_at_Age <- Simulation$LifeHistory$Maturity_at_Age
Post_Spawning_Mortality <- Simulation$LifeHistory$Post_Spawning_Mortality

age_df <- data.frame(Age=Ages, Maturity=Maturity_at_Age,
                     PSM=Post_Spawning_Mortality)

ggplot(age_df, aes(x=Age, y=Maturity)) +
  geom_line() +
  theme_bw() +
  labs(x='Age (months)', y='Probability of Spawning')

ggsave('img/Simulation_Testing/Prob_Spawning.png', width=4, height=4)





# Simulation Testing Scenarios ----
n_months_vector <- c(12, 24, 36, 48, 60)
data_types_vector <- c('CAW',
                       'CAW+Effort',
                       'CAW+Catch',
                       'CAW+Effort+Index',
                       'CAW+Catch+Index')
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

calc_Fref_RE <- function(sim, ll, utilpow=0.4) {
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

  est <- (tail(assess$rep$F_m,12)/Est$F_m)
  om <- (OM_DF$F_mort/OM$F_m)

  data.frame(RE=median((est-om)/om),
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

  # SB_SB0_RE <- lapply(1:nsim, calc_SB_SB0_RE, ll=ll) %>% do.call('rbind', .)
  # SB_SB0_RE$Scenario <- Scenario

  relF_RE <- lapply(1:nsim, calc_Fref_RE, ll=ll) %>% do.call('rbind', .)
  relF_RE$Scenario <- Scenario

  results_list[[i]] <- bind_rows(F_RE, SPR_RE, relF_RE)
}

results_df <- do.call('rbind', results_list)
results_df$n_months <- factor(results_df$n_months)

results_df <- results_df %>% filter(Var!='Fref')
ggplot(results_df, aes(x=n_months, y=RE, fill=Scenario)) +
  facet_grid(Var~Data_types, scales = 'free') +
  geom_hline(yintercept = 0, linetype=2) +
  geom_boxplot() +
  theme_bw() +
  labs(x='Number of Months', y='Relative Error')

ggsave('img/Simulation_Testing/Relative_Error.png')

