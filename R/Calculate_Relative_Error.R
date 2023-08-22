calc_F_RE <- function(sim, ll) {

  assess <- ll$assess[[sim]]
  if(class(assess$opt) =='try-error')
    return(NULL)
  n_months <- ll$n_months
  OM <- ll$Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% assess$Data$Month_ind)
  Year <- sort(assess$Data$Year)
  Month <- assess$Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   OM=OM$F_mort,
                   Estimate=assess$rep$F_m)

  df <- df %>% tail(12)
  # data.frame(RE=(median(df$Estimate)-median(df$OM))/median(df$OM),
  #            n_months=ll$n_months, Data_types=ll$Data_types, Var='F')
  data.frame(RE=median((df$Estimate-df$OM)/df$OM),
             n_months=ll$n_months, Data_types=ll$Data_types, Var='F')

}


calc_SB_SB0_RE <- function(sim, ll) {
  assess <- ll$assess[[sim]]
  if(class(assess$opt) =='try-error')
    return(NULL)

  n_months <- ll$n_months
  OM <- ll$Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% assess$Data$Month_ind)
  Year <- sort(assess$Data$Year)
  Month <- assess$Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   OM=OM$SB_fished/OM$SB_unfished_eq,
                   Estimate=assess$rep$SB_m/assess$rep$SB0_m)

  df <- df %>% tail(12)
  data.frame(RE=median((df$Estimate-df$OM)/df$OM),
             n_months=ll$n_months, Data_types=ll$Data_types, Var='SB/SB0')
}


calc_SPR_RE <- function(sim, ll) {
  assess <- ll$assess[[sim]]
  if(class(assess$opt) =='try-error')
    return(NULL)
  n_months <- ll$n_months
  OM <- ll$Simulation$Time_Series %>%
    filter(Sim==sim, Month_ind %in% assess$Data$Month_ind)
  Year <- sort(assess$Data$Year)
  Month <- assess$Data$Month
  Year_Month <- paste(Year, Month, sep="_")

  df <- data.frame(Year=Year,
                   Month=Month,
                   Year_Month=Year_Month,
                   OM=OM$SPR,
                   Estimate=assess$rep$SPR)

  df <- df %>% tail(12)
  data.frame(RE=median((df$Estimate-df$OM)/df$OM),
             n_months=ll$n_months, Data_types=ll$Data_types, Var='SPR')

}

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
