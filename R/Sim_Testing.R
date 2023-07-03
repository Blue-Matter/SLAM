#' Run a Simulation Test
#'
#' @param x Scenario number. Row of `grid`
#' @param grid data.frame with simulation grid. Columns `n_months` and `data_types`
#' @param Simulation An object of class `Simulation`
#' @param Sampling  An object of class `Sampling`
#'
#' @return A list
#' @export
#'
Sim_Test <- function(x, grid, Simulation, Sampling) {

  # Generate Data
  Sampling$n_recent_months <- grid$n_months[x]
  Sampled_Data <- Generate_Data(Simulation, Sampling = Sampling)

  nsim <- max(Simulation$At_Age_Time_Series$Sim)
  # Do assessment
  outlist <- list()
  for (i in 1:nsim) {
    message(i, '/', nsim)
    Data <- Import(Sampled_Data, sim=i, Data_types = grid$data_types[x])

    Parameters <- Initialize_Parameters(Data)

    assess <- Assess(Data, Parameters)

    outlist[[i]] <- assess
  }

  out <- list()
  out$assess <- outlist
  out$n_months <-  grid$n_months[x]
  out$Data_types <- grid$data_types[x]
  out$Simulation <- Simulation
  out$Sampling <- Sampling
  out
}




#' Title
#'
#' @param sim
#' @param ll
#'
#' @return
#' @export
#'
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

  df <- df %>% tail(n_months)
  # data.frame(RE=(median(df$Estimate)-median(df$OM))/median(df$OM),
  #            n_months=ll$n_months, Data_types=ll$Data_types, Var='F')
  data.frame(RE=median((df$Estimate-df$OM)/df$OM),
             n_months=ll$n_months, Data_types=ll$Data_types, Var='F')

}

#' Title
#'
#' @param sim
#' @param ll
#'
#' @return
#' @export
#'
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

  df <- df %>% tail(n_months)
  data.frame(RE=median((df$Estimate-df$OM)/df$OM),
             n_months=ll$n_months, Data_types=ll$Data_types, Var='SB/SB0')
}

#' Title
#'
#' @param sim
#' @param ll
#'
#' @return
#' @export
#'
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

  df <- df %>% tail(n_months)
  data.frame(RE=median((df$Estimate-df$OM)/df$OM),
             n_months=ll$n_months, Data_types=ll$Data_types, Var='SPR')

}

