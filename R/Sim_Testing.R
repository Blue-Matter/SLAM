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

