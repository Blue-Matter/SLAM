
read_pars <- function(Indata, list, name, ncol) {
  if (ncol==1) ind <- 2
  if (ncol>1) ind <- 2:(2+ncol-1)
  list[[name]] <- as.numeric(Indata[Indata$Parameter==name,ind])
  list
}


#' Import the Life History Parameters from a CSV file
#'
#' @param csv The name of the CSV file containing the parameters
#' @param dir The directory the CSV is in (defaults to current working directory)
#' @export
#'
Import_LifeHistory <- function(csv='Parameters.csv', dir=getwd()) {
  LifeHistory <- list()

  Indata <- read.csv(file.path(dir, csv))

  LifeHistory <- read_pars(Indata, LifeHistory, 'maxage', 1)
  LifeHistory$Ages <- 0:LifeHistory$maxage
  LifeHistory$nAge <- length(LifeHistory$Ages)
  LifeHistory <- read_pars(Indata, LifeHistory, 'Weight_Age_Mean', LifeHistory$nAge)
  LifeHistory <- read_pars(Indata, LifeHistory, 'Weight_Age_SD', LifeHistory$nAge)
  LifeHistory <- read_pars(Indata, LifeHistory, 'M_at_Age', LifeHistory$nAge)
  LifeHistory <- read_pars(Indata, LifeHistory, 'Maturity_at_Age',   LifeHistory$nAge)
  LifeHistory <- read_pars(Indata, LifeHistory, 'Post_Spawning_Mortality', LifeHistory$nAge)

  LifeHistory <- read_pars(Indata, LifeHistory, 'R0_bar', 1)
  LifeHistory <- read_pars(Indata, LifeHistory, 'R0_m', 12)
  LifeHistory <- read_pars(Indata, LifeHistory, 'sigmaR', 1)
  LifeHistory <- read_pars(Indata, LifeHistory, 'steepness', 1)

  LifeHistory
}


#' Import the Exploitation Parameters from a CSV file
#' @rdname Import_LifeHistory
#' @export
Import_Exploitation <- function(csv='Parameters.csv', dir = getwd()) {
  Exploitation  <- list()
  Indata <- read.csv(file.path(dir, csv))
  Exploitation <- read_pars(Indata, Exploitation, 'SA50', 1)
  Exploitation <- read_pars(Indata, Exploitation, 'SA95', 1)
  Exploitation <- read_pars(Indata, Exploitation, 'nyears', 1)
  nyears <- Exploitation$nyears
  Exploitation <- read_pars(Indata, Exploitation, 'Effort_Annual_Mean', nyears)
  Exploitation <- read_pars(Indata, Exploitation, 'Effort_Annual_SD', nyears)
  Exploitation <- read_pars(Indata, Exploitation, 'Effort_Month_Mean', 12)
  Exploitation <- read_pars(Indata, Exploitation, 'Effort_Month_SD', 12)
  Exploitation <- read_pars(Indata, Exploitation, 'q', 1)

  Exploitation
}

#' Import the Data Parameters from a CSV file
#' @rdname Import_LifeHistory
#' @export
Import_Data <- function(csv='Parameters.csv', dir = getwd()) {
  Data <- list()

  Indata <- read.csv(file.path(dir, csv))

  Data <- read_pars(Indata, Data, 'n_recent_months', 1)
  Data <- read_pars(Indata, Data, 'Rel_Sample_Month', 12)
  Data <- read_pars(Indata, Data, 'CPUE_CV', 1)
  Data <- read_pars(Indata, Data, 'Catch_CV', 1)
  Data <- read_pars(Indata, Data, 'Effort_CV', 1)
  Data <- read_pars(Indata, Data, 'Weight_Bin_Max', 1)
  Data <- read_pars(Indata, Data, 'Weight_Bin_Width', 1)
  Data <- read_pars(Indata, Data, 'CAW_Annual_Sample_Size', 1)
  Data <- read_pars(Indata, Data, 'CAW_Annual_ESS', 1)

  Data

}
