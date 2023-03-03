
#' Example Life History Parameters
#'
#'
#' @format ## `Life History`
#' A named list:
#' \describe{
#'   \item{Stock_Name}{A name for the }
#'   \item{Species}{}
#'   \item{maxage}{}
#' }
"LifeHistory"



#' Get the file path of the system files for the SLAM package
#'
#' @param file The name of the file to locate
#' @return The location of the file included in the SLAM package
#' @export
#'
#' @examples
#' File_Location()
File_Location <- function(file='Simulation_Parameters.csv') {
  if (is.null(file))
    return(system.file(package="SLAM"))
  system.file(file, package="SLAM")
}

read_pars <- function(Indata, list, name, ncol, is.numeric=TRUE) {
  if (ncol==1) ind <- 2
  if (ncol>1) ind <- 2:(2+ncol-1)
  value <- trimws(Indata[Indata$Parameter==name,ind])
  if (is.numeric) value <- as.numeric(value)
  list[[name]] <- value
  list
}


#' Import the Life History, Exploitation, and Sampling Parameters
#'
#' The Life History, Exploitation, and Sampling Parameters can imported from a correctly
#' formatted CSV file. See `Simulation_Parameters.csv` included with the `SLAM` package
#' (found at `File_Location()`) for an example.
#'
#' @param csv The path and name of the CSV file containing the parameters
#' @return A named list with the life-history, exploitation, and sampling parameters
#' @export
Import_Parameters <- function(csv=NULL) {
  if (is.null(csv))
    csv <- system.file("Simulation_Parameters.csv", package="SLAM")
  out <- list()
  out$LifeHistory <- Import_LifeHistory(csv)
  out$Exploitation <- Import_Exploitation(csv)
  out$Sampling <- Import_Sampling(csv)
  out
}

#' Import the Life History Parameters from a CSV file
#' @rdname Import_Parameters
#' @return A named list with the life history parameters
#' @export
Import_LifeHistory <- function(csv=NULL) {
  if (is.null(csv))
    csv <- system.file("Simulation_Parameters.csv", package="SLAM")
  Indata <- read.csv(csv)

  LifeHistory <- list()
  LifeHistory <- read_pars(Indata, LifeHistory, 'Stock_Name', 1, FALSE)
  LifeHistory <- read_pars(Indata, LifeHistory, 'Species', 1, FALSE)
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
#' @rdname Import_Parameters
#' @return A named list with the exploitation parameters
#' @export
Import_Exploitation <- function(csv=NULL) {
  if (is.null(csv))
    csv <- system.file("Simulation_Parameters.csv", package="SLAM")
  Indata <- read.csv(csv)

  Exploitation  <- list()
  Exploitation <- read_pars(Indata, Exploitation, 'Fleet_Name', 1, FALSE)
  Exploitation <- read_pars(Indata, Exploitation, 'SA50', 1)
  Exploitation <- read_pars(Indata, Exploitation, 'SA95', 1)
  Exploitation <- read_pars(Indata, Exploitation, 'nmonths', 1)

  if (ncol(Indata)>=Exploitation$nmonths) {
    Exploitation <- try(read_pars(Indata, Exploitation, 'Effort_Pattern', Exploitation$nmonths), silent=TRUE)
    if (inherits(Exploitation, 'try-error')) {
      stop('Could not import information from row `Effort_Pattern`. Is it numerical values of length `nmonths`?')
    }
    Exploitation$Effort_month <- Exploitation$Effort_Pattern
    Exploitation$Effort_current <- Exploitation$Effort_Pattern[length(Exploitation$Effort_Pattern)]
  } else {
    Exploitation <- read_pars(Indata, Exploitation, 'Effort_pattern', 1, FALSE)
    Exploitation <- read_pars(Indata, Exploitation, 'Effort_current', 1)
    Exploitation <- read_pars(Indata, Exploitation, 'Effort_cv', 1)
  }

  Exploitation <- read_pars(Indata, Exploitation, 'q', 1)
  Exploitation <- read_pars(Indata, Exploitation, 'q_cv', 1)
  Exploitation <- read_pars(Indata, Exploitation, 'HARA_power', 1)
  Exploitation
}


#' Import the Sampling Parameters from a CSV file
#' @rdname Import_Parameters
#' @export
Import_Sampling <- function(csv=NULL) {
  if (is.null(csv))
    csv <- system.file("Simulation_Parameters.csv", package="SLAM")
  Indata <- read.csv(csv)

  Data <- list()
  Data <- read_pars(Indata, Data, 'Sampling_Name', 1, FALSE)
  Data <- read_pars(Indata, Data, 'n_recent_months', 1)
  Data <- read_pars(Indata, Data, 'Rel_Sample_Month', 12)
  Data <- read_pars(Indata, Data, 'CPUE_CV', 1)
  Data <- read_pars(Indata, Data, 'Catch_CV', 1)
  Data <- read_pars(Indata, Data, 'Effort_CV', 1)
  Data <- read_pars(Indata, Data, 'Weight_Bin_Max', 1)
  Data <- read_pars(Indata, Data, 'Weight_Bin_Width', 1)
  Data <- read_pars(Indata, Data, 'CAW_Annual_Sample_Size', 1)
  Data <- read_pars(Indata, Data, 'CAW_Annual_ESS', 1)
  Data <- read_pars(Indata, Data, 'CAA_Annual_Sample_Size', 1)
  Data <- read_pars(Indata, Data, 'CAA_Annual_ESS', 1)
  Data

}
