assignGroups <- function(df, var="mean") {
  group <- 1
  df$groups <- NA
  for (i in 1:nrow(df)) {
    if (i == 1) {
      if (!is.na(df[[var]][i])) {
        df$groups[i] <- group + 1
      }
    } else {
      if (!is.na(df[[var]][i])) {
        df$groups[i] <- group
      } else if (is.na(df[[var]][i]) & !is.na(df[[var]][i-1])) {
        group <- group+1
      }
    }
  }
  df$x <- 1:nrow(df)
  df
}

dropLastNAs <- function(df, var='mean') {
  ind <- df[[var]] %>% is.na() %>% rev()
  ind <- min(which(!ind))
  head(df, -ind+1)
}

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}



# Generate monthly recruitment pattern
GenMonthlyRec <- function(mu=4, sigma=2) {
  n <- length(mu)
  rec <- list()
  for (i in 1:n) {
    tt <- dnorm(0:11, 5.5, sigma[i])
    tt <- tt/sum(tt)
    tt <- rep(tt,5)
    rec_pattern <- rep(0, 12)
    rec_pattern[floor(mu[i])] <- tt[6]
    j <- floor(mu[i])+1
    k <- 1

    while(j<=12) {
      rec_pattern[j] <- tt[6+k]
      j <- j+1
      k <- k+1
    }
    for (j in 1:(floor(mu[i])-1)) {
      rec_pattern[j] <- tt[6+k]
      k <- k+1
    }
    rec[[i]] <- rec_pattern
  }
  rec <- do.call('rbind', rec)
  apply(rec, 2, sum)/sum(apply(rec, 2, sum))
}

# Calculate selectivity at length
calSelL <- function(L, LF5, LFS, Vmaxlen, Linf) {

  sigma_asc <- (LFS-LF5)/sqrt(-log(0.05,2))
  sigma_dec <- (Linf-LFS)/sqrt(-log(Vmaxlen,2))

  sl <- rep(0, length(L))
  ind1 <- which(L<=LFS)
  ind2 <- which(L>LFS)
  sl[ind1] <-2^(-((L[ind1]-LFS)/sigma_asc)^2)
  sl[ind2] <- 2^(-((L[ind2]-LFS)/sigma_dec)^2)
  sl
}


#' Title
#'
#' @param scen
#'
#' @return
#' @export
#'
#' @examples
Load_Scenario <- function(scen=1) {
  Scenario_Parameters <- as.matrix(Scenario_Parameters)
  vars <- Scenario_Parameters[,1]
  vals <- Scenario_Parameters[,scen+1]
  List <- list()
  for (i in 1:length(vars)) {
    if (vars[i]!='effortPattern') {
      List[[vars[i]]] <- as.numeric(vals[i])
    } else {
      List[[vars[i]]] <- vals[i]
    }

  }
  List
}


h2CR <- function(h) {
  (4*h)/(1-h)
}

CR2h <- function(CR) {
  CR/(CR+4)
}

BH_SRR <- function(R0, h, SB, SBpR) {
  (4*R0*h*SB)/(SBpR*R0*(1-h)+(5*h-1)*SB)
}



logit <- function(p) log(p/(1-p))
inv.logit <- function(x) exp(x)/(1+exp(x))


SimData <- function(Pop, CAL_ESS=200, Effort_SD=0.2, CPUE_SD=0.2,
                    Fmeanprior=NULL,
                    fitCPUE=1, fitEffort=1,
                    use_Frwpen=0,
                    use_R0rwpen=0,
                    use_years=5) {

  # Populate data
  nyears <- dim(Pop$Number)[2]/12
  sim.timesteps <- nyears*12
  data.timesteps <- (sim.timesteps-use_years*12+1):sim.timesteps
  nts <- length(data.timesteps)
  nages <- length(Pop$Len_at_Age)

  Data <- list()
  Data$Len_Age <- Pop$Len_at_Age
  Data$SD_Len_at_Age <- Pop$SD_Len_at_Age
  Data$Wt_at_Age <- Pop$Wt_at_Age
  Data$Mat_at_Age <- Pop$Mat_at_Age
  Data$M_at_Age <- Pop$M_at_Age
  Data$phi_at_Age <- Pop$phi_at_Age
  Data$Len_Bins <- Pop$Len_Bins
  Data$Len_Mids <- Pop$Len_Mids
  Data$CAL <- Pop$CAL_samp[,data.timesteps]
  Data$CAL_ESS <- rep(CAL_ESS, nts)
  Data$Effort <- Pop$Effort_ind[data.timesteps]/mean(Pop$Effort_ind[data.timesteps],na.rm=TRUE)
  Data$Effort_SD <- rep(Effort_SD, nts)
  Data$CPUE <- Pop$Index[data.timesteps]/mean(Pop$Index[data.timesteps])
  Data$CPUE_SD <- rep(CPUE_SD, nts)

  if (!is.null(Fmeanprior)) {
    Data$use_Fmeanprior <- 1
    Data$F_meanprior <- Fmeanprior
  } else {
    Data$use_Fmeanprior <- 0
    Data$F_meanprior <- c(0.1,0.1)
  }

  Data$Fit_Effort <- fitEffort
  Data$Fit_CPUE <- fitCPUE

  Data$use_Frwpen <- use_Frwpen
  Data$use_R0rwpen <- use_R0rwpen
  Data$Fit_Effort <- fitEffort
  Data$Fit_CPUE <- fitCPUE

  out <- list()
  out$Pop <- Pop
  out$Data <- Data
  out
}

