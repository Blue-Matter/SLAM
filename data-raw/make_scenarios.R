

# ---- Make Scenario Parameter Files ----
library(usethis)
library(SLAM)
library(dplyr)

# Scenarios:

Monthly_Recruitment_Pattern <- c('Constant', 'Pulse')
Data_n_months <- c(1, 3, 12, 60, 120, 360)
Data_types <- c('CAW', 'CAW+Index', 'CAW+Effort', 'Index+Effort', 'CAW+Index+Effort')
Conditions <- c('Idealized', 'Process+Observation Error')

Scenario_Grid <- expand.grid(Monthly_Recruitment_Pattern=Monthly_Recruitment_Pattern,
                             Data_n_months=Data_n_months,
                             Data_types=Data_types,
                             Conditions=Conditions,
                             stringsAsFactors = FALSE)

names <- apply(Scenario_Grid, 1, paste0, collapse="_")
names <- gsub(" ", "", names, fixed = TRUE)

Scenario_Grid$Name <- names

Scenario_Grid_drop <- Scenario_Grid %>% filter(Data_n_months==1 , Data_types !="CAW")

Scenario_Grid <- Scenario_Grid[!(Scenario_Grid$Name %in% Scenario_Grid_drop$Name),]

make_scenario_data <- function(i, Scenario_Grid) {

  Scenario <- Scenario_Grid[i,]
  LifeHistory <- Import_LifeHistory(dir='inst')
  Exploitation <- Import_Exploitation(dir='inst')

  Data <- Import_Data(dir='inst')

  Rec_Pattern_Pars <- switch(Scenario$Monthly_Recruitment_Pattern,
                             'Constant'=list(1, 1000),
                             'Pulse'=list(6.5, 1),
                             'Diffuse'=list(6.5,3),
                             'Bi-Modal'=list(c(3.5,8.5), c(1,1)))

  LifeHistory$R0_m <- Generate_Monthly_Recruitment(Rec_Pattern_Pars[[1]], Rec_Pattern_Pars[[2]])

  Data$n_recent_months <- Scenario$Data_n_months

  Data_types <- strsplit(Scenario$Data_types, '\\+')[[1]]

  Condition_Parameters <- switch(Scenario$Conditions,
                                 'Idealized'=list(sigmaR=0.01,
                                                  Effort_Month_SD=0.01,
                                                  CPUE_CV=0.01,
                                                  Catch_CV=0.01,
                                                  Effort_CV=0.01,
                                                  CAW_Annual_Sample_Size=10000,
                                                  CAW_Annual_ESS=10000),
                                 'Process+Observation Error'=list(sigmaR=0.6,
                                                                  Effort_Month_SD=0.1,
                                                                  CPUE_CV=0.3,
                                                                  Catch_CV=0.3,
                                                                  Effort_CV=0.3,
                                                                  CAW_Annual_Sample_Size=1000,
                                                                  CAW_Annual_ESS=500))

  nm1 <- intersect(names(Condition_Parameters), names(LifeHistory))
  LifeHistory <- modifyList(LifeHistory, Condition_Parameters[nm1])

  nm1 <- intersect(names(Exploitation), names(Condition_Parameters))
  Exploitation <- modifyList(Exploitation, Condition_Parameters[nm1])

  nm1 <- intersect(names(Data), names(Condition_Parameters))
  Data <- modifyList(Data, Condition_Parameters[nm1])

  month_opt_F <- calculate_optimal_fishing(LifeHistory, Exploitation, opt_type=1, utilpow=0.3)$F_m

  Exploitation$Effort_Month_Mean <- month_opt_F/sum(month_opt_F)

  out <- list(LifeHistory=LifeHistory,
              Exploitation=Exploitation,
              Data=Data)

  name <- Scenario$Name
  assign(name, out)

  do.call("use_data", list(as.name(name), overwrite = TRUE))
}

sapply(1:nrow(Scenario_Grid), function(i)
  make_scenario_data(i=i, Scenario_Grid = Scenario_Grid))




usethis::use_data(Scenario_Grid, overwrite = TRUE)






