
LifeHistory <- Import_LifeHistory(dev=TRUE)
Exploitation <- Import_Exploitation(dev=TRUE)
Sampling <- Import_Sampling(dev=TRUE)

LifeHistory_Pulse <- LifeHistory
LifeHistory_Pulse$R0_m <- Generate_Monthly_Recruitment(6.5, 2.5)

LifeHistory_No_Error <- LifeHistory
LifeHistory_No_Error$sigmaR <- 0.0001
LifeHistory_No_Error$R0_bar <- 1

LifeHistory_Pulse_No_Error <- LifeHistory_Pulse
LifeHistory_Pulse_No_Error$sigmaR <- 0.0001
LifeHistory_Pulse_No_Error$R0_bar <- 1


Exploitation_No_Error <- Exploitation
Exploitation_No_Error$q_cv <- 0.0001
Exploitation_No_Error$Effort_cv <- 0.0001


Perfect_Sampling <- SLAM::Sampling
Perfect_Sampling$CPUE_CV <- 0.0001
Perfect_Sampling$Catch_CV <- 0.0001
Perfect_Sampling$Effort_CV <- 0.0001
Perfect_Sampling$CAW_Annual_ESS <- 1E6
Perfect_Sampling$CAW_Annual_Sample_Size <- 1E6
Perfect_Sampling$CAA_Annual_ESS <- 1E6
Perfect_Sampling$CAA_Annual_Sample_Size <- 1E6


usethis::use_data(LifeHistory, overwrite = TRUE)
usethis::use_data(LifeHistory_Pulse, overwrite = TRUE)
usethis::use_data(Exploitation, overwrite = TRUE)
usethis::use_data(Sampling, overwrite = TRUE)

usethis::use_data(LifeHistory_No_Error, overwrite = TRUE)
usethis::use_data(Exploitation_No_Error, overwrite = TRUE)
usethis::use_data(Perfect_Sampling, overwrite = TRUE)
