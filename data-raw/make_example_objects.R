
LifeHistory <- Import_LifeHistory()
Exploitation <- Import_Exploitation()
Sampling <- Import_Sampling()

usethis::use_data(LifeHistory, overwrite = TRUE)
usethis::use_data(Exploitation, overwrite = TRUE)
usethis::use_data(Sampling, overwrite = TRUE)
