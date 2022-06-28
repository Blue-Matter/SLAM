## code to prepare `DATASET` dataset goes here
# Scenario_Parameters <- readxl::read_excel('data-raw/Scenario_Parameters.xlsx')

casestudydata <- readRDS('data-raw/casestudydata.rds')
usethis::use_data(casestudydata, overwrite = TRUE)
# usethis::use_data(Scenario_Parameters, casestudydata, overwrite = TRUE)

