## code to prepare `DATASET` dataset goes here
Scenario_Parameters <- readxl::read_excel('data-raw/Scenario_Parameters.xlsx')
usethis::use_data(Scenario_Parameters, overwrite = TRUE)

casestudydata <- readRDS('data-raw/casestudydata.rds')
usethis::use_data(casestudydata)

