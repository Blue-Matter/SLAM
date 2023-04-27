library(dplyr)

dir <- 'G:/My Drive/1_PROJECTS/Biospherics_Indo_Octopus/Data'

data <- readxl::read_excel(file.path(dir, 'FOR ADRIAN_Update.xlsx'))

data$Village %>% unique() %>% sort()
# Darawa not there?

data$Fisher %>% unique()
data$`Only fishing?` %>% unique()

# Only fishing? = No
# How can we calculate CPUE?
data %>% group_by(`Only fishing?`) %>%
  summarise(n=n())

# Effort - large values 30 +?
data_effort <- data %>% filter(`Fishing Duration (hour)`!='Unrecorded')
data_effort$`Fishing Duration (hour)` %>% as.numeric() %>% unique() %>%  sort()

data_effort$`Fishing Duration (hour)` %>% as.numeric() %>% hist()




data$`Fishing Duration (hour)` %>% table()

xtabs((data_effort$`Fishing Duration (hour)`) ~ data_effort$`Only fishing?`)
tt = data %>% filter(`Fishing Duration (hour)`!='Unrecorded', `Only fishing?`=='No')
tt$`Fishing Duration (hour)` %>% as.numeric() %>% sum()
