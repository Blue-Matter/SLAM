library(SLAM)
library(ggplot2)
library(flextable)

Data.Dir <- 'G:/My Drive/1_PROJECTS/Biospherics_Indo_Octopus/Data'

Case_Study_Data <- read.csv(file.path(Data.Dir, 'Case_Study_Data.csv'))

Case_Study_Sites <- read.csv(file.path(Data.Dir, 'Case_Study_Sites.csv'))


# 1. Summarize Data - number of months of data ----
summary_list <- list()
for (i in 1:length(Data_Villages)) {
  village <- Data_Villages[i]
  summary_list[[i]] <- Summarize_Data(village, Case_Study_Data, Case_Study_Sites)
}
summary_df <- do.call('rbind', summary_list)


table_df <- summary_df %>% group_by(Sites) %>% summarize(n_site=unique(n_sites),
                                                        First_date=min(Date),
                                                        Last_date=max(Date),
                                                        n_months=unique(n_months),
                                                        missing_months=unique(missing_months),
                                                        CAW_nsamp=round(mean(n_CAW, na.rm=TRUE),0),
                                                        Effort=round(mean(Effort, na.rm=TRUE),0),
                                                        Catch=round(mean(Catch, na.rm=TRUE),0))

table_df$First_date <- format(as.Date(table_df$First_date), "%Y-%m")
table_df$Last_date <- format(as.Date(table_df$Last_date), "%Y-%m")


ft <- flextable(table_df)
ft %>% save_as_docx( path = "Results/Case_Studies/Tables/Summary_Table.docx")

# 2. Filter Villages with at least 12 months of CAW Data
data_df <- table_df %>% filter((n_months-missing_months)>=12)


# 3. Process Data for each Site
data_list <- list()
for (i in 1:length(data_df$Sites)) {
  Sites <- data_df$Sites[i]
  data_list[[i]] <- Process_Data(Sites, Case_Study_Data)
}


# 4. Make Data Objects
save.dir <- file.path(Data.Dir, 'Casestudies')
for (i in 1:length(data_df$Sites)) {
  tt <- Make_Data_Objects(i, data_list, Case_Study_Sites, save.dir=save.dir)
}



# 5. Fit Assessment Model
data_files <- list.files(save.dir)

i <- 9
data <- readRDS(file.path(save.dir, data_files[i]))
data$Year
data$Month

data$CAW

data$CAW_ESS[data$CAW_ESS>50] <- 50
data$Effort_Mean
data$Effort_SD
data$Index_Mean
data$Index_SD

MyAssess <- Assess(data, Est_Seasonal = TRUE)


MyAssess$rep$nll_joint

plot(MyAssess$rep$F_m, type='l')
plot(MyAssess$rep$SPR, type='l')
plot(MyAssess$rep$R0_m, type='l')

Report(MyAssess)

# 6. Figures




## CAW
dflist <- list()
for (i in 1:length(data_list2)) {
  dflist[[i]] <- data_list2[[i]]$CAW
}
df_CAW <- do.call('rbind', dflist)
village <- unique(df_CAW$Village)[3]
df_CAW2 <- df_CAW %>% filter(Village==village)
df_CAW2 <- df_CAW2 %>% group_by(Year, Month_Name) %>%
  mutate(Freq=Count/sum(Count))


ggplot(df_CAW2, aes(x=Bin, y=Freq)) +
  facet_grid(Year~Month_Name) +
  geom_bar(stat='identity') +
  theme_bw() +
  labs(x='Weight class (kg)', y='Frequency')



## Effort
dflist <- list()
for (i in 1:length(data_list2)) {
  dflist[[i]] <- data_list2[[i]]$Effort
}
df_Effort <- do.call('rbind', dflist)

ggplot(df_Effort, aes(x=Date2, y=Effort)) +
  facet_grid(~Village) +
  geom_point() +
  geom_line() +
  expand_limits(y=0) +
  theme_bw()


# 5. Make Data objects



data <- Make_Data_Objects(2, data_list2)

data$Fit_Effort <- 1
data$CAW_ESS[data$CAW_ESS>200] <- 200

MyAssess <- Assess(data, Est_Seasonal = TRUE)


rep <- MyAssess$obj$report()
rep$nll_joint
MyAssess$obj$gr()

rep <- MyAssess$rep
rep$F_m
rep$F_minit



rep$Mean_monthly_F

rep$nll_joint

data$CAW[1,]
rep$predCAW[,1]

length(rep$CAWnll)


rep$F_minit
rep$F_m

rep$Fa_init
rep$Za_init
rep$N_m



MyAssess$obj$hessian

MyAssess$opt

Report(MyAssess)



Data <- Import(Example_Data()[1])

Data




















