library(SLAM)
library(ggplot2)
library(flextable)
require(maps)
require(mapdata)
library(ggthemes)
library(ggrepel)

Data.Dir <- 'G:/My Drive/1_PROJECTS/Biospherics_Indo_Octopus/Data'

Case_Study_Data <- read.csv(file.path(Data.Dir, 'Case_Study_Data.csv'))

Case_Study_Sites <- read.csv(file.path(Data.Dir, 'Case_Study_Sites.csv'))

Data_Villages <- Case_Study_Data$Village %>% unique() %>% sort()
length(Data_Villages)


# 1. Aggregate Village Data into Sites ----
village_list <- list()
for (i in 1:length(Data_Villages)) {
  village <- Data_Villages[i]
  village_list[[i]] <- Aggregate_Data(village, Case_Study_Data, Case_Study_Sites)
}
village_df <- do.call('rbind', village_list)

Sites <- village_df$Sites %>% unique()
Sites %>% length()
village_df$Village  %>% unique() %>% length()

# 2. Summarize Data - number of months of data ----
summary_list <- list()
for (i in seq_along(Sites)) {
  sites <- Sites[i]
  summary_list[[i]] <- Summarize_Data(sites, Case_Study_Data, Case_Study_Sites)
}

summary_df <- do.call('rbind', summary_list)

summary_df$Site %>% unique()
summary_df$Site %>% unique() %>% length()

villages <- summary_df$Site %>% strsplit(., ',') %>% unlist() %>% trimws() %>% unique()
length(villages)
Data_Villages[!Data_Villages %in% villages]

# 3. Filter for less than 2 months of sequential missing CAW data ----
subset_df <- function(df, ...) {

  if (nrow(df)<12) return(NULL)

  # dummy value for NA
  for (i in 2:(nrow(df)-1)) {
    if (is.na(df$n_CAW[i]) & !is.na(df$n_CAW[i+1]) & !is.na(df$n_CAW[i-1]))
              df$n_CAW[i] <- Inf

  }

  # filter to longest continuous data set
  group <- 1
  df$group <- group
  df$run <- 1
  df$keep <- TRUE
  for (i in 2:nrow(df)) {
    if (!is.na(df$n_CAW[i]) & !is.na(df$n_CAW[i-1])) {
      df$run[i] <- df$run[i-1] + 1
      df$group[i] <- group
      df$keep[i] <- TRUE
    } else {
      df$run[i] <- 0
      group <- group+1
      df$group[i] <- group
      df$keep[i] <- FALSE
    }
  }
  tt <- df %>% group_by(group, keep) %>% summarize(rr=max(run))
  top <- tt %>% ungroup() %>% filter(keep==TRUE, rr==max(rr))
  df <- df %>% filter(group==top$group)

  # replace dummy value with NA again
  df$n_CAW[!is.finite(df$n_CAW)] <- NA
  df$n_months <- nrow(df)
  df$run <- df$group <- df$keep <- NULL
  if (nrow(df)<12) return(NULL)
  df
}


summary_df2 <- summary_df %>% group_by(Site) %>% group_map(., subset_df, .keep = TRUE) %>%
  do.call('rbind', .)

summary_df2$Site %>% unique()
summary_df2$Site %>% unique() %>% length()


# 4. Filter for at 12 months of CAW data ----
# already done above
summary_df2 %>% filter(n_months<12) %>% reframe(unique(Site))

summary_df3 <- summary_df2 %>% filter(n_months>=12)

summary_df3$Site %>% unique() %>% length()




# 5. Summarize n months of data

summary_df3 <- summary_df3 %>% group_by(Site) %>%
  mutate(n_CAW_months=sum(!is.na(n_CAW)),
         n_Effort_months=sum(is.finite(Effort)),
         n_Catch_months=sum(is.finite(Catch)))

table_df <- summary_df3 %>% group_by(Site) %>% summarize(n_site=unique(n_sites),
                                                        First_date=min(Date),
                                                        Last_date=max(Date),
                                                        n_months=unique(n_months),
                                                        n_CAW_months=unique(n_CAW_months),
                                                        n_Effort_months=unique(n_Effort_months),
                                                        n_Catch_months=unique(n_Catch_months),
                                                        CAW_nsamp=round(mean(n_CAW, na.rm=TRUE),0),
                                                        Effort=round(mean(Effort, na.rm=TRUE),0),
                                                        Catch=round(mean(Catch, na.rm=TRUE),0))

table_df$First_date <- format(as.Date(table_df$First_date), "%Y-%m")
table_df$Last_date <- format(as.Date(table_df$Last_date), "%Y-%m")
table_df$Date <- paste(table_df$First_date, table_df$Last_date, sep=" -- ")
table_df <- table_df %>% mutate(Number=1:nrow(table_df)) %>%
  select(Number, Site, n_village=n_site, Date, n_months, n_CAW_months,
         n_Effort_months,  n_Catch_months, CAW_nsamp, Effort, Catch)

table_df <- table_df %>% arrange(Site)

ft <- flextable(table_df)
ft %>% save_as_docx(path = "Results/Case_Studies/Tables/Summary_Table.docx")


# Number of villages
data_villages <- table_df$Site %>% strsplit(., ',') %>% unlist() %>% unique() %>% trimws()
data_villages %>% length()

# Map of Sites

Sites_df <- Case_Study_Sites %>% select(Village, Combined, Latitude, Longitude,
                            Province=Province.Provinsi) %>%
  filter(Village %in% data_villages)

missing_villages <- data_villages[!data_villages %in% Sites_df$Village]
Sites_df$Combined <- factor(Sites_df$Combined)


ylim <- Sites_df$Latitude %>% range()
xlim <- Sites_df$Longitude %>% range()
world <- map_data('world')

Sites_df$Village == 'Tetandara (Arubura)'

map <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group),
               fill = "lightgreen", color = "black") +
  coord_fixed(1.3) +
  geom_point(data = Sites_df,aes(x = Longitude, y = Latitude, color=Combined),
             size=2, alpha = 0.5, show.legend = F) +
  geom_label_repel(data=Sites_df, aes(x = Longitude, y = Latitude, label=Village,
                                      color=Combined),
                  max.overlaps = 25) +
  xlim(95.011, 141.02) + ylim(-11.10718, 5.907130) +
  # xlim(xlim) + ylim(ylim) +
  guides(color='none') +
  labs(x='Longitude', y='Latitude')

map
# ggsave('img/Case_Studies/map.png', map, width=8, height=6.15)



# 3. Process Data for each Site ----
data_list <- list()
for (i in 1:length(table_df$Site)) {
  Sites <- table_df$Site[i]
  Use_Dates <- table_df$Date[i]

  Use_Dates <- lubridate::ymd(paste0(strsplit(Use_Dates, '--')[[1]], '-01'))

  data_list[[i]] <- Process_Data(Sites, Case_Study_Data, Use_Dates)
}



# 4. Make Data Objects  ----
save.dir <- file.path(Data.Dir, 'Casestudies')

for (i in 1:length(table_df$Site)) {
  tt <- Make_Data_Objects(i, data_list, Case_Study_Sites, save.dir=save.dir)
}

# 5. Fit Assessment Model  ----
data_files <- list.files(save.dir, pattern='.rdata')

# Write Data CSVs
for (i in seq_along(data_files)) {
  Data <- readRDS(file.path(save.dir,data_files[[i]]))

  csvfile <- paste0(tools::file_path_sans_ext(data_files[[i]]), '.csv')

  Write_Data2CSV(Data, csvfile, dir=save.dir, overwrite = TRUE)

}


for (i in seq_along(data_files)) {
  data <- readRDS(file.path(save.dir, data_files[i]))

  Parameters <- Initialize_Parameters(data, sigmaF_m=0.3)

  DoAssess <- Assess(data, Parameters = Parameters)

  maxit <- 0
  while(!all(DoAssess$chk) & maxit<10) {
    maxit <- maxit+1
    Parameters$ls50 <- log(exp(Parameters$ls50) * rlnorm(1, 0, 0.1))
    DoAssess <- Assess(data, Parameters = Parameters)

  }

  nm <- strsplit(data$Metadata$Value[1], ',') %>% unlist() %>% trimws() %>%
    paste(collapse = '_')
  nm <- paste0(nm, '.rdata')
  saveRDS(DoAssess, file.path('Results/Case_Studies/Assessments',nm))

}


# 6. Plot  Data  ----
assess_files <- list.files('Results/Case_Studies/Assessments')
assess_files <- assess_files[!assess_files=="Tetandara (Arubara).rdata"]

assess_list <- list()
for (i in seq_along(assess_files)) {
  assess_list[[i]] <- readRDS(file.path('Results/Case_Studies/Assessments', assess_files[i]))
}

## Plot F and SPR
F_SPR_list <- list()
for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data
  pred_FM <- assess_list[[i]]$rep$F_m
  pred_SPR <- assess_list[[i]]$rep$SPR

  F_SPR_list[[i]] <- data.frame(Site=data$Metadata$Value[1],
                                Date=lubridate::my(paste(data$Month, data$Year, sep="-")),
                                F=pred_FM,
                                SPR=pred_SPR)

}
F_SPR_df <- do.call('rbind', F_SPR_list)

p1 <- ggplot(F_SPR_df, aes(x=Date, y=F)) +
  facet_wrap(~Site, ncol=4, scales='free_y') +
  expand_limits(y=c(0,1)) +
  geom_line() + theme_bw() +
  labs(y='Fishing Mortality (F)') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        strip.text=element_text(size=8))
p1

ggsave('img/Case_Studies/F.png', p1, width=10, height=6)


p1 <- ggplot(F_SPR_df, aes(x=Date, y=SPR)) +
  facet_wrap(~Site, ncol=4) +
  expand_limits(y=c(0,1)) +
  geom_line() + theme_bw() +
  labs(y='Spawning Potential Ratio (SPR)') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        strip.text=element_text(size=8))
p1

ggsave('img/Case_Studies/SPR.png', p1, width=10, height=6)



## Plot Seasonal Recruitment
R0_list <- list()

for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data
  df <- get_seasonal_recruitment(assess_list[[i]])
  df$Site <- data$Metadata$Value[1]

  R0_list[[i]] <- df

}
R0_df <- do.call('rbind', R0_list)
R0_df$Month_Name <- factor(R0_df$Month_Name, levels=month.abb, ordered = TRUE)
R0_df$Month <- factor(R0_df$Month, levels=1:12, ordered = TRUE)

p1 <- ggplot(R0_df, aes(x=Month_Name, y=Recruitment, group=1)) +
  facet_wrap(~Site) +
  geom_line() +
  expand_limits(y=0) +
  theme_bw() +
  labs(x='Month', y='Relative Annual Recruitment') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        strip.text = element_text(size=6))

ggsave('img/Case_Studies/SeasonalR0.png', p1, width=8, height=5)

## Plot Selectivity
select_list <- list()
for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data
  df <- data.frame(Age=0:14, Selectivity=assess_list[[i]]$rep$selA)
  df$Site <- data$Metadata$Value[1]
  select_list[[i]] <- df
}
select_df <- do.call('rbind', select_list)

p1 <- ggplot(select_df, aes(x=Age, y=Selectivity)) +
  facet_wrap(~Site) +
  geom_line() +
  expand_limits(y=0) +
  theme_bw() +
  labs(x='Age', y='Selectivity') +
  theme(strip.text = element_text(size=6))


ggsave('img/Case_Studies/Select.png', p1, width=8, height=5)




## Plot CAW Data
for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data

  nbin <- length(data$Weight_Mids)
  df <- data.frame(Year=rep(data$Year, each=nbin),
                   Month=rep(data$Month, each=nbin),
                   Bin=data$Weight_Mids)

  df$interpolated <- data$iterpolated$interpolated

  df2 <- data.frame(Year=data$Year,
                    Month=data$Month,
                    Nsamp=round(data$CAW_nsamp,0))

  df$Weight <- as.vector(t(data$CAW))
  df$Fit <-  as.vector(assess_list[[i]]$rep$predCAW)

  df <- left_join(df, df2,by = join_by(Year, Month))

  ninter <- sum(df$interpolated)
  df$interpolated <- factor(df$interpolated, levels=c(FALSE, TRUE), ordered = TRUE)

  p <- ggplot(df) +
    facet_grid(Month~Year) +
    geom_bar(aes(x=Bin, y=Weight, fill=interpolated), stat='identity') +
    geom_line(aes(x=Bin, y=Fit)) +
    labs(title=data$Metadata$Value[1]) +
    geom_text(aes(x=-Inf, y=Inf, label=Nsamp, hjust = 0, vjust = 1)) +
    theme_bw() +
    labs(x='Weight Bin (kg)', y='Count', fill='Interpolated') +
    scale_fill_manual(values=c('black', 'grey')) +
    guides(fill='none')

  n_year <- df$Year %>% unique() %>% length()

  nm <- paste0('img/Case_Studies/CAW/', i, '.png')
  height <- 11
  width <- n_year * 3
  ggsave(nm, p, width=width, height=height)

}


# Make plot of overall weight comps
plist <- list()
for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data

  nbin <- length(data$Weight_Mids)
  df <- data.frame(Year=rep(data$Year, each=nbin),
                   Month=rep(data$Month, each=nbin),
                   Bin=data$Weight_Mids)

  df$interpolated <- data$iterpolated$interpolated

  df2 <- data.frame(Year=data$Year,
                    Month=data$Month,
                    Nsamp=round(data$CAW_nsamp,0))

  df$Weight <- as.vector(t(data$CAW))
  df$Fit <-  as.vector(assess_list[[i]]$rep$predCAW)

  df <- left_join(df, df2,by = join_by(Year, Month))

  ninter <- sum(df$interpolated)
  df$interpolated <- factor(df$interpolated, levels=c(FALSE, TRUE), ordered = TRUE)

  df <- df %>% group_by(Bin) %>% summarise(Weight=sum(Weight))
  df$Site <- data$Metadata$Value[1]
  df$Weight <- df$Weight/sum(df$Weight)
  plist[[i]] <- df
}
pdf <- do.call('rbind', plist)

p <- ggplot(pdf) +
  facet_wrap(~Site) +
  geom_bar(aes(x=Bin, y=Weight), stat='identity') +
    theme_bw() +
    labs(x='Weight Bin (kg)', y='Count') +
    guides(fill='none') +
  theme(strip.text = element_text(size=8))
p
ggsave('img/Case_Studies/CAW_overall.png', p, width=10, height=6)



## Plot Fit to Relative Catch
catch_list <- list()
for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data
  predC <- assess_list[[i]]$rep$predCB

  catch_list[[i]] <- data.frame(Site=data$Metadata$Value[1],
                                Date=lubridate::my(paste(data$Month, data$Year, sep="-")),
                                Observed=data$Catch_Mean,
                                Predicted=predC/mean(predC[!is.na(data$Catch_Mean)]))
}

catch_df <- do.call('rbind', catch_list) %>%
  tidyr::pivot_longer(., cols=3:4)

catch_df <- catch_df %>% filter(is.finite(value)==TRUE)
p1 <- ggplot(catch_df, aes(x=Date, y=value, color=name)) +
  facet_wrap(~Site, ncol=4) +
  expand_limits(y=0) +
  geom_line() + theme_bw() +
  labs(color='Legend', x='Date', y='Relative Catch') +
  theme(strip.text=element_text(size=6))

ggsave('img/Case_Studies/Catch.png', p1, width=12, height=7)


## Plot Fit to Index
index_list <- list()
for (i in seq_along(assess_list)) {
  data <- assess_list[[i]]$Data
  predIndex <- assess_list[[i]]$rep$stpredIndex

  index_list[[i]] <- data.frame(Site=data$Metadata$Value[1],
                                Date=lubridate::my(paste(data$Month, data$Year, sep="-")),
                                Observed=data$Index_Mean,
                                Predicted=predIndex)
}

index_df <- do.call('rbind', index_list) %>%
  tidyr::pivot_longer(., cols=3:4)

index_df <- index_df %>% filter(is.finite(value)==TRUE)
p1 <- ggplot(index_df, aes(x=Date, y=value, color=name)) +
  facet_wrap(~Site) +
  expand_limits(y=0) +
  geom_line() + theme_bw() +
  labs(color='Legend', x='Date', y='Index of Abundance')

ggsave('img/Case_Studies/Index.png', p1, width=12, height=5)




utilpow <- 0.4
ref_list <- list()

for (i in seq_along(assess_list)) {
  assess <-  assess_list[[i]]

  season_R0 <- get_seasonal_recruitment(assess)

  Est <- calculate_optimal_fishing(R0_m=season_R0$Recruitment,
                                   steepness=assess$Data$h,
                                   Weight_Age_Mean=assess$Data$Weight_Age_Mean,
                                   Maturity_at_Age=assess$Data$Maturity_at_Age,
                                   M_at_Age=assess$Data$M_at_Age,
                                   Post_Spawning_Mortality=assess$Data$Post_Spawning_Mortality,
                                   sel_at_age=assess$rep$selA,
                                   opt_type=1, utilpow=utilpow)


  average_monthly_F <- calc_average_monthly_F(assess)
  average_monthly_F$Optimal <- Est$F_m
  average_monthly_F$Variable <- "F"
  average_monthly_F <- average_monthly_F %>%
    tidyr::pivot_longer(., cols=c(2, 4), names_to = 'Name', values_to = 'Value')


  average_monthly_SPR <- calc_average_monthly_SPR(assess)
  average_monthly_SPR$Optimal <- Est$SPR
  average_monthly_SPR$Variable <- "SPR"
  average_monthly_SPR <- average_monthly_SPR %>%
    tidyr::pivot_longer(., cols=c(2, 4), names_to = 'Name', values_to = 'Value')

  df <- bind_rows(average_monthly_F, average_monthly_SPR)

  df$Site <- assess_list[[i]]$Data$Metadata$Value[1]
  ref_list[[i]] <- df

}

ref_df <- do.call('rbind', ref_list)
ref_df$Month_Name <- factor(ref_df$Month_Name, levels=month.abb, ordered = TRUE)
ref_df$Month <- factor(ref_df$Month, levels=1:12, ordered = TRUE)

p1 <- ggplot(ref_df , aes(x=Month_Name, y=Value, color=Name, group=Name)) +
  facet_grid(Variable~Site) +
  geom_line() +
  expand_limits(y=0) +
  theme_bw() +
  labs(x='Month', color='Legend') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        strip.text = element_text(size=8))

p1 <- ggplot(ref_df %>% filter(Variable=='F') , aes(x=Month_Name, y=Value, color=Name, group=Name)) +
  facet_wrap(~Site) +
  geom_line() +
  expand_limits(y=0) +
  theme_bw() +
  labs(x='Month', y='Fishing Mortality', color='Legend') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        strip.text = element_text(size=6))

ggsave('img/Case_Studies/F_Refs.png', p1, width=10, height=8)


p1 <- ggplot(ref_df %>% filter(Variable=='SPR') , aes(x=Month_Name, y=Value, color=Name, group=Name)) +
  facet_wrap(~Site) +
  geom_line() +
  expand_limits(y=0) +
  theme_bw() +
  labs(x='Month', y='Spawning Potential Ratio', color='Legend') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        strip.text = element_text(size=6))

ggsave('img/Case_Studies/SPR_Refs.png', p1, width=10, height=8)



