
library(ggplot2)
library(dplyr)
library(lubridate)
library(ggthemes)

source('R/functions.r')

# ---- Indonesia Case Study Data ----
Data <- read.csv('G:/My Drive/1_PROJECTS/Biospherics_Indo_Octopus/Data/Forkani_Darawa_Octopus_Master.xlsx - Raw Data.csv')

# --- Format Dates ---
Data_mod <- Data
# Format Dates
Data_mod$Date <- lubridate::mdy(Data$Date)
Data_mod$Time.Fishing <- suppressWarnings(lubridate::hm(Data$Time.Fishing)) %>%
  as.numeric() / 3600 # in hours
Data_mod$Month_Name <- month.abb[Data_mod$Month]

Years <- rep(min(Data_mod$Year):max(Data_mod$Year), each=12)
Months <- 1:12
Dates_DF <- data.frame(Year=Years, Month=Months)
Dates_DF$Year_Month <- paste(Dates_DF$Year, Dates_DF$Month, sep="_")
Dates_DF$Date <- lubridate::my(paste(Dates_DF$Month, Dates_DF$Year, sep="-"))

# --- Meta-data ---
date.range <- Data_mod$Date %>% range()
n.collector <- Data_mod$Data.collector %>% tolower() %>% unique() %>% length()
n.fisher <- Data_mod$Fisher %>% tolower() %>% unique() %>% length()
n.sites <- Data_mod$Site1 %>% tolower() %>% unique() %>% length()
data.frame(First_Date=date.range[1], Last_Date=date.range[2],
           N.collector=n.collector,
           N.fisher=n.fisher,
           N.sites=n.sites)

Data_mod$Site1 %>% tolower() %>% unique() %>% sort()

# --- Effort Data ----
Effort_Month_Fisher <- Data_mod %>%
  group_by(Year, Month_Name, Month, Fisher) %>%
  distinct(Year, Month, Day, Fisher, Time.In, Time.Fishing) %>%
  group_by(Year, Month_Name, Month, Fisher) %>%
  summarise(Effort=sum(Time.Fishing, na.rm=TRUE), .groups='keep')

Effort <- Effort_Month_Fisher %>% group_by(Year, Month, Month_Name) %>%
  filter(is.na(Effort)==FALSE) %>%
  summarise(Effort=sum(Effort), .groups='keep')
Effort$Date <- lubridate::my(paste(Effort$Month, Effort$Year, sep="-"))

Effort <- left_join(Dates_DF, Effort)
Effort <- assignGroups(Effort, 'Effort')
Effort <- dropLastNAs(Effort, 'Effort')

ggplot(Effort, aes(x=Date, y=Effort, group=groups)) +
  geom_point() +
  stat_summary(fun=sum, geom="line") +
  labs(x="", y='Total Effort (hours)') +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  coord_cartesian(clip = 'off') +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_blank())
ggsave('Figures/Data/Effort.png', width=6)

Effort %>% arrange(Year, Month_Name) %>% filter(Year==2017)
Effort %>% arrange(Year, Month_Name) %>% filter(Year==2018)
Effort %>% arrange(Year, Month_Name) %>% filter(Year==2019)
Effort %>% arrange(Year, Month_Name) %>% filter(Year==2020)
Effort %>% arrange(Year, Month_Name) %>% filter(Year==2021)
Effort %>% arrange(Year, Month_Name) %>% filter(Year==2022)

# ---- Catch Data ----

Catch_Month_Fisher <- Data_mod %>%
  group_by(Year, Month, Month_Name, Fisher) %>%
  summarize(Catch=sum(Individual.Weight), .groups='keep')

Catch_Month <- Data_mod %>%
  group_by(Year, Month, Month_Name) %>%
  summarize(Catch=sum(Individual.Weight, na.rm=TRUE), .groups='keep')

Catch_Month <- left_join(Dates_DF, Catch_Month)
Catch_Month <- assignGroups(Catch_Month, 'Catch')
Catch_Month <- dropLastNAs(Catch_Month, 'Catch')

ggplot(Catch_Month, aes(x=Date, y=Catch, group=groups)) +
  geom_point() +
  stat_summary(fun=sum, geom="line") +
  labs(x="", y='Catch (kg)') +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  coord_cartesian(clip = 'off') +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.background=element_blank())

ggsave('Figures/Data/TotalCatch.png', width=6)

# ---- CPUE -----

CPUE_Month_Fisher <- left_join(Catch_Month_Fisher, Effort_Month_Fisher,
                               by=c("Year", "Month", "Fisher",'Month_Name'))

CPUE_Month_Fisher <- CPUE_Month_Fisher %>% filter(Effort>0)


CPUE_Month <- CPUE_Month_Fisher %>%
  mutate(CPUE=Catch/Effort, lCPUE=log(CPUE))

CPUE_Month$Year <- factor(CPUE_Month$Year, ordered = TRUE)
CPUE_Month$Month <- factor(CPUE_Month$Month, ordered = TRUE)

CPUE_Month$Year_Month <- paste(CPUE_Month$Year, CPUE_Month$Month, sep="_")

CPUE_Month$Year_Month <- factor(CPUE_Month$Year_Month,
                                levels=CPUE_Month$Year_Month %>% unique())

CPUE_Month <- CPUE_Month %>% filter(is.na(CPUE)==FALSE)

mod <- lm(lCPUE~Fisher+Year_Month, data=CPUE_Month)

CPUE_Month$Estimated_CPUE <- exp(predict(mod)+sd(summary(mod)$residuals)^2/2)
CPUE_Month$Est_CPUE975 <- exp(predict(mod)+predict(mod, se.fit=T)$se*1.96+sd(summary(mod)$residuals)^2/2)
CPUE_Month$Est_CPUE025 <- exp(predict(mod)-predict(mod, se.fit=T)$se*1.96+sd(summary(mod)$residuals)^2/2)

Est_CPUE <- aggregate(CPUE_Month$Estimated_CPUE~ CPUE_Month$Year_Month, FUN=mean, na.rm=T)
Est_CPUE025 <- aggregate(CPUE_Month$Est_CPUE025~ CPUE_Month$Year_Month, FUN=mean, na.rm=T)
Est_CPUE975 <- aggregate(CPUE_Month$Est_CPUE975~ CPUE_Month$Year_Month, FUN=mean, na.rm=T)
SD <- predict(mod, se.fit=T)
logSD <- aggregate(SD$se.fit~ CPUE_Month$Year_Month, FUN=mean, na.rm=T)
DFa <- data.frame(Year_Month=Est_CPUE$`CPUE_Month$Year_Month`,
                  mean=Est_CPUE$`CPUE_Month$Estimated_CPUE`,
                  lower=Est_CPUE025$`CPUE_Month$Est_CPUE025`,
                  upper=Est_CPUE975$`CPUE_Month$Est_CPUE975`,
                  logSD=logSD$`SD$se.fit`)


DF2 <- left_join(Dates_DF, DFa, by="Year_Month")
DF2 <- assignGroups(DF2)
DF2 <- dropLastNAs(DF2)
CPUE_DF <- DF2

ggplot(DF2, aes(x=Date, y=mean, group=groups)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='lightgray') +
  geom_point() +
  stat_summary(fun=mean, geom="line") +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b-%Y",
               expand=c(0,0)) +
  coord_cartesian(clip = 'off') +
  labs(x='', y='Mean CPUE') +
  theme_clean() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        plot.background=element_blank())
ggsave('Figures/Data/CPUE.png', width=6)




# ---- Weight Composition Data -----
SizeData <- Data_mod %>% select(Year, Month, Day, Month_Name, Sex, Mantle.Length, Individual.Weight,
                                Fisher, Data.collector) %>%
  mutate(Length=as.numeric(Mantle.Length),
         Weight=as.numeric(Individual.Weight)) %>%
  filter(is.na(Length)==FALSE, is.na(Weight)==FALSE,
         Weight>0, Length>0)

SizeData$Month_Name <- factor(SizeData$Month_Name, ordered = TRUE, levels = month.abb[1:12])

SizeData$Sex[SizeData$Sex == 'Betina'] <- 'Female'
SizeData$Sex[SizeData$Sex == 'Jantan'] <- 'Male'
SizeData$Sex[SizeData$Sex == 'Tidak Tercatat'] <- 'Not Recorded'

ggplot(SizeData, aes(x=Length, y=Weight)) +
  facet_grid(~Sex) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())

ggsave('Figures/Data/Length_Weight.png', width=6)

ggplot(SizeData, aes(x=Length, y=Weight, color=Data.collector)) +
  facet_grid(Year~Month_Name) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())
ggsave('Figures/Data/Length_Weight_month_year_collector.png', width=9, height=7)


# apply correction to weight records
prob <- SizeData %>% filter(Year %in% 2021, Month %in% c(1:3,8))
ggplot(prob, aes(x=Length, y=Weight, color=Data.collector)) +
  facet_grid(Year~Month_Name) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())
ggsave('Figures/Data/Length_Weight_month_year_collector_sub.png', width=6, height=4)

# January
ggplot(prob %>% filter(Month==1), aes(x=Length, y=Weight, color=Data.collector)) +
  facet_grid(Year~Month_Name) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())

Sub <- prob %>% filter(Month==1, Length<9, Weight>1.5)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()

SizeData_correct <- SizeData
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10

Sub <- prob %>% filter(Month==1, Length<13, Weight>2.5)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10

Sub <- prob %>% filter(Month==1, Length<5, Weight>1)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10


# February
ggplot(prob %>% filter(Month==2), aes(x=Length, y=Weight, color=Data.collector)) +
  facet_grid(Year~Month_Name) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())

Sub <- prob %>% filter(Month==2, Length<5, Weight>1.5)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10

Sub <- prob %>% filter(Month==2, Length<=11, Weight>2)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10


# March
ggplot(prob %>% filter(Month==3), aes(x=Length, y=Weight, color=Data.collector)) +
  facet_grid(Year~Month_Name) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())

Sub <- prob %>% filter(Month==3, Length<=11, Weight>=2)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10


# August
ggplot(prob %>% filter(Month==8), aes(x=Length, y=Weight, color=Data.collector)) +
  facet_grid(Year~Month_Name) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())

Sub <- prob %>% filter(Month==8, Length==10, Weight>3)
rows <- plyr::match_df(SizeData,Sub) %>% rownames()
SizeData_correct[rownames(SizeData_correct)%in%rows,'Weight'] <- SizeData[rownames(SizeData)%in%rows,'Weight']/10


ggplot(SizeData_correct, aes(x=Length, y=Weight)) +
  facet_grid(~Sex) +
  geom_point(alpha=0.5)+
  labs(x='Length (cm)', y='Weight (kg)') +
  theme_clean() +
  theme(plot.background=element_blank())

ggsave('Figures/Data/Length_Weight_correct.png', width=6)




Weight_Range <- range(SizeData_correct$Weight)
hist(SizeData_correct$Weight)
by <- 0.1
Weight_Bins <- seq(0, to=Weight_Range[2], by=by)
Weight_Mids <- seq(0.5*by, by=by, length.out=length(Weight_Bins)-1)

WeightData_Binned <- SizeData_correct %>% group_by(Year, Month, Month_Name, Sex) %>%
  summarise(bin = cut(Weight, breaks = Weight_Bins), .groups = 'keep') %>%
  summarise(Count=as.numeric(table(bin)), Bin=Weight_Mids, .groups = 'keep')

maxWeightBin <- 2
Weight_Mids <- Weight_Mids[Weight_Mids<=maxWeightBin]
Weight_Bins <- Weight_Bins[1:(length(Weight_Mids)+1)]

WeightData_Binned_2 <- WeightData_Binned %>% group_by(Year, Month, Sex) %>%
  mutate(n2=sum(Count[Bin>=maxWeightBin])) %>%
  mutate(Count2=replace(Count, Bin>=maxWeightBin, n2[Bin>=maxWeightBin])) %>%
  filter(Bin<=maxWeightBin)

WeightData_Binned_2$Month_Name <- factor(WeightData_Binned_2$Month_Name, ordered = TRUE, levels = month.abb[1:12])
WeightData_Binned_2$Year <- factor(WeightData_Binned_2$Year)

WeightData_Binned_2 <- WeightData_Binned_2 %>% group_by(Year, Month_Name, Sex) %>%
  mutate(n=sum(Count), Count=Count/sum(Count))

text <- WeightData_Binned_2 %>% distinct(Year, Month_Name, Sex, n)

WeightData_Binned_2$BinLabel <- WeightData_Binned_2$Bin
WeightData_Binned_2$BinLabel <- factor(WeightData_Binned_2$BinLabel)
levels(WeightData_Binned_2$BinLabel)[length(levels(WeightData_Binned_2$BinLabel))] <- paste0(levels(WeightData_Binned_2$BinLabel)[length(levels(WeightData_Binned_2$BinLabel))], "+")

ggplot(WeightData_Binned_2 %>% filter(Sex=='Female'),
       aes(x=BinLabel, y=Count)) +
  facet_grid(Month_Name~Year) +
  geom_bar(stat='identity') +
  theme_clean() +
  geom_text(data=text%>%filter(Sex=='Female'),
            aes(x=-Inf, y=Inf, label=n, hjust=0, vjust=1)) +
  labs(x='Weight Classes (kg)', y='Relative Frequency') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(plot.background=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))
ggsave('Figures/Data/Weight_comp_Female.png', width=9, height=7)

ggplot(WeightData_Binned_2 %>% filter(Sex=='Male'), aes(x=BinLabel, y=Count)) +
  facet_grid(Month_Name~Year) +
  geom_bar(stat='identity') +
  theme_clean() +
  geom_text(data=text%>%filter(Sex=='Male'),
            aes(x=-Inf, y=Inf, label=n, hjust=0, vjust=1)) +
  labs(x='Weight Classes (kg)', y='Relative Frequency') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(plot.background=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))
ggsave('Figures/Data/Weight_comp_Male.png', width=9, height=7)


ggplot(WeightData_Binned_2 %>% filter(Sex=='Not Recorded'), aes(x=BinLabel, y=Count)) +
  facet_grid(~Month_Name) +
  geom_bar(stat='identity') +
  theme_clean() +
  geom_text(data=text%>%filter(Sex=='Not Recorded'),
            aes(x=-Inf, y=Inf, label=n, hjust=0, vjust=1)) +
  labs(x='Weight Classes (kg)', y='Relative Frequency') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(plot.background=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))


text2 <- text %>% group_by(Year, Month_Name) %>%
  summarize(n=sum(n), Sex=NA)

ggplot(WeightData_Binned_2 %>% filter(Sex!='Not Recorded'), aes(x=BinLabel, y=Count, fill=Sex)) +
  facet_grid(Month_Name~Year) +
  geom_bar(stat='identity') +
  theme_clean() +
  geom_text(data=text2,
            aes(x=-Inf, y=Inf, label=n, hjust=0, vjust=1),size=3) +
  labs(x='Weight Classes (kg)', y='Relative Frequency') +
  theme(plot.background=element_blank(),
        legend.position="bottom")


# Length Comp
Length_Range <- range(SizeData_correct$Length)
by <- 1
Length_Bins <- seq(0, to=Length_Range[2], by=by)
Length_Mids <- seq(0.5*by, by=by, length.out=length(Length_Bins)-1)

LengthData_Binned <- SizeData_correct %>% group_by(Year, Month, Month_Name, Sex) %>%
  summarise(bin = cut(Length, breaks = Length_Bins), .groups = 'keep') %>%
  summarise(Count=as.numeric(table(bin)), Bin=Length_Mids, .groups = 'keep')


maxLengthBin <- 15
Length_Mids <- Length_Mids[Length_Mids<=maxLengthBin]
Length_Bins <- Length_Bins[1:(length(Length_Mids)+1)]

LengthData_Binned_2 <- LengthData_Binned %>% group_by(Year, Month, Sex) %>%
  mutate(n2=sum(Count[Bin>=maxLengthBin])) %>%
  mutate(Count2=replace(Count, Bin>=maxLengthBin, n2[Bin>=maxLengthBin])) %>%
  filter(Bin<=maxLengthBin)


LengthData_Binned_2$Month_Name <- factor(LengthData_Binned_2$Month_Name, ordered = TRUE, levels = month.abb[1:12])
LengthData_Binned_2$Year <- factor(LengthData_Binned_2$Year)

LengthData_Binned_2 <- LengthData_Binned_2 %>% group_by(Year, Month_Name, Sex) %>%
  mutate(n=sum(Count), Count=Count/sum(Count))

text <- LengthData_Binned_2 %>% distinct(Year, Month_Name, Sex, n)

LengthData_Binned_2$BinLabel <- LengthData_Binned_2$Bin
LengthData_Binned_2$BinLabel <- factor(LengthData_Binned_2$BinLabel)
levels(LengthData_Binned_2$BinLabel)[length(levels(LengthData_Binned_2$BinLabel))] <-
  paste0(levels(LengthData_Binned_2$BinLabel)[length(levels(LengthData_Binned_2$BinLabel))], "+")

ggplot(LengthData_Binned_2 %>% filter(Sex=='Female'),
       aes(x=BinLabel, y=Count)) +
  facet_grid(Month_Name~Year) +
  geom_bar(stat='identity') +
  theme_clean() +
  geom_text(data=text%>%filter(Sex=='Female'),
            aes(x=-Inf, y=Inf, label=n, hjust=0, vjust=1)) +
  labs(x='Length Classes (cm)', y='Relative Frequency') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(plot.background=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))
ggsave('Figures/Data/Length_comp_Female.png', width=9, height=7)


ggplot(LengthData_Binned_2 %>% filter(Sex=='Male'),
       aes(x=BinLabel, y=Count)) +
  facet_grid(Month_Name~Year) +
  geom_bar(stat='identity') +
  theme_clean() +
  geom_text(data=text%>%filter(Sex=='Male'),
            aes(x=-Inf, y=Inf, label=n, hjust=0, vjust=1)) +
  labs(x='Length Classes (cm)', y='Relative Frequency') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(plot.background=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))
ggsave('Figures/Data/Length_comp_Male.png', width=9, height=7)


# ---- Make List for Data used in Assessment ----

data <- list()

maxage <- 14
Ages <- 0:maxage
n_ages <- length(Ages)

# ---- Weight-at-Age ----

# Herwig et al. (2012)
# power function
# digitized points from Fig 1 with https://automeris.io/WebPlotDigitizer/

age_weight <- read.csv('data-raw/age_weight.csv', header=FALSE)
names(age_weight) <- c('Age', 'Weight')

model <- lm(log(Weight)~log(Age), data=age_weight)

lm_eqn <- function(df){
  m <- model
  eq <- bquote(italic(W) == .(a) %.% italic(A) ^ .(b)* ","~~italic(r)^2~"="~.(r2),
                   list(a = format(unname(exp(coef(m))[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}

pwr <- function(x)
  exp(predict(model, newdata = data.frame(Age=x)))

ggplot(age_weight, aes(x=Age, y=Weight)) +
  geom_point() +
  stat_function(fun = pwr) +
  geom_text(x = 150, y = 2, label = lm_eqn(age_weight), parse = TRUE) +
  labs(x="Age (days)", y="Weight (kg)") +
  theme_clean() +
  theme(plot.background=element_blank())
ggsave('Figures/LifeHistory/Weight_Age.png', width=4, height=4)

# calculate weight-at-age in months
s <- summary(model)
sd <- sd(s$residuals)
df <- data.frame(Age_month=0:maxage,
                 age_day=seq(0, by=30.41667, length.out=maxage+1)
                 )
df$meanW <- pwr(x=df$age_day)
# 95% CI
df$lower <- exp(log(df$meanW)-0.5*sd^2 - 1.96*sd)
df$upper <- exp(log(df$meanW)-0.5*sd^2 + 1.96*sd)

ggplot(df, aes(x=Age_month, y=meanW)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill='lightgray') +
  geom_line() +
  coord_cartesian(ylim=c(0, 4)) +
  geom_text(x = 150, y = 2, label = lm_eqn(age_weight), parse = TRUE) +
  labs(x="Age (month)", y="Weight (kg)") +
  theme_clean() +
  theme(plot.background=element_blank())

data$Weight_Age <- df$meanW
data$Weight_Age_SD <- rep(sd, length(data$Weight_Age))


# ---- Length-at-Age ----

length_weight <- SizeData_correct %>% select(Weight, Length)
mod <- lm(log(Length)~log(Weight), data=length_weight)
summary(mod)
coef(mod)
exp(coef(mod)[1])

pwr2 <- function(x)
  exp(predict(mod, newdata = data.frame(Weight=x), type='response'))

ggplot(length_weight, aes(x=Weight, y=Length)) +
  geom_point(alpha=0.1) +
  stat_function(fun = pwr2)

pred <- pwr2(data$Weight_Age)


# ---- Mortality-at-Age ----

data$M_at_Age <- rep(0.1, n_ages)

# ---- Maturity-at-Age ----

data$Mat_at_Age <- c(rep(0,11), 0.05, 0.95, rep(1, n_ages-13))

df_mat <- data.frame(Age=Ages, Spawn=data$Mat_at_Age )

ggplot(df_mat, aes(x=Age, y=Spawn)) +
  geom_line() +
  labs(x="Age (month)", y="Probability of Spawning") +
  theme_clean() +
  theme(plot.background=element_blank())
ggsave('Figures/LifeHistory/Prob_Spawning.png', width=4, height=4)


# ---- Post-Spawning-Mortality-at-Age ----
data$PSM_at_Age <- data$Mat_at_Age

# ---- Steepness -----
CR <- exp(0.2*log(2.4E5)) # Van Heukelem 1973 & Hilborn and Walters 2021
h <- CR2h(CR)
data$h <- round(h,2)


# ---- Relative Effort Time-Series ----
data$Effort <- Effort$Effort/mean(Effort$Effort, na.rm=TRUE)
data$Effort_SD <- rep(0.2, length(data$Effort))

# ---- CPUE Time-Series ----
data$CPUE <- CPUE_DF$mean/mean(CPUE_DF$mean, na.rm=TRUE)
data$CPUE_SD <- CPUE_DF$logSD

# ---- Weight Composition ----
data$WghtBins <- Weight_Bins
data$WghtMids <- Weight_Mids
nbins <- length(data$WghtMids)
nmonths <- length(data$Effort)

CAW <- matrix(0, nrow=nbins, ncol=nmonths)
for (i in 1:nrow(Effort)) {
  y <- Effort$Year[i]
  m <- Effort$Month[i]
  w <- WeightData_Binned_2 %>% filter(Year==y, Month==m) %>%
    group_by(Bin) %>%
    summarize(n=sum(Count2))
  if (length(w$n)>0)
    CAW[,i] <- w$n
}

CAW_ESS <- apply(CAW, 2, sum)
CAW_ESS[CAW_ESS>200] <- 200

data$CAW <- CAW
data$CAW_ESS <- CAW_ESS


# ---- Length Composition ----
# no used

# ---- Save Data Object ----
saveRDS(data, 'data-raw/casestudydata.rds')

