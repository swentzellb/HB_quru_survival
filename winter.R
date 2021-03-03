# Sage Wentzell-Brehme
# January 31, 2021

# Looking at trends in winter climate data at Hubbard Brook Experimental Forest

# load libraries
library(tidyverse)
library(lubridate)
library(lme4)
library(ggplot2)


# set working directory
setwd("~/Thesis/HB_quru_survival")

#Hubbard Brook temperature data
HB_temp <- read.csv("data/HBEF_air_temp_daily_1957-2019.csv")

#snow & soil frost data
HB_snow <- read.csv("data/HBEF_snowcourse_1956-2020.csv")

# Clean missing data values to NA (-99 = missing)
HB_snow$snow_depth[HB_snow$snow_depth==-99] <- NA
HB_snow$swe[HB_snow$swe==-99] <- NA
HB_snow$frost_depth[HB_snow$frost_depth==-99] <- NA
HB_snow$frost_pct[HB_snow$frost_pct==-99] <- NA

#check that NA values are correct
summary(HB_snow)

#subset the data to the HQ site & 2011-2020 winter time period
HB_snowHQ <- HB_snow %>%
  filter(Site == "STAHQ", 
         Date >= as.Date("2011-06-06"))

#subset to Watershed 1
HB_snowW1 <- HB_snow %>%
  filter(Site %in% c("STA1", "STA2", "STA3"), 
         Date >= as.Date("2011-06-06"))

#subset to Watershed 6
HB_snowW6 <- HB_snow %>%
  filter(Site %in% c("STA9", "STA10", "STA11"), 
         Date >= as.Date("2011-06-06"))



table(HB_snow$Site)


#another supplemental table showing # of data points in each yr
# FIX THIS
Fig1 <- HB_snowHQ %>%
   select(WINTER, snow_depth) %>%
   filter(na.rm=TRUE) %>%
   group_by(WINTER) %>%
   summarise(count = n())

Fig2 <- HB_snowW1 %>%
  select(WINTER, snow_depth) %>%
  filter(na.rm=TRUE) %>%
  group_by(WINTER) %>%
  summarise(count = n())


snowdepthsum <- HB_snowHQ %>%
  select(WINTER, snow_depth) %>%
  filter(na.rm=TRUE) %>%
  group_by(WINTER) %>%
  summarise(count = n())

frostdepthsum <- HB_snowHQ %>%
  select(WINTER, frost_depth) %>%
  filter(na.rm=TRUE) %>%
  group_by(WINTER) %>%
  summarise(count = n())

frostpctsum <- HB_snowHQ %>%
  select(WINTER, frost_pct) %>%
  filter(na.rm=TRUE) %>%
  group_by(WINTER) %>%
  summarise(count = n())



#plot snow depth by frost depth
ggplot(data = HB_snow)+
  geom_point(mapping=aes(x=snow_depth, y=frost_depth))

#summarize climate variables for each winter period at HQ site
sumHQ <- HB_snowHQ %>%
  dplyr::group_by(WINTER)%>%
  summarize(medfrostdepth = median(frost_depth, na.rm=TRUE),
            sumfrostdepth = sum(frost_depth>=50, na.rm=TRUE),
            sumsnowdepth = sum(snow_depth>=200, na.rm=TRUE), 
            medsnowdepth = median(snow_depth, na.rm=TRUE),
            meansnowdepth = mean(snow_depth, na.rm=TRUE),
            sdsnowdepth = sd(snow_depth, na.rm=TRUE),
            medfrostpct = median(frost_pct, na.rm=TRUE))

#summarize climate variables for each winter period at Watershed 1
sumW1 <- HB_snowW1 %>%
  dplyr::group_by(WINTER)%>%
  summarize(medfrostdepth = median(frost_depth, na.rm=TRUE),
            sumfrostdepth = sum(frost_depth>=50, na.rm=TRUE),
            sumsnowdepth = sum(snow_depth>=200, na.rm=TRUE), 
            medsnowdepth = median(snow_depth, na.rm=TRUE),
            meansnowdepth = mean(snow_depth, na.rm=TRUE),
            sdsnowdepth = sd(snow_depth, na.rm=TRUE),
            medfrostpct = median(frost_pct, na.rm=TRUE))

#summarize climate variables for each winter period at Watershed 6
sumW6 <- HB_snowW6 %>%
  dplyr::group_by(WINTER)%>%
  summarize(medfrostdepth = median(frost_depth, na.rm=TRUE),
            sumfrostdepth = sum(frost_depth>=50, na.rm=TRUE),
            sumsnowdepth = sum(snow_depth>=200, na.rm=TRUE), 
            medsnowdepth = median(snow_depth, na.rm=TRUE),
            meansnowdepth = mean(snow_depth, na.rm=TRUE),
            sdsnowdepth = sd(snow_depth, na.rm=TRUE))

#check correlation of summary variables between HQ and Watershed 1 sampling sites
cor.test(sumHQ$medsnowdepth, sumW1$medsnowdepth)
cor.test(sumHQ$meansnowdepth, sumW1$meansnowdepth)
cor.test(sumHQ$medfrostdepth, sumW1$medfrostdepth)
cor.test(sumHQ$sumsnowdepth, sumW1$sumsnowdepth)
cor.test(sumHQ$sumfrostdepth, sumW1$sumfrostdepth)
cor.test(sumHQ$medfrostpct, sumW1$medfrostpct)

#check correlation of summary variables between HQ and Watershed 6 sampling sites
cor.test(sumHQ$medsnowdepth, sumW6$medsnowdepth)
cor.test(sumHQ$meansnowdepth, sumW6$meansnowdepth)
cor.test(sumHQ$medfrostdepth, sumW6$medfrostdepth)
cor.test(sumHQ$sumsnowdepth, sumW6$sumsnowdepth)
cor.test(sumHQ$sumfrostdepth, sumW6$sumfrostdepth)




#plot the median frost depth by winter
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=medfrostdepth))

#plot the median snow depth by winter
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=medsnowdepth))

#plot the sd snow depth by winter
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=sdsnowdepth))

#plot distribution of snow depth by winter
boxplot(snow_depth ~ WINTER, data = HB_snowHQ)

# linear regression of 
lm_1 <- lm(meansnowdepth ~ WINTER, data=sumHQ)
summary(lm_1)
glance(lm_1)

table(HB_snow$Date)

# look at varying winter lengths
# test <- HB_snow %>%
#   group_by(WINTER)%>%
#   summarize(winterstart = min(Date, na.rm=TRUE),
#             winterend = max(Date, na.rm=TRUE))


################################################
################################################
# add in random effect of survival data

# add in additional column for singular year rather than interval for WINTER
random$WINTER <- seq(2013,2020)
#
ran <- random %>%
  mutate(WINTER = seq(2013,2020)) %>%
  dplyr::select(WINTER, condval, condsd)

#join the random effects data table by the summarized winter climate variables from the HQ site
ranef_winter <- left_join(random, sumHQ, by = "WINTER")

#summarize climate variables at HQ site from 2012 to 2020
sumHQtotal <- sumHQ %>%
  summarize(medfrostdepth = median(medfrostdepth, na.rm=TRUE),
            medsumfrostdepth = median(sumfrostdepth, na.rm=TRUE),
            medsumsnowdepth = median(sumsnowdepth, na.rm=TRUE), 
            medsnowdepth = median(medsnowdepth, na.rm=TRUE))


##############################################
#plots of winter climate variables by year
##############################################

#plot of number of days with frost depth > 50 mm
ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=WINTER, y=sumfrostdepth))+
  geom_hline(yintercept=7, linetype="dotted")+
  theme_bw()

#plot of median frost depth by year
ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=WINTER, y=medfrostdepth))+
  theme_bw()

ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=WINTER, y=medsnowdepth))+
  geom_hline(yintercept=182.88, linetype="dotted")+
  theme_bw()

#plot of mean snow depth in mm by year
ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=WINTER, y=meansnowdepth))+
  theme_bw()+
  ylim(0,300)

#plot of number of days with snow depth > 200 mm by year
ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=WINTER, y=sumsnowdepth))+
  geom_hline(yintercept=9, linetype="dotted")+
  theme_bw()


################################################
# plots of ranef vs winter climate variables
################################################

ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=medfrostdepth, y=condval))+
  theme_bw()+
  ylab("random effect values")

ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=medsnowdepth, y=condval))+
  theme_bw()+
  ylab("random effect values")+
  xlab("median snow depth (mm)")

ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=sumsnowdepth, y=condval))+
  theme_bw()+
  ylab("random effect values")

ggplot(data = ranef_winter)+
  geom_point(mapping=aes(x=sumfrostdepth, y=condval))+
  theme_bw()+
  ylab("random effect values")

###########################################################################
#correlation between winter climate variables and random effect of survival
cor(ranef_winter$condval, ranef_winter$medsnowdepth)
cor(ranef_winter$condval, ranef_winter$sumfrostdepth)



###########################################################
#linear models of winter climate variables effects on ranef

# linear regression of ranef values by median frost depth 
lm1 <- lm(condval ~ medfrostdepth, data=ranef_winter)
summary(lm1)

# linear regression of ranef values by median snow depth 
lm2 <- lm(condval ~ medsnowdepth, data=ranef_winter)
summary(lm2)

# linear regression of ranef values by sum of snow depth > 200 mm
lm3 <- lm(condval ~ sumsnowdepth, data=ranef_winter)
summary(lm3)

# linear regression of ranef values by sum of frost depth > 50 mm
lm4 <- lm(condval ~ sumfrostdepth, data=ranef_winter)
summary(lm4)

# none of the winter climate variables are significantly impacting the ranef of survival 
















#######################################################
### Temperature Data


#summarize the data
summary(HB_temp)
table(HB_temp$STA)

#assign the date column to a date variable format
HB_temp$date <- as.Date(HB_temp$date)

#add a month column
HB_temp <- mutate(HB_temp, month = format(date, "%m"), 
                  year = year(date))

#subset the data to only look at the winter months - Dec, Jan, Feb
list_of_values <- c("01", "02", "12")
HBwinter <- HB_temp %>%
  filter(month %in% list_of_values)


#filter to just the temp data from the HQ or Headquarters plot
#create a new column for winter (ex: Dec 2018-Feb 2019 = winter 2019)
HBwinter <-  HBwinter %>%
  mutate(HBwinter, winter = ifelse(month == "12", year + 1, year))%>%
  filter(STA=="HQ")

#HBwinter$winter <- as.numeric(HBwinter$winter)
#HBwinter$winter[HBwinter$month== "12"] <- HBwinter$year + 1
#  mutate(HBwinter, winter = year)

#summarize the data by winter to look at the number of days that stay below freezing
#or have a maximum temp below 0 degrees Celsius
cold <- HBwinter %>%
  group_by(winter)%>%
  summarize(below = sum(MAX<=0))


#summarize data by winter to look at average temp and variation
#including average, maximum, and minimum daily temps
HBsum <- HBwinter %>%
  group_by(winter) %>%
  summarize(mean = mean(AVE),
            sdave = sd(AVE),
            minmean = mean(MIN), 
            sdmin = sd(MIN), 
            maxmean = mean(MAX), 
            sdmax = sd(MAX))

# plotting variation in average daily temp over each winter
ggplot(data = HBsum) + 
  geom_point(mapping = aes(x = winter, y = sdave)) +
  labs(x = "Date", y = "Variation in Average Daily Temp")

# linear regression of standard deviation of average daily temp over each winter ~ year
lm1 <- lm(sdmin ~ winter, data=HBsum)
summary(lm1)

# Figure 1
# plot the winter average daily minimum temp by year  
ggplot(data = HBsum, mapping = aes(x = winter, y = minmean)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(x = "Date", y = "Average Minimum Temp (C)") +
  theme_bw()

# linear regression of average minimum daily temp each winter ~ year
lm2 <- lm(minmean ~ winter, data=HBsum)
summary(lm2)

#Figure 2
# plotting variation in minimum daily temp over each winter
ggplot(data = HBsum, mapping = aes(x = winter, y = sdmin)) + 
  geom_point() +
  geom_smooth(method=lm)+
  labs(x = "Date", y = "Variation in Minimum Daily Temp")+
  theme_bw()

# linear regression of standard deviation of minimum daily temp each winter ~ year
lm3 <- lm(sdmin ~ winter, data=HBsum)
summary(lm3)

# plotting variation in max daily temp over each winter
ggplot(data = HBsum, mapping = aes(x = winter, y = sdmax)) + 
  geom_point() +
  geom_smooth(method=lm)

# linear regression of standard deviation of maximum daily temp each winter ~ year
lm4 <- lm(sdmax ~ winter, data=HBsum)
summary(lm4)

# Figure 3
# plot the number of days below freezing (0 degrees Celsius) each winter
ggplot(data = cold, mapping = aes(x = winter, y=below)) + 
  geom_point()+
  geom_smooth(method=lm)+
  labs(x = "Date", y = "# of Below Freezing Days") +
  theme_bw()

# linear regression of number of days below freezing each winter ~ year
lm5 <- lm(below ~ winter, data=cold)
summary(lm5)

