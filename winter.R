# Sage Wentzell-Brehme
# January 31, 2021

# Looking at trends in winter climate data at Hubbard Brook Experimental Forest

# load libraries
library(tidyverse)
library(lubridate)
library(lme4)


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

#subset the data to the HQ site 
HB_snowHQ <- HB_snow %>%
  filter(Site == "STAHQ")

#subset to recent years
HB_snow93 <- HB_snow %>%
  filter(Date >= as.Date("1993-01-01"))

#plot snow depth by frost depth
ggplot(data = HB_snow)+
  geom_point(mapping=aes(x=snow_depth, y=frost_depth))

#summarize climate variables for each winter period
summary <- HB_snowHQ %>%
  dplyr::group_by(WINTER)%>%
  summarize(medfrostdepth = median(frost_depth, na.rm=TRUE),
            sumfrostdepth = sum(frost_depth>=100, na.rm=TRUE),
            sumsnowdepth = sum(snow_depth>=200, na.rm=TRUE), 
            medsnowdepth = median(snow_depth, na.rm=TRUE),
            meansnowdepth = mean(snow_depth, na.rm=TRUE),
            sdsnowdepth = sd(snow_depth, na.rm=TRUE))

summary(HB_snowHQ)
summary(HB_snow93)

#plot the median frost depth by winter
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=medfrostdepth))

#plot the number of days with frost depth >= 100 mm
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=sumfrostdepth))

#plot the median snow depth by winter
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=medsnowdepth))

#plot the sd snow depth by winter
ggplot(data = summary)+
  geom_point(mapping=aes(x=WINTER, y=sdsnowdepth))

#plot distribution of snow depth by winter
boxplot(snow_depth ~ WINTER, data = HB_snowHQ)




ggplot(data = HB_snow93)+
  geom_point(mapping=aes(x=WINTER, y=sumfrostdepth))

# linear regression of 
lm_1 <- lm(meansnowdepth ~ WINTER, data=summary)
summary(lm_1)

table(HB_snow$Date)

test <- HB_snow %>%
  group_by(WINTER)%>%
  summarize(winterstart = min(Date, na.rm=TRUE),
            winterend = max(Date, na.rm=TRUE))


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

