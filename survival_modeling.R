# Sage Wentzell-Brehme
# January 22, 2021
# Senior Thesis
# Data analysis
# QURU sdlg survival at Hubbard Brook valley wide plots 
# between 2011 and 2020

# Reordering spreadsheet of seedling survival data



#load packages
library(tidyverse)
library(tidyr)
library(ggplot2)
library(lme4)


#set working directory
setwd("~/Thesis/HB_quru_survival/data")

#read in csv files
Seed <- read_csv("HB_VW_Oak_Trans_Survival_2020_CLEAN.csv") #updated sdlg spreadsheet w/2020 data
Seed$VW <- paste0("p", Seed$VW) # add p to beginning of VW column

#########################################################################
# Data Cleaning
#########################################################################

#add PlotTag column to seed dataset to match coordinate dataset
Seed$PlotTag <- Seed$VW

#separate into Point and Plot Direction - N/S columns
#PlotDir = Plot Direction - N/S of transect line
Seed <- Seed %>% separate(Point, c("Point", "PlotDir"), sep="_")

#separate into Distance along transect and Transect Direction - E/W columns
#by taking the last character in the string
#Distance - along transect line
#TransDir - Transect Direction to the East or West
Seed <- Seed %>% separate(Point, c("Distance", "TransDir"), sep = -1)

#separate and remove point type, plot off transect or on center line of transect, from distance along transect
Seed <- Seed %>% separate(Distance, c("PointType", "Distance"), sep = 1)


# make unique column for seedling ID
seedIDs <- mutate(Seed, seedID = paste(VW,SdlgNum,TransDir,Distance,sep="-"))

# test to make sure there aren't duplicates
length(seedIDs$seedID) #1375
length(unique(seedIDs$seedID)) #1369 - there are 6 seedlings in here twice

# find seedling duplicates: 1. count/tag duplicates based on seedID, 
# 2. join duplicate tags back to seedIDs data frame
# 3. filter to duplicate seedlings
find_duplicates <- seedIDs %>% 
  group_by(seedID) %>%
  summarize(count = n(),
            dups = ifelse(count>1,T,F)) %>%
  right_join(seedIDs) %>%
  filter(dups == T)

# We can decide what to do with duplicates later - 
# for now I'll keep rows where Surv20 not equal to 0
# 1. remove_dups gets row numbers (index) where duplicates & Surv 2020 = 0 (keep 1/PD)
# 2. make duplicate rows NA in data frame and filter out NA rows
remove_dups <- which(seedIDs$seedID %in% find_duplicates$seedID & seedIDs$Stat20 == 0)
seedIDs[remove_dups,] <- rep(NA, ncol(seedIDs))
seedIDs <- filter(seedIDs, !is.na(seedID))

# test to make sure there aren't duplicates again
length(seedIDs$seedID) #1369
length(unique(seedIDs$seedID)) #1369 - we're good

# move edited data frame back to Seed name, remove some extraneous columns
Seed <- dplyr::select(seedIDs, -SdlgSpp, -SdlgNum, -PointType, -Distance, 
                      -TransDir, -PlotDir,-Along, -TransEW, -FromTrans, -NSline)

Seed$Stat20 <- as.numeric(Seed$Stat20) # fix Stat20 class
Seed$Stat16 <- as.numeric(Seed$Stat16) # fix Stat16 class
# replace PD and NF seedling status to 0
Seed[Seed=="PD"] <- 0
Seed[Seed=="NF"] <- 0
Seed[Seed=="nf"] <- 0


# Calculate & format the survival interval
# 1. Select only seedID and columns with survival data
seedSurv <- dplyr::select(Seed, seedID, grep("Stat",colnames(Seed)))

# 2. Loop over each year and calculate the survival over the interval
# & 3. Build a data frame where year interval & survival are columns
for(i in 2:(ncol(seedSurv)-1)){
  # Calculate survival over the 2-summer interval, 1 = survived, 0 = died, NA = other
  surv_interval <- ifelse(seedSurv[,i]==1 & seedSurv[,(i+1)]==1, 1, # Survived
                          ifelse(seedSurv[,i]==0 & seedSurv[,(i+1)]==1, 0, # Died
                                 ifelse(seedSurv[,i]==1 & is.na(seedSurv[,(i+1)]), NA, # Born
                                        ifelse(is.na(seedSurv[,i]), NA, NA)))) # Who knows
  
  # Make data frame with seedID, interval code, survival 
  surv_df <- data.frame(seedID = seedSurv[,1], 
                        interval = paste0("20",as.numeric(str_sub(colnames(surv_interval),5,6))-1,"-",
                                          "20",str_sub(colnames(surv_interval),5,6)),
                        survival = as.vector(surv_interval))
  
  # Build a long data frame
  if(i == 2){
    seedSurvInt <- surv_df
  } else {
    seedSurvInt <- bind_rows(seedSurvInt, surv_df)
  }
}

#4. Format leaf number, damage, height, & ratio live branches to interval format
seedLvs <- dplyr::select(Seed, seedID, grep("Lvs",colnames(Seed))) 
seedLvs[,2:10] <- lapply(seedLvs[,2:10], as.numeric)
seedLvs <- pivot_longer(seedLvs, !seedID, names_to = "interval", 
                        values_to = "leafNumber") %>%
  mutate(interval = paste0("20",as.numeric(str_sub(interval,4,5))-1,"-",
                           "20",str_sub(interval,4,5))) %>%
  filter(interval != "2011-2012") # leaf number started in 2012

seedHgt <- dplyr::select(Seed, seedID, grep("Hgt",colnames(Seed))) 
seedHgt[,2:ncol(seedHgt)] <- lapply(seedHgt[,2:ncol(seedHgt)], as.numeric)
seedHgt <- pivot_longer(seedHgt, !seedID, names_to = "interval", 
                        values_to = "seedHeight") %>%
  mutate(interval = paste0("20",as.numeric(str_sub(interval,4,5))-1,"-",
                           "20",str_sub(interval,4,5))) 

seedDmg <- dplyr::select(Seed, seedID, grep("Damage",colnames(Seed))) 
seedDmg[,2:ncol(seedDmg)] <- lapply(seedDmg[,2:ncol(seedDmg)], as.numeric)
seedDmg <- pivot_longer(seedDmg, !seedID, names_to = "interval", 
                        values_to = "seedDamage") %>%
  mutate(interval = paste0("20",as.numeric(str_sub(interval,7,8))-1,"-",
                           "20",str_sub(interval,7,8))) 

seedBrch <- dplyr::select(Seed, seedID, grep("Brch",colnames(Seed))) %>%
  dplyr::select(-Brch15) 
seedBrch[,2:ncol(seedBrch)] <- lapply(seedBrch[,2:ncol(seedBrch)], as.numeric)
seedBrch <- seedBrch %>%
  mutate(BrchLvD20 = LiveBrch20/(LiveBrch20+DeadBrch20),
         BrchLvD19 = LiveBrch19/(LiveBrch19+DeadBrch19),
         BrchLvD18 = LiveBrch18/(LiveBrch18+DeadBrch18)) %>%
  dplyr::select(seedID, BrchLvD20, BrchLvD19, BrchLvD18)
seedBrch <- pivot_longer(seedBrch, !seedID, names_to = "interval", 
                         values_to = "brchLvD") %>%
  mutate(interval = paste0("20",as.numeric(str_sub(interval,8,9))-1,"-",
                           "20",str_sub(interval,8,9))) 

# 5. Combine everything together into a data frame for the stat model
seedInterval <- left_join(seedSurvInt, seedLvs) %>%
  left_join(seedDmg) %>%
  left_join(seedBrch) %>%
  filter(!is.na(survival)) # remove lines where survival is NA


# Find the seedling age in each interval based on birth year
seedsToAge <- unique(seedInterval$seedID)

for(s in 1:length(seedsToAge)){
  
  # Get seedID data from giant df
  seedID_df <- seedInterval[which(seedInterval$seedID == seedsToAge[s]),] %>%
    separate(interval, into = c("start_int","end_int"), sep = "-",
             remove = FALSE)
  
  # Find birth year
  seedID_birthyear <- Seed$BirthYear[which(Seed$seedID == seedsToAge[s])]
  
  # If birth year is NA, assume first year of census is birth year
  # 70 seedlings in data fall into this category
  if(is.na(seedID_birthyear)){
    
    # Assign min/max ages from first census
    min_age <- 1
    max_age <- as.numeric(max(seedID_df$start_int)) - as.numeric(min(seedID_df$start_int)) + 1
    ages <- rev(seq(min_age, max_age))
    
    # Test to make sure seedID measured for each year alive
    if(length(ages) == nrow(seedID_df)){
      seedID_df$yearsAlive <- ages
      
    } else { # Match the ages with the years found in the census
      yearsShouldBeAlive <- seq(as.numeric(min(seedID_df$start_int)), 
                                as.numeric(max(seedID_df$start_int)))
      yearsAlive_df <- data.frame(yearsShouldBeAlive,
                                  ages = ages)
      seedID_df$yearsAlive<- yearsAlive_df$age[yearsAlive_df$yearsShouldBeAlive %in% seedID_df$start_int]
      #print(paste("s =",s,"seedling not in census every year alive"))
      #print(seedID_df)
    }
  } else { # Else use the birth year to assign age, n = 1166 seedlings
    
    # Find min/max ages from birth year
    min_age <- as.numeric(min(seedID_df$start_int)) - seedID_birthyear + 1
    max_age <- as.numeric(max(seedID_df$start_int)) - seedID_birthyear + 1
    ages <- rev(seq(min_age, max_age))
    
    # Test to make sure seedID was recorded in each year
    if(length(ages) == nrow(seedID_df)){
      seedID_df$yearsAlive <- ages
    } else {
      # Match the ages with the years found in the census
      yearsShouldBeAlive <- seq(as.numeric(min(seedID_df$start_int)), 
                                as.numeric(max(seedID_df$start_int)))
      yearsAlive_df <- data.frame(yearsShouldBeAlive,
                                  ages = ages)
      seedID_df$yearsAlive<- yearsAlive_df$age[yearsAlive_df$yearsShouldBeAlive %in% seedID_df$start_int]
      #print(paste("s =",s,"seedling not in census every year alive"))
      #print(seedID_df)
    }
  }
  # Compile data frame with seed ages
  if(s == 1) {
    seedIntervalAge <- seedID_df
  } else {
    seedIntervalAge <- bind_rows(seedIntervalAge, seedID_df)
  }
}

# Save data frame with cleaned and aggregated data
write_csv(seedIntervalAge, "seedIntervalAge_tidy.csv")


#########################################################
#########################################################
#Statistical Analysis
#generalized linear models

#change variable classes to factor
seedInterval$seedDamage <- as.factor(seedInterval$seedDamage)
seedInterval$interval <- as.factor(seedInterval$interval)

#Create a model for survival by Leaf number, Damage, Live/Dead Branch ratio 
Mod_1 <- glm(survival ~ leafNumber + seedDamage + brchLvD + yearsAlive, data=seedIntervalAge, 
             family = quasibinomial(logit))
summary(Mod_1)

# Find model RMSE
RSS <- c(crossprod(Mod_1$residuals))
MSE <- RSS / length(Mod_1$residuals)
RMSE <- sqrt(MSE)

#Create a model for survival by Leaf number, Damage, Live/Dead Branch ratio 
Mod_1_ <- glm(survival ~ leafNumber + seedDamage + brchLvD, data=seedInterval, 
             family = quasibinomial(logit))
summary(Mod_1_)

#Model for survival by Leaf Number
Mod_2 <- glm(survival ~ leafNumber, data=seedInterval, family=quasibinomial(logit))
summary(Mod_2)

#Model for survival by Leaf Damage
Mod_3 <- glm(survival ~ seedDamage, data=seedInterval, family=quasibinomial(logit))
summary(Mod_3)
TukeyHSD(aov(Mod_3)) #compare categories

#Model for survival by Live/Dead Branch Ratio
Mod_4 <- glm(survival ~ brchLvD, data=seedInterval, family=quasibinomial(logit))
summary(Mod_4)

#################################################################
#Add in year as a random effect

#Model for survival by Leaf number with Year interval as a random effect
M1 <- glmer(survival ~ leafNumber + yearsAlive + (1|interval), 
            data=seedIntervalAge, family="binomial")
summary(M1)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M1))) #residual sum of squares
MSE <- RSS / length(residuals(M1)) #mean squared error
RMSE <- sqrt(MSE)

#Model for survival by Leaf number with Year interval as a random effect
M1_ <- glmer(survival ~ leafNumber + (1|interval), data=seedInterval, family="binomial")
summary(M1_)
ranef(M1_)
#AIC = 478.6

#Model for survival by Leaf number and Live/Dead Branch ratio with Year interval as a random effect
M2 <- glmer(survival ~ leafNumber + brchLvD + (1|interval), data=seedInterval, 
            family = "binomial")
summary(M2)
ranef(M2)
#AIC = 61.8


#Model for survival by leaf number
M3 <- glm(survival ~ leafNumber, data = seedInterval, family ="binomial")
summary(M3)
#AIC = 513

#Model for survival by branch ratio with Year interval as a random effect
M4 <- glmer(survival ~ brchLvD + (1|interval), data = seedInterval, family ="binomial")
summary(M4)
ranef(M4)
#AIC = 163

#Model for survival by Leaf Damage
M5 <- glm(survival ~ seedDamage, data=seedInterval, family="binomial")
summary(M5)
stepAIC(M5)
TukeyHSD(aov(M5)) # compare categories

M6 <- glmer(survival ~ seedDamage + (1|interval), data=seedInterval, family="binomial")
summary(M6)
###################################################
##Correlation tests

#correlation between leaf number and branch ratio
cor.test(seedInterval$leafNumber, seedInterval$brchLvD)
#correlated p-value <0.01 and cor = 0.17

#correlation between leaf number and survival
cor.test(seedInterval$leafNumber, seedInterval$survival)
#correlated p-value <0.01 and cor = 0.49

#correlation between branch ratio and survival
cor.test(seedInterval$brchLvD, seedInterval$survival)
#correlated p-value <0.01 and cor = 0.1


