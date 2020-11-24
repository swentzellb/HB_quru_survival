# Sage Wentzell-Brehme
# June 15, 2020
# Summer Science Research Program
# Initial Hubbard Brook historical Red Oak seedling survival and health data (2011-2019)
# UPDATED Sept 8, 2020
# Thesis data analysis - QURU sdlg survival at Hubbard Brook valley wide plots 

#UPDATED Sept 30, 2020 to include 2020 data

#load packages
library(tidyverse)
library(tidyr)
library(ggplot2)
library(raster)
library(lme4)


#set working directory
setwd("~/Thesis/HB_quru_survival/data")

#read in csv files
Seed <- read_csv("HB_VW_Oak_Trans_Survival_2020_CLEAN.csv") #updated sdlg spreadsheet w/2020 data
Seed$VW <- paste0("p", Seed$VW) # add p to beginning of VW column
PlotCoord <- read.csv("../Plot_Table_Basic.csv")

summary(Seed) #summarize the data
table(is.na(Seed$Stat19)) # determine how many NAs are present

#########################################################################
# Data Cleaning
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

#################################################################
#################################################################
# Calculate Survival between 2018 and 2019
#################################################################

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

# replace PD and NF seedling status to 0
Seed[Seed=="PD"] <- 0
Seed[Seed=="NF"] <- 0
Seed[Seed=="nf"] <- 0
Seed$Stat20 <- as.numeric(Seed$Stat20) # fix Stat20 class
Seed$Stat16 <- as.numeric(Seed$Stat16) # fix Stat16 class

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

#Statistical Analysis
#generalized linear models
seedInterval$interval <- as.integer(seedInterval$interval)

#Create a model for survival by Leaf number, Damage, Live/Dead Branch ratio 
#with year interval as a random effect
#interval can't be a character or a factor (but does it mean anything when I make it an integer?)
Mod_1 <- glm(survival ~ leafNumber + seedDamage + brchLvD + (1|interval), data=seedInterval, 
             family = quasibinomial(logit))
summary(Mod_1)

#Model for survival by Leaf Number
Mod_2 <- glm(survival ~ leafNumber, data=seedInterval, family=quasibinomial(logit))
summary(Mod_2)

#Model for survival by Leaf Damage

#Model for survival by Live/Dead Branch Ratio





Mod1 <- glm(Surv_18_19 ~ HQdist, data = Seed_Dist, family=quasibinomial(logit))
summary(Mod1)


Mod2 <- lm(Abund19 ~ HQdist, data = Seed_Dist)
summary(Mod2) #significant p-value

Mod3 <- glm(Surv_18_19 ~ MeanHgt18 + MeanDamage18 + MeanLvs18, data=Seed_Bin,
            family=quasibinomial(logit))
summary(Mod3)

#correlation tests

cor.test(Seed_Bin$MeanDamage18, Seed_Bin$MeanLvs18)

#comparing distance to abundance
cor.test(Seed_Dist$HQdistances, Seed_Dist$Abund19) #p-value <0.01

#comparing distance to sdlg survival
cor.test(Seed_Dist$HQdist, Seed_Dist$Surv_18_19) #p-value 0.054





 
# #change PD and NAs in year columns to 0 
# Seed$Stat20[Seed$Stat20=="PD"| is.na(Seed$Stat20)] <- 0
# Seed$Stat19[Seed$Stat19=="PD"| is.na(Seed$Stat19)] <- 0
# Seed$Stat18[Seed$Stat18=="PD" | is.na(Seed$Stat18)] <- 0 
# Seed$Stat17[Seed$Stat17=="PD" | is.na(Seed$Stat17)] <- 0 
# Seed$Stat16[Seed$Stat16=="PD" |Seed$Stat16=="NF" |is.na(Seed$Stat16)] <- 0 
# Seed$Stat15[Seed$Stat15=="PD" |Seed$Stat15=="NF" |is.na(Seed$Stat15)] <- 0 

# # change PD and NAs in 2014 column to 0
# Seed$Stat14test[Seed$Stat14=="PD" |Seed$Stat14=="NF" |is.na(Seed$Stat14)] <- 0 
# 
# # change PD and NAs in 2013 column to 0
# Seed$Stat13test[Seed$Stat13=="PD" |Seed$Stat13=="NF" |is.na(Seed$Stat13)] <- 0 
# 
# # change PD and NAs in 2012 column to 0
# Seed$Stat12test[Seed$Stat12=="PD" |Seed$Stat12=="NF" |is.na(Seed$Stat13)] <- 0 

# #check all NAs have been removed
# table(is.na(Seed$Stat17test))
# table(Seed$Stat16test)
# table(is.na(Seed$Stat18test))
# table(Seed$Stat13test)

# # create test variables equal to seedling status variables for each yr
# Seed$Stat20test <- Seed$Stat20
# Seed$Stat19test <- Seed$Stat19
# Seed$Stat18test <- Seed$Stat18
# Seed$Stat17test <- Seed$Stat17
# Seed$Stat16test <- Seed$Stat16
# Seed$Stat15test <- Seed$Stat15
# Seed$Stat14test <- Seed$Stat14
# Seed$Stat13test <- Seed$Stat13
# Seed$Stat12test <- Seed$Stat12
# Seed$Stat11test <- Seed$Stat11
# 
# #change PD and NAs in 2020 column to 0 
# Seed$Stat20test[Seed$Stat20=="PD"| is.na(Seed$Stat20)] <- 0
# 
# # change PD and NAs in 2019 column to 0
# Seed$Stat19test[Seed$Stat19=="PD"| is.na(Seed$Stat19)] <- 0
# 
# # change PD and NAs in 2018 column to 0
# Seed$Stat18test[Seed$Stat18=="PD" | is.na(Seed$Stat18)] <- 0 
# 
# # change PD and NAs in 2017 column to 0
# Seed$Stat17test[Seed$Stat17=="PD" | is.na(Seed$Stat17)] <- 0 
# 
# # change PD and NAs in 2016 column to 0
# Seed$Stat16test[Seed$Stat16=="PD" |Seed$Stat16=="NF" |is.na(Seed$Stat16)] <- 0 
# 
# # change PD and NAs in 2015 column to 0
# Seed$Stat15test[Seed$Stat15=="PD" |Seed$Stat15=="NF" |is.na(Seed$Stat15)] <- 0 
# 
# # change PD and NAs in 2014 column to 0
# Seed$Stat14test[Seed$Stat14=="PD" |Seed$Stat14=="NF" |is.na(Seed$Stat14)] <- 0 
# 
# # change PD and NAs in 2013 column to 0
# Seed$Stat13test[Seed$Stat13=="PD" |Seed$Stat13=="NF" |is.na(Seed$Stat13)] <- 0 
# 
# # change PD and NAs in 2012 column to 0
# Seed$Stat12test[Seed$Stat12=="PD" |Seed$Stat12=="NF" |is.na(Seed$Stat13)] <- 0 
# 
# 
# #check all NAs have been removed
# table(is.na(Seed$Stat17test))
# table(Seed$Stat16test)
# table(is.na(Seed$Stat18test))
# table(Seed$Stat13test)
# 
# #change to numeric variables
# Seed$Stat20test <- as.numeric(as.character(Seed$Stat20test))
# Seed$Stat19test <- as.numeric(as.character(Seed$Stat19test))
# Seed$Stat18test <- as.numeric(as.character(Seed$Stat18test))
# Seed$Stat17test <- as.numeric(as.character(Seed$Stat17test))
# Seed$Stat16test <- as.numeric(as.character(Seed$Stat16test))
# Seed$Stat15test <- as.numeric(as.character(Seed$Stat15test))
# Seed$Stat14test <- as.numeric(as.character(Seed$Stat14test))
# Seed$Stat13test <- as.numeric(as.character(Seed$Stat13test))
# Seed$Stat12test <- as.numeric(as.character(Seed$Stat12test))
# 
# Seed$Lvs20 <- as.numeric(as.character(Seed$Lvs20))
# Seed$Lvs19 <- as.numeric(as.character(Seed$Lvs19))
# Seed$Lvs18 <- as.numeric(as.character(Seed$Lvs18))
# Seed$Lvs17 <- as.numeric(as.character(Seed$Lvs17))
# Seed$Lvs16 <- as.numeric(as.character(Seed$Lvs16))
# Seed$Lvs15 <- as.numeric(as.character(Seed$Lvs15))
# Seed$Lvs14 <- as.numeric(as.character(Seed$Lvs14))
# Seed$Lvs13 <- as.numeric(as.character(Seed$Lvs13))
# Seed$Lvs12 <- as.numeric(as.character(Seed$Lvs12))
# 
# #create transition survival variables
# Seed$Surv20 <- Seed$Stat20test
# Seed$Surv19 <- Seed$Stat19test
# Seed$Surv18 <- Seed$Stat18test
# Seed$Surv17 <- Seed$Stat17test
# Seed$Surv16 <- Seed$Stat16test
# Seed$Surv15 <- Seed$Stat15test
# Seed$Surv14 <- Seed$Stat14test
# Seed$Surv13 <- Seed$Stat13test
# Seed$Surv12 <- Seed$Stat12test
# 
# #remove new seedlings in the second year
# Seed$Surv20[Seed$Stat19test==0] <- 0
# Seed$Surv19[Seed$Stat18test==0] <- 0
# Seed$Surv18[Seed$Stat17test==0] <- 0
# Seed$Surv17[Seed$Stat16test==0] <- 0
# Seed$Surv16[Seed$Stat15test==0] <- 0
# Seed$Surv15[Seed$Stat14test==0] <- 0
# Seed$Surv14[Seed$Stat13test==0] <- 0
# Seed$Surv13[Seed$Stat12test==0] <- 0
# 
# 
# #Seed$Surv12[Seed$Stat11test==0] <- 0
# 
# #check new variables
# table(Seed$Surv17)
# table(Seed$Stat17test)
# table(Seed$Stat16test)
# table(Seed$Stat20)

# #Create a smaller dataset
# Seed2 <- Seed %>%
#   mutate(Sdlg = paste(VW, SdlgNum, sep="_")) %>%
#   dplyr::select(Sdlg, Surv20, Surv19, Surv18, Surv17, Surv16, Surv15,
#          Surv14, Surv13, Lvs20, Lvs19, Lvs18, Lvs17, Lvs16, Lvs15, 
#          Lvs14, Lvs13, Lvs12, Damage20, Damage19, Damage18, Damage17,
#          Damage16, Damage15)
# 
# 
# #transform the table with pivot_longer
# Seed2 %>%
#   pivot_longer(Surv20:Damage15,
#                names_to = "interval",
#                names_pattern = "\*(\d{2})",
#                values_to = c("Survival", "Leaves", "Damage")
#   ) # returns an error: '\*' is an unrecognized escape in character string starting ""\*"
# 
# #this sort of works? it returns a table at least 
# #something is definitely still off though
# # copied from the tidyr example
# Seed2 %>% 
#  pivot_longer(Surv20:Damage15,
#                 names_to = c(".value", "set"),
#                 names_pattern = "(.)(.)",
#                 values_to = c("Survival", "Leaves", "Damage")
#    )
# 












#Calculate survival between 2018 and 2019 across all VW plots
Surv_18_19 <- Seed %>%
  filter(Stat18test==1) %>% #filter to only those seedlings alive in 2018 - this removes new seedlings from 2019
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE) #sum all live seedlings in 2019
            /sum(Stat18test==1, na.rm=TRUE), # divide by all live seedlings from 2018
            Count18 = sum(Stat18test==1),
            Count19 = sum(Stat19test==1)) # add column for count of live seedlings at plot in 2019

#Calculate survival between 2017 and 2018 across all VW plots
Surv_17_18 <- Seed %>%
  filter(Stat17test==1) %>% #filter to only those seedlings alive in 2017 - this removes new seedlings from 2018
  summarize(Surv_17_18 = sum(Stat18test==1, na.rm=TRUE) #sum all live seedlings in 2018
            /sum(Stat17test==1, na.rm=TRUE), # divide by all live seedlings from 2017
            Count17 = sum(Stat17test==1),
            Count18 = sum(Stat18test==1)) # add column for count of live seedlings at plot in 2018

#Calculate survival between 2017 and 2019 across all VW plots
Surv_17_19 <- Seed %>%
  filter(Stat17test==1) %>% #filter to only those seedlings alive in 2017 - this removes new seedlings from 2018 and 19
  summarize(Surv_17_19 = sum(Stat19test==1, na.rm=TRUE) #sum all live seedlings in 2019
            /sum(Stat17test==1, na.rm=TRUE), # divide by all live seedlings from 2017
            Count17 = sum(Stat17test==1),
            Count19 = sum(Stat19test==1)) # add column for count of live seedlings at plot in 2019

#Calculate survival between 2017 and 2019 across 10 hemlock plots
Surv_17_19_hem <- Seed %>%
  filter(Stat17test==1)%>% 
  filter(VW == "p354" | VW == "p352" | VW =="p351" | VW =="p337"|VW=="p338"
         |VW=="p317"|VW=="p322"|VW=="p379"|VW=="p378"|VW=="p396")%>%
  summarize(Surv = sum(Stat19test==1, na.rm=TRUE) 
          /sum(Stat17test==1, na.rm=TRUE),
          Count17 = sum(Stat17test==1),
          Count19 = sum(Stat19test==1))


########################################################################################################
########################################################################################################
#Binning
# Combine PlotTag & TransDir into one unique transect ID
Seed <- Seed %>%
  mutate(PlotTrans = paste(VW, TransDir, sep = "_"))

# For each unique transect ID (PlotTrans), bin to 20m
TransBins <- seq(0, 100, by = 20)

# Make column with distance changed to bin number
Seed <- Seed %>%
  mutate(DistBin = findInterval(Distance, TransBins)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))


###############################################################
# summarize seedling survival between 2018 and 2019 to 20m bins
Seed_Bin <- Seed %>%
  filter(Stat18test==1) %>%
  group_by(PlotTrans, DistBin) %>%
  summarize(Surv_18_19 = sum(Stat19test, na.rm=TRUE)/sum(Stat18test, na.rm=TRUE),
            Abund18 = sum(Stat18test), 
            Abund19 = sum(Stat19test), 
            MeanHgt18 = mean(as.numeric(Hgt18), na.rm=TRUE),
            MeanLvs18 = mean(as.numeric(Lvs18), na.rm=TRUE),
            MeanDamage18 = mean(Damage18, na.rm=TRUE)) 

#summarize seedling survival between 2018 and 2019 to plot transect level
Seed_Bin2 <- Seed %>%
  filter(Stat18test==1) %>%
  group_by(PlotTag, PlotTrans, TransDir) %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test), 
            Abund19 = sum(Stat19test),
            MeanDamage19 = mean(Damage19, na.rm=TRUE))

# subset 20m survival data to 10 VW transects with Enviromental data
Seed_Bin_Enviro <- Seed %>%
  filter(Stat18test==1) %>% 
  filter(VW == "p354" | VW == "p352" | VW =="p351" | VW =="p337"|VW=="p338" |VW=="p317"|VW=="p322"|VW=="p379"|VW=="p378"|VW=="p396")
#################################################################################
#################################################################################
#Add in Distance data - converting from UTM 

# fix column name
colnames(PlotCoord)[1] <- "PlotTag"

# define location of road near Pierce lab
USFS <- c(283404.12, 4868995.17)

#Adjust Coordinates by transect side (add 50 for W and -50 for E)
#create dataset with PlotTrans as identifier
test <- Seed_Bin2 %>%
  dplyr::select(TransDir, PlotTag, PlotTrans) 

#add another column equal to the transect direction
test$Adj <- test$TransDir

#change E/W to adding or subtracting 50 meters
test$Adj[test$TransDir=="E"] <- -50
test$Adj[test$TransDir=="W"]<- 50
test$Adj <- as.numeric(test$Adj) # make a numeric variable

#join the Plot coordinate dataset with the adjustments to east distances
CoordAdj <- left_join(PlotCoord, test, by="PlotTag")

#set NAs equal to 0
CoordAdj$Adj[is.na(CoordAdj$Adj)] <- 0

#adjust Easting coordinate based on transect side
CoordAdj <- CoordAdj %>% dplyr::mutate(UTM.EastingAdj = UTM.Easting + Adj)

#make smaller dataset
plot_locs <- cbind(CoordAdj$UTM.EastingAdj, CoordAdj$UTM.Northing)

#find distance between that point and others in the valley
plot_distances <- pointDistance(USFS, plot_locs, lonlat=FALSE)

# attach back to main dataset
CoordAdj$HQdist <- plot_distances

#make the dataset smaller - only have PlotTag (key variable), distances
CoordAdj <- CoordAdj %>%
  dplyr::select(PlotTag, PlotTrans, HQdist)%>%
  group_by(PlotTrans, DistBin) %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test==1), 
            Abund19 = sum(Stat19test)) 

# total seedling survival across subset 10 VW transects (Enviro data)
Seed_Bin_Enviro_test <- Seed %>%
  filter(Stat18test==1) %>% 
  filter(VW == "354" | VW == "352" | VW =="351" | VW =="337"|VW=="338"
         |VW=="317"|VW=="322"|VW=="379"|VW=="378"|VW=="396") %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)
            /sum(Stat18test==1, na.rm=TRUE))

#PLOT 4
#plot of 2019 seedling abundance by plot distance to entrance of valley

#join Seedling Abundance at Plot Transect level with distances to valley entrance
Seed_Dist <- left_join(Seed_Bin2, CoordAdj, by="PlotTrans")

#create plot
ggplot(Seed_Dist, aes(x=HQdist, y=Abund19))+
  geom_point(color = "royalblue3", size=2.5)+
  theme_light()+
  labs(x= "Distance to valley entrance (m)",
       y = "Number of seedlings") +
  theme(axis.title = element_text(size=12))
  theme(axis.ticks = element_text(size=11))
