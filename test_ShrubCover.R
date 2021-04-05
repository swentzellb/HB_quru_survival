# Sage Wentzell-Brehme
# September 14, 2020
# Thesis
# Combining Hubbard Brook historical Red Oak seedling survival and health data (2011-2020)
# with environmental data from 2019 (and 2020 shrub cover and substrate data)

# Combining: Hubbard Brook historical QURU seedling survival and health data (2011-2019)
# and Hubbard Brook environmental data from 2019 on 10 subset transects in hemlock-dom areas

#Figures

#load packages
library(tidyverse)
library(tidyr)
library(lme4)
library(MASS)
library(raster)
library(ggplot2)

#set working directory
setwd("~/Thesis/HB_quru_survival/data")

#read in csv files
Seed <- read.csv("HB_VW_Oak_Trans_Survival_2020_CLEAN.csv", header=TRUE,  fileEncoding = "UTF-8-BOM") # file encoding fixes the strange symbol added in column title
Seed$VW <- paste0("p", Seed$VW) # add p to beginning of VW column (same as PlotTag)
# historical survival data

Cover_Sub <- read.csv("HB_VW_Oak_Trans_ShrubCover_2020_CLEAN.csv", header=TRUE,  fileEncoding = "UTF-8-BOM") 
#refers to canopy and shrub cover and substrate
Seed <- read_csv("HB_VW_Oak_Trans_Survival_2020_CLEAN.csv") #updated sdlg spreadsheet w/2020 data

#####################################################################################################################
#####################################################################################################################
###Data Cleaning - Cover Substrate

#separate out and create new identifying columns
Cover_Sub<- Cover_Sub %>% 
  mutate(Point2 = Point) %>%
  separate(Point2, c("Point2", "PlotDir"), sep="_") %>%
  separate(Point2, c("Distance", "TransDir"), sep = -1) %>%
  separate(Distance, c("PointType", "Distance"), sep = 1) %>%
  mutate(PlotPoint = paste(PlotTag, Point, sep = "_")) #unique point ID

##########################################################
##Binning
# For each unique transect ID (PlotTrans), bin to 20m
TransBins <- seq(0, 100, by = 20)

# create new bins and add identifying columns
Cover_Sub <- Cover_Sub %>%
  mutate(DistBin = findInterval(Distance, TransBins)) %>% # bin to 20m
  mutate(PlotTrans = paste(PlotTag, TransDir, sep = "_")) %>% 
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_")) #unique bin ID
  
#remove p and change PlotTag name to VW
 Cover_Sub <- Cover_Sub %>%
   separate(PlotTag, c("PlotTag", "VW"), sep=1)


#########################################################################
#########################################################################
#Data Cleaning for Historical seedling dataset
 
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
 
 # Combine PlotTag & TransDir into one unique transect ID
 Seed <- Seed %>%
   mutate(PlotTrans = paste(VW, TransDir, sep="_")) %>%
   mutate(DistBin = findInterval(Distance, TransBins)) %>%
   mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_")) # add a unique ID for each bin
 
 #################################################################
 #################################################################
 # Calculate Survival across intervals
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
 d <- rep(NA, ncol(seedIDs))
 d_df <- do.call("rbind", replicate(6, d, simplify = FALSE))
 seedIDs[remove_dups,] <- d_df
 seedIDs <- filter(seedIDs, !is.na(seedID))
 
 # test to make sure there aren't duplicates again
 length(seedIDs$seedID) #1369
 length(unique(seedIDs$seedID)) #1369 - we're good
 
 # move edited data frame back to Seed name, remove some extraneous columns
 # Seed <- dplyr::select(seedIDs, -SdlgSpp, -SdlgNum, -PointType, -Distance, 
                       # -TransDir, -PlotDir,-Along, -TransEW, -FromTrans, -NSline)
 
 Seed <- dplyr::select(seedIDs, -SdlgSpp, -SdlgNum, -PointType, -Distance, 
                        -TransDir, -PlotDir,-Along, -TransEW, -FromTrans, -NSline)
 
 # replace PD and NF seedling status to 0
 Seed[Seed=="PD"] <- "0"
 Seed[Seed=="NF"] <- "0"
 Seed[Seed=="nf"] <- "0"
 Seed$Stat20 <- as.numeric(Seed$Stat20) # fix Stat20 class
 Seed$Stat16 <- as.numeric(Seed$Stat16) # fix Stat16 class
 
 #######################################################################################################
 #######################################################################################################
 # Calculate Survival between 2018 and 2019
 #######################################################################################################
 
 # create test variables equal to seedling status variables for each yr
 Seed$Stat19test <- Seed$Stat19
 Seed$Stat18test <- Seed$Stat18
 
 # change PD and NAs in 2019 column to 0
 Seed$Stat19test[Seed$Stat19=="PD"| is.na(Seed$Stat19test)] <- 0
 
 # change PD and NAs in 2018 column to 0
 Seed$Stat18test[is.na(Seed$Stat18test) | Seed$Stat18=="PD"] <- 0 
 
 #check all NAs have been removed
 table(is.na(Seed$Stat18test))
 
 # change class to numeric
 Seed$Stat18test <- as.numeric(Seed$Stat18test)
 Seed$Stat19test <- as.numeric(Seed$Stat19test)
 
 #check correct arrangement - should be 852=0
 table(Seed$Stat19test)
 
 ##################################################################
 ## # same for survival between 2019 and 2020
 ##################################################################
 
 # create test variables equal to seedling status variables for each yr
 Seed$Stat20test <- Seed$Stat20
 
 # change PD and NAs in 2020 column to 0
 Seed$Stat20test[Seed$Stat20=="PD"| Seed$Stat20=="NF"|Seed$Stat20=="nf"|is.na(Seed$Stat20test)] <- 0
 
 #check all NAs have been removed
 table(is.na(Seed$Stat20test))
 
 # change class to numeric
 Seed$Stat20test <- as.numeric(Seed$Stat20test)
 
 #check correct arrangement - should be 
 table(Seed$Stat20test)
 
 #############################################################################
 # individual seedling survival status from 2018 to 2019
 #############################################################################
 
 # add another column to Seed dataset equal to stat19
 Seed$Stat18_19 <- Seed$Stat19test
 
 # removes seedlings that were not present in 2018 (set to NA)
 Seed$Stat18_19[Seed$Stat18test==0]<- NA
 
 table(is.na(Seed$Stat18_19)) #should be 654 false
 
 # add another column to Seed dataset equal to stat19
 Seed$Stat19_20 <- Seed$Stat20test
 
 # removes seedlings that were not present in 2018 (set to NA)
 Seed$Stat19_20[Seed$Stat19test==0]<- NA
 
 table(is.na(Seed$Stat19_20)) 
 
 
 
 Seed_Bin <- Seed %>%
   # filter(Stat18test==1) %>%
   group_by(PlotTrans, DistBin) %>%
   summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
             Abund18 = sum(Stat18test==1), 
             Abund19 = sum(Stat19test, na.rm=TRUE)) %>%
   mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))
 
 
 
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
 
 
 # 5. Combine everything together into a data frame for the stat model
 seedInterval <- seedSurvInt %>%
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
 

 
 # summarize seedling survival between 2018 and 2019 to 20m bins
 Seed_Bin <- Seed %>%
 #  filter(Stat18test==1) %>%
   group_by(PlotTrans, DistBin) %>%
   summarize(propSurv = sum(survival==1)/n()) %>%
   mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))
 
 # using seedInterval data frame which includes seedling data, 
 # indiv characteristics
 # & survival of seedlings across year intervals
 
 Enviro <- c("2018-2019", "2019-2020")

 
 seedSurvival <- seedIntervalAge %>%
   filter(interval %in% Enviro) %>%
   group_by(seedID)
   summarise(propSurv = sum(survival==1)/n())
 
 test <- seedIntervalAge %>%
   filter(interval=="2019-2020")
 
 ggplot(seedSurvival, aes(x=interval, y=propSurv))+
   geom_bar(stat="identity", fill="skyblue4")+
   ylim(0,1)+
   xlab("Survival Year Interval")+
   ylab("Proportion Seedlings Survived")+
   theme_bw()+
   theme(axis.text.x = element_text(angle=19))
 
 
 
 
 
 
 

#separate into Point and Plot Direction (N/S of transect) columns
#PlotDir = Plot Direction - N/S of transect line
Seed <- Seed %>% separate(Point, c("Point", "PlotDir"), sep="_") %>% 
  separate(Point, c("Distance", "TransDir"), sep = -1) %>% #Distance along transect and Transect Direction (E/W)
  separate(Distance, c("PointType", "Distance"), sep = 1) #remove point type (C/P for center line or point to side)

#Define variables as numeric for analysis
Seed$Lvs18 <- as.numeric(Seed$Lvs18)
Seed$Hgt18 <- as.numeric(Seed$Hgt18)

##########################################################################################
#Binning

# Combine PlotTag & TransDir into one unique transect ID
Seed <- Seed %>%
  mutate(PlotTrans = paste(VW, TransDir, sep="_")) %>%
  mutate(DistBin = findInterval(Distance, TransBins)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_")) # add a unique ID for each bin

####################################################################
##New Variables

# add another column equal to PlotTag - PlotTag is primary key to join with another dataset
Seed$PlotTag <- Seed$VW

#make age variable for seedling age in 2018
Seed$Age18 <- (2018 - Seed$BirthYear) + 1

#######################################################################################################
#######################################################################################################
# Calculate Survival between 2018 and 2019
#######################################################################################################

# create test variables equal to seedling status variables for each yr
Seed$Stat19test <- Seed$Stat19
Seed$Stat18test <- Seed$Stat18

# change PD and NAs in 2019 column to 0
Seed$Stat19test[Seed$Stat19=="PD"| is.na(Seed$Stat19test)] <- 0

# change PD and NAs in 2018 column to 0
Seed$Stat18test[is.na(Seed$Stat18test) | Seed$Stat18=="PD"] <- 0 

#check all NAs have been removed
table(is.na(Seed$Stat18test))

# change class to numeric
Seed$Stat18test <- as.numeric(Seed$Stat18test)
Seed$Stat19test <- as.numeric(Seed$Stat19test)

#check correct arrangement - should be 852=0
table(Seed$Stat19test)

##################################################################
## # same for survival between 2019 and 2020
##################################################################

# create test variables equal to seedling status variables for each yr
Seed$Stat20test <- Seed$Stat20

# change PD and NAs in 2020 column to 0
Seed$Stat20test[Seed$Stat20=="PD"| Seed$Stat20=="NF"|Seed$Stat20=="nf"|is.na(Seed$Stat20test)] <- 0

#check all NAs have been removed
table(is.na(Seed$Stat20test))

# change class to numeric
Seed$Stat20test <- as.numeric(Seed$Stat20test)

#check correct arrangement - should be 
table(Seed$Stat20test)

#############################################################################
# individual seedling survival status from 2018 to 2019
#############################################################################

# add another column to Seed dataset equal to stat19
Seed$Stat18_19 <- Seed$Stat19test

# removes seedlings that were not present in 2018 (set to NA)
Seed$Stat18_19[Seed$Stat18test==0]<- NA

table(is.na(Seed$Stat18_19)) #should be 654 false

# add another column to Seed dataset equal to stat19
Seed$Stat19_20 <- Seed$Stat20test

# removes seedlings that were not present in 2018 (set to NA)
Seed$Stat19_20[Seed$Stat19test==0]<- NA

table(is.na(Seed$Stat19_20)) 

##############################
################################
## Find total survival between 2019 and 2020
Surv_19_20 <- Seed %>%
  filter(Stat19test==1) %>% #filter to only those seedlings alive in 2019 - this removes new seedlings from 2020
  summarize(Surv_18_19 = sum(Stat20test==1, na.rm=TRUE) #sum all live seedlings in 2020
            /sum(Stat19test==1, na.rm=TRUE), # divide by all live seedlings from 2019
            Count19 = sum(Stat19test==1),
            Count20 = sum(Stat20test==1)) # add column for count of live seedlings at plot in 2020


Surv_18_19 <- Seed %>%
  filter(Stat18test==1) %>% #filter to only those seedlings alive in 2018 - this removes new seedlings from 2019
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE) #sum all live seedlings in 2019
            /sum(Stat18test==1, na.rm=TRUE), # divide by all live seedlings from 2018
            Count18 = sum(Stat18test==1),
            Count19 = sum(Stat19test==1)) # add column for count of live seedlings at plot in 2019


################################################################################################################
#############################################################################################################
#### Binning
############################################################################################################
# summarize seedling survival between 2018 and 2019 to 20m bins
Seed_Bin <- Seed %>%
  filter(Stat18test==1) %>%
  group_by(PlotTrans, DistBin) %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test==1), 
            Abund19 = sum(Stat19test, na.rm=TRUE)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))

#summarize seedling survival between 2018 and 2019 to plot transect level
Seed_Bin2 <- Seed %>%
  filter(Stat18test==1) %>%
  group_by(PlotTag, PlotTrans, TransDir) %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test), 
            Abund19 = sum(Stat19test),
            MeanDamage19 = mean(Damage19, na.rm=TRUE))

#summarize seedling survival between 2018 and 2019 to transect level
Seed_Bin3 <- Seed %>%
  filter(Stat18test==1) %>%
  group_by(PlotTag) %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test), 
            Abund19 = sum(Stat19test),
            MeanDamage18 = mean(Damage18, na.rm=TRUE),
            MeanLvs18 = mean(Lvs18, na.rm=TRUE), 
            MeanAge = mean(Age18, na.rm=TRUE))

# summarize 2018 to 2019 seedling survival to 10 VW transects with Enviromental data
Seed_Bin_Enviro <- Seed %>%
  dplyr::filter(Stat18test==1) %>% 
  dplyr::filter(VW == "p354" | VW == "p352" | VW =="p351" | VW =="p337"|VW=="p338"
         |VW=="p317"|VW=="p322"|VW=="p379"|VW=="p378"|VW=="p396") %>%
  group_by(PlotTrans, DistBin) %>%
  summarize(Surv_18_19 = sum(Stat19test, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test), 
            Abund19 = sum(Stat19test)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))%>%
  dplyr::select(PlotTransBin, Abund18, Abund19, Surv_18_19, PlotTrans)

#summarize Canopy Cover densio, Shrub Cover and Available Substrate to 20m bins
Shrub_Bin <- Cover_Sub %>%
  group_by(PlotTag, PlotTrans, DistBin)%>%
  summarize(ShrubCover_mean = mean(ShrubCover, na.rm=TRUE))%>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))

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
  dplyr::select(PlotTag, PlotTrans, HQdist)

#################################################################
#######################################################################################
###Joining
#join Cover subset data with Seedling Survival dataset (binned at 20m)
Cover_Surv <- left_join(Shrub_Bin, Seed_Bin_Enviro, by="PlotTransBin")

#join Seedling Abundance at Plot Transect level with distances to valley entrance
Seed_Dist <- left_join(Seed_Bin2, CoordAdj, by="PlotTrans")

##############################################################################
#join new canopy_cover_substrate datasets (binned to 20m level) with survival

#subset cover and survival dataset
Cover_Surv2 <- Cover_Surv %>%
  dplyr::select(ShrubCover_mean, PlotTransBin, Surv_18_19, Sub_mean, Abund18)

# create new dataset with canopy cover, other environmental measures, and survival data
Can_Cover_Surv <- left_join(Cover_Surv2, Can_2, by = "PlotTransBin")


#################################################################
##################################################################
#correlation tests
#comparing Canopy Cover and Shrub Cover
cor.test(Cover_Surv$CanCover_mean, Cover_Surv$ShrubCover_mean) #p-value <0.001

#comparing Canopy Cover and Available substrate
cor.test(Cover_Surv$CanCover_mean, Cover_Surv$Sub_mean) #p-value 0.2

#comparing canopy cover to survival 
cor.test(Cover_Surv$CanCover_mean, Cover_Surv$Surv_18_19) #p-value 0.2

#comparing shrub cover to survival
cor.test(Cover_Surv$ShrubCover_mean, Cover_Surv$Surv_18_19) #p-value 0.1

#comparing available substrate to abundance
cor.test(Cover_Surv$Sub_mean, Cover_Surv$Abund19) #p-value 0.49

#comparing canopy cover photos and shrub cover
cor.test(Can_Cover_Surv$ShrubCover_mean, Can_Cover_Surv$PerCnpyOpenTotal) #p-value 0.6

# number of leaves and age in 2018
cor.test(Seed$Lvs18, Seed$Age18) #correlated

#canopy photos and seedling abundance in 2018
cor.test(Can_Cover_Surv$PerCnpyOpenTotal, Can_Cover_Surv$Abund18)

#canopy photos and shrub cover 
cor.test(Can_Cover_Surv$PerCnpyOpenTotal, Can_Cover_Surv$ShrubCover_mean)

# height in 2018 and age in 2018 
cor.test(Seed$Hgt18, Seed$Age18) # very correlated


#######################################################################################
#######################################################################################
#Statistical Modeling
#######################################################################################
#modeling plot level (20m bin) sdlg survival

#model w/ densiometer canopy, shrub cover, and substrate 
Mod1 <- glm(Surv_18_19 ~ CanCover_mean + ShrubCover_mean 
            + Sub_mean, data = Cover_Surv, family=quasibinomial(logit))#generalized linear model with quasibinomial family
summary(Mod1)# show model output

#model w/ densiometer canopy and substrate 
Mod2 <- glm(Surv_18_19 ~ CanCover_mean + Sub_mean, data = Cover_Surv, family=quasibinomial(logit))
summary(Mod2)

#model w/ canopy photos and substrate cover
Mod3 <- glm(Surv_18_19 ~ PerCnpyOpenTotal + Sub_mean, data = Can_Cover_Surv,
            family=quasibinomial(logit))
summary(Mod3)

#model w/canopy photos, available substrate, and shrub cover
Mod4 <- glm(Surv_18_19 ~ PerCnpyOpenTotal + Sub_mean + ShrubCover_mean, data = Can_Cover_Surv,
            family=quasibinomial(logit))
summary(Mod4)

#model w/canopy photos and shrub cover
Mod5 <- glm(Surv_18_19 ~ PerCnpyOpenTotal + ShrubCover_mean, data=Can_Cover_Surv,
            family = quasibinomial(logit))
summary(Mod5)

#model w/canopy photos and abundance = density dependent mortality
Mod6 <- glm(Surv_18_19 ~ PerCnpyOpenTotal + Abund18, data=Can_Cover_Surv,
            family = quasibinomial(logit))
summary(Mod6)


#model w/Distance to entrance
Mod7 <- glm(Surv_18_19 ~ HQdist, data = Seed_Dist, family=quasibinomial(logit))
summary(Mod7)

#model abundance by distance to entrance
Mod8 <- lm(Abund19 ~ HQdist, data = Seed_Dist)
summary(Mod8) #significant p-value

#model survival w/shrub cover, abundance, and canopy cover
Mod9 <- glm(Surv_18_19 ~ ShrubCover_mean + Abund18 + PerCnpyOpenTotal, data=Can_Cover_Surv,
            family = quasibinomial(logit))
summary(Mod9)

#########################################################################
## Individual Sdlg Survival Modeling

#make variables numeric
Seed$Lvs18 <- as.numeric(Seed$Lvs18)
Seed$Hgt18 <- as.numeric(Seed$Hgt18)

#make smaller dataset
Seed2 <- Seed %>% 
  dplyr::select(PlotTag, Age18, Hgt18, Stat18_19, Stat19test, Stat18test,
                Damage18, Lvs18) %>%
  filter(!is.na(Stat18_19))

# initial model - Age, Height 18, Damage rating 18
M1 <- glm(Stat18_19 ~ Age18 + Damage18 + Hgt18, data = Seed2, family="binomial")
summary(M1)

# model 2 - Age, Leaf number 18, Damage rating 18
M2 <- glm(Stat18_19 ~ Age18 + Damage18 + Lvs18, data = Seed2, family="binomial")
summary(M2)

# model 3 - Age, Leaf number 18, Damage rating 18 (as category)
M3 <- glm(Stat18_19 ~ Age18 + as.factor(Damage18) + Lvs18, data = Seed2, family="binomial")
summary(M3)
stepAIC(M3) # check model fit
#used this model in HB presentation - July 2020 + ESA pres - August 2020


#Test model fit - compare single variable models
M4 <- glm(Stat18_19 ~ as.factor(Damage18), data = Seed2, family ="binomial")
summary(M4)
TukeyHSD(aov(M4)) # compare categories 

M5 <- glm(Stat18_19 ~ Lvs18, data = Seed2, family ="binomial")
summary(M5)

M6 <- glm(Stat18_19 ~ Age18, data = Seed2, family ="binomial")
summary(M6)

###############################################################################
################################################################
### Model predictions for individual sdlg survival

##########visualizing: binomial plots with prediction line
# specify formula
log_form <- y ~ 1 / (1 + exp(b * (x - a)))
##### age plot
log_model <- nls(formula = log_form, data = list(x = Seed2$Age18, y = Seed2$Stat18_19), start = list(a = 2, b = 1),control=list(maxiter=500))
cor(Seed2$Stat18_19, predict(log_model))
coef(log_model)
# add column with age predictions to dataset
Seed2$predAge <- 1 / (1 + exp(coef(log_model)["b"] * (Seed2$Age18 - coef(log_model)["a"])))

##### Lvs plot 
log_model_2 <- nls(formula = log_form, data = list(x = Seed2$Lvs18, y = Seed2$Stat18_19), start = list(a = 1.2, b = 0.3),control=list(maxiter=500))
cor(Seed2$Stat18_19, predict(log_model_2))
coef(log_model_2)
# add column with leaves predictions to dataset
Seed2$predLvs <- 1 / (1 + exp(coef(log_model_2)["b"] * (Seed2$Lvs18 - coef(log_model_2)["a"])))

##############################################################################################
##############################################################################################
## Plotting

#PLOT 1 
#plot of mean shrub cover at 20m bins by percent seedling survival 2018-2019
p1 <- ggplot(Cover_Surv, aes(x=ShrubCover_mean, y=Surv_18_19))+
  geom_point(size=3)+
  labs(x="Shrub Cover (%)",
       y="Seedling survival (proportion)")+
  theme_classic()+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
plot(p1)
#ggsave("Shrub_Surv1819.jpeg", p1, width=6, height=5)


#PLOT 2
#plot of mean available substrate (litter) by seedling abundance in 2019
ggplot(Cover_Surv, aes(x=Sub_mean, y=Abund19))+
  geom_point()+
  labs(x="% Available Substrate",
       y="Number of Seedlings")+
  theme_light()

#PLOT 3
#plot of 2019 seedling abundance by plot distance to entrance of valley
#create plot
p3 <- ggplot(Seed_Dist, aes(x=HQdist, y=Abund19))+
  geom_point(color = "royalblue3", size=3)+
  theme_classic()+
  labs(x= "Distance to HBEF entrance (m)",
       y = "Number of seedlings")+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
plot(p3)

#ggsave("Dist_Abund19.jpeg", p3, width=6, height=4)

#PLOT 5
#plot of seedling survival (18/19) at plot transect with distance to valley entrance
#create plot
p4 <- ggplot(Seed_Dist, aes(x=HQdist, y=Surv_18_19))+
  geom_point(size=3)+
  theme_classic()+
  labs(x= "Distance to HBEF entrance (m)",
       y = "Seedling survival (proportion)")+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14)) 
plot(p4)


#ggsave("Dist_Surv1819.jpeg", p4, width=6, height=4)


#PLOT 5
#plot of canopy cover (photos) at 20m bins by % sdlg survival 2018-2019
p5 <- ggplot(Can_Cover_Surv, aes(x=PerCnpyOpenTotal, y=Surv_18_19))+
  geom_point()+
  labs(x="Light Availability (%)",
       y="Seedling Survival (proportion)")+
  theme_classic()+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))

plot(p5)


#####################################################
## Plot Indiv Sdlg Survival by Characteristic variables 

#Plot 6
#plot age in 2018 by survival status
p6 <- ggplot(Seed2, aes(x=Age18, y=Stat18_19))+
  scale_y_continuous(breaks=c(0, 0.5, 1.0))+
  geom_jitter( height=0.1, color="turquoise4")+
  theme_classic()+
  geom_line(aes(y=predAge), size=1)+
  labs(y="Survival", x="Age")+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
plot(p6)

#save plot
#ggsave("Age_Surv1819.jpeg", p6, width=6, height=4)

#Plot 7
#plot number of leaves by survival status
p7<- ggplot(Seed2, aes(x=Lvs18, y=Stat18_19))+
  scale_y_continuous(breaks=c(0, 0.5, 1.0))+
  geom_jitter( height=0.1, color="turquoise4")+
  theme_classic()+
  geom_line(aes(y=predLvs), size=1)+
  labs(y="Survival", x="Number of leaves")+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
plot(p7)

# save plot
#ggsave("LeafNum_Surv1819.jpeg", p8, width=6, height=4)

#PLOT 8
# plot damage rating by sdlg survival proportion 2018/2019
#create dataset to plot damage rating and survival 
Damage <- Seed2 %>%
  filter(Stat18test==1)%>%
  group_by(Damage18) %>%
  summarize(Surv_18_19 = sum(Stat19test, na.rm=TRUE)/sum(Stat18test, na.rm=TRUE))%>%
  filter(!is.na(Damage18))

#create plot 
Damage$Damage18 <- as.factor(Damage$Damage18)
ggplot(Damage, aes(x=Damage18, y=Surv_18_19))+
  geom_bar(stat="Identity")+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Leaf damage rating", y="Seedling survival (proportion)")+
  theme_classic()+
  scale_x_discrete(breaks =c("0", "1", "2", "3", "4"), 
                   labels = c("0%", "1-25%", "26-50%","51-75%","76-100%"))+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))

#plot damage by survival summarized to plot level
ggplot(Seed_Bin3, aes(x=MeanDamage18, y=Surv_18_19))+
  geom_point()

#plot leaf number by survival summarized to plot level
ggplot(Seed_Bin3, aes(x=MeanLvs18, y=Surv_18_19))+
  geom_point()

#plot age by survival summarized to plot level
ggplot(Seed_Bin3, aes(x=MeanAge, y=Surv_18_19))+
  geom_point()

#plot abundance by survival summarized to 20m level
ggplot(Can_Cover_Surv, aes(x=Abund18, y=Surv_18_19))+
  geom_point()

#plot damage in 18 by survival status
ggplot(Seed2, aes(x=Damage18, y=Stat18_19))+
  geom_jitter()

##################################################################
##################################################################
###### Test Plots
### Bar plots for Lvs and Age
table(is.na(Lvs$Lvs18))

Lvs <- Seed2 %>%
  filter(Stat18test==1)%>% # filter out blank rows and outliers for Lvs18
  group_by(Lvs18) %>%
  summarize(Surv_18_19 = sum(Stat19test, na.rm=TRUE)/sum(Stat18test, na.rm=TRUE))

#plot leaf number by proportion survived 2018/2019            
ggplot(Lvs, aes(x=Lvs18, y=Surv_18_19))+
  geom_bar(stat="Identity")+
  theme_light()+
  labs(x="Number of leaves", y="Seedling survival (proportion)")

#create data subset for Age
Age <- Seed2 %>%
  filter(Stat18test==1, Age18!=0)%>%
  group_by(Age18) %>%
  summarize(Surv_18_19 = sum(Stat19test, na.rm=TRUE)/sum(Stat18test, na.rm=TRUE))

table(Seed2$Age)

#plot age by proportion survival 2018/2019
ggplot(Age, aes(x=Age18, y=Surv_18_19))+
  geom_bar(stat="Identity")+
  coord_cartesian(xlim=c(0,20))+
  labs(x="Age (years)", y="Seedling Survival (proportion")+
  theme_light()

#plot age by sdlg abundance
ggplot(Seed2, aes(x=Age18, y=Stat18_19))+
  geom_bar(stat="Identity")

#plot of total abundances of mature tree species - BA Prism data
ggplot(Prism_Abund, aes(x=reorder(Species, -Abund), y=Abund))+
  geom_bar(stat="identity")
