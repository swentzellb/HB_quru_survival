# Sage Wentzell-Brehme
# June 18, 2020
# Summer Science Research Program
# Combining Hubbard Brook historical Red Oak seedling survival and health data (2011-2019)
# with environmental data from 2019

# updated Sept 8, 2020
# thesis data analysis
# preliminary data cleaning and subsetting for 
# Combining: Hubbard Brook historical QURU seedling survival and health data (2011-2019)
# and Hubbard Brook environmental data from 2019 on 10 subset transects in hemlock-dom areas

#load packages
library(tidyverse)
library(tidyr)
library(lme4)
library(MASS)
library(raster)

#set working directory
setwd("~/Thesis/HB_quru_survival/data")

#read in csv files
Seed <- read.csv("HB_VW_Oak_Trans_Survival_CLEAN.csv") 
Seed$VW <- paste0("p", Seed$VW) # add p to beginning of VW column (same as PlotTag)
# historical survival data
Litter <- read.csv("HB_VW_Oak_Trans_Enviro_Litter_CLEAN.csv") 
#refers to leaf litter layer
Cover_Sub <- read.csv("HB_VW_Oak_Trans_Enviro_Cover_Substrate_CLEAN.csv") 
#refers to canopy and shrub cover and substrate
Prism <- read.csv("HB_VW_Oak_Trans_Enviro_Prism_CLEAN.csv")
#refers to basal area prism - used to measure mature tree community
Can <- read.csv("HB_VW_Oak_Trans_CanPhotos_CLEAN.csv")
#refers to canopy photos 
Can$Plot <- paste0("p", Can$Plot) # add p to beginning of VW column (same as PlotTag)
PlotCoord <- read.csv("Plot_Table_Basic.csv")
#refers to plot coordinate - or distance dataset

################################################################################################################
################################################################################################################
###Data Cleaning - Litter 

#separate into Point and Plot Direction - N/S columns
#PlotDir = Plot Direction - N/S of transect line
Litter <- Litter %>% 
  separate(Point, c("Point", "PlotDir"), sep="_")

#separate into Distance along transect and Transect Direction - E/W columns
#by taking the last character in the string
#Distance - along transect line
#TransDir - Transect Direction to the East or West
Litter <- Litter %>% separate(Point, c("Distance", "TransDir"), sep = -1)

#separate and remove point type, plot off transect or on center line of transect, from distance along transect
Litter <- Litter %>% separate(Distance, c("PointType", "Distance"), sep = 1)

View(Litter)

#Binning
# Combine PlotTag & TransDir into one unique transect ID
Litter <- Litter %>%
  mutate(PlotTrans = paste(PlotTag, TransDir, sep = "_"))

# For each unique transect ID (PlotTrans), bin to 20m
TransBins <- seq(0, 100, by = 20)

# Make column with distance changed to bin number
Litter <- Litter %>%
  mutate(DistBin = findInterval(Distance, TransBins))

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

# create new bins and add identifying columns
Cover_Sub <- Cover_Sub %>%
  mutate(DistBin = findInterval(Distance, TransBins)) %>% # bin to 20m
  mutate(PlotTrans = paste(PlotTag, TransDir, sep = "_")) %>% 
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_")) #unique bin ID
  
#remove p and change PlotTag name to VW
# Cover_Sub <- Cover_Sub %>%
#   separate(PlotTag, c("PlotTag", "VW"), sep=1)

#######################################################################################################################
#######################################################################################################################
###Data Cleaning - Prism

#separate into Point and Plot Direction - N/S columns
Prism<- Prism %>% 
  mutate(Point2 = Point) %>%
  separate(Point2, c("Point2", "PlotDir"), sep="_") %>%
  separate(Point2, c("Distance", "TransDir"), sep = -1) %>%
  separate(Distance, c("PointType", "Distance"), sep = 1) %>%
  mutate(PlotPoint = paste(PlotTag, Point, sep = "_")) #unique point ID
#View(Prism)

########################################################################
#Binning - Prism
# Combine PlotTag & TransDir into one unique transect ID
Prism <- Prism %>%
  mutate(PlotTrans = paste(PlotTag, TransDir, sep = "_")) %>%
  mutate(DistBin = findInterval(Distance, TransBins)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_")) #unique bin ID

#########################################################################
#########################################################################
#Data Cleaning for Historical seedling dataset

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

#############################################################################
# individual seedling survival status from 2018 to 2019
#############################################################################

# add another column to Seed dataset equal to stat19
Seed$Stat18_19 <- Seed$Stat19test

# removes seedlings that were not present in 2018 (set to NA)
Seed$Stat18_19[Seed$Stat18test==0]<- NA

table(is.na(Seed$Stat18_19)) #should be 654 false

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

##############################################################################################
###############################################################################################
## Data Cleaning - canopy photos

#clean canopy photo dataset
Can <- Can %>%
  dplyr::select(Plot, Point, PerCnpyOpenTotal) %>%
  mutate(Distance = substr(Point, start=2, stop=3))%>%
  mutate(DistBin = findInterval(Distance, TransBins),
         Point2 = Point, 
         PlotPoint = paste(Plot, Point, sep="_")) %>% #add unique ID for point
  separate(Point2, c("x", "TransDir"), sep = -1)%>%
  mutate(PlotTransBin = paste(Plot, TransDir, DistBin, sep="_")) #add unique ID for bins

#subset smaller dataset for canopy photos
Can_2 <- Can %>% 
  dplyr::select(PlotTransBin, PerCnpyOpenTotal, PlotPoint)

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
Cover_Sub_Bin <- Cover_Sub %>%
  group_by(PlotTag, PlotTrans, DistBin)%>%
  summarize(CanCover_mean = mean(CanopyCover, na.rm=TRUE),
            ShrubCover_mean = mean(ShrubCover, na.rm=TRUE), 
            Sub_mean = mean(Litter, na.rm=TRUE))%>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))

#################################################################
#######################################################################################
###Joining
#join Cover subset data with Seedling Survival dataset (binned at 20m)
Cover_Surv <- left_join(Cover_Sub_Bin, Seed_Bin_Enviro, by="PlotTransBin")

#join Seedling Abundance at Plot Transect level with distances to valley entrance
Seed_Dist <- left_join(Seed_Bin2, CoordAdj, by="PlotTrans")

##############################################################################
#join new canopy_cover_substrate datasets (binned to 20m level) with survival

#subset cover and survival dataset
Cover_Surv2 <- Cover_Surv %>%
  dplyr::select(ShrubCover_mean, PlotTransBin, Surv_18_19, Sub_mean, Abund18)

# create new dataset with canopy cover, other environmental measures, and survival data
Can_Cover_Surv <- left_join(Cover_Surv2, Can_2, by = "PlotTransBin")

#################################################################################
##################################################################################
###Basal Area Prism data

#summarize all the total tree abundances by species
Prism_Abund <- Prism %>%
  summarise(TSCA = sum(TSCA, na.rm=TRUE),
            FAGR = sum(FAGR, na.rm=TRUE),
            ACRU = sum(ACRU, na.rm=TRUE),
            BEAL = sum(BEAL, na.rm=TRUE),
            PIRU = sum(PIRU, na.rm=TRUE),
            BEPA = sum(BEPA, na.rm=TRUE),
            FRAM = sum(FRAM, na.rm=TRUE),
            ACSA = sum(ACSA, na.rm=TRUE),
            PIST = sum(PIST, na.rm=TRUE),
            ACSP = sum(ACSP, na.rm=TRUE),
            ACPE = sum(ACPE, na.rm=TRUE),
            POGR = sum(POGR, na.rm=TRUE))

#reshape the dataframe
Prism_Abund <-gather(Prism_Abund, `TSCA`, `FAGR`, `ACRU`,`BEAL`,
                   `PIRU`, `BEPA`, `FRAM`, `ACSA`, `PIST`, `ACSP`,
                   `ACPE`, `POGR`,
                   key = "Species", value = "Abund")

#plot of total abundances of mature tree species - BA Prism data
ggplot(Prism_Abund, aes(x=reorder(Species, -Abund), y=Abund))+
  geom_bar(stat="identity")

Prism_Bin <- Prism %>%
  group_by(PlotTransBin)%>%
  summarize(TSCA = sum(TSCA, na.rm=TRUE),
            FAGR = sum(FAGR, na.rm=TRUE),
            ACRU = sum(ACRU, na.rm=TRUE),
            BEAL = sum(BEAL, na.rm=TRUE),
            PIRU = sum(PIRU, na.rm=TRUE),
            BEPA = sum(BEPA, na.rm=TRUE),
            FRAM = sum(FRAM, na.rm=TRUE),
            ACSA = sum(ACSA, na.rm=TRUE),
            PIST = sum(PIST, na.rm=TRUE),
            ACSP = sum(ACSP, na.rm=TRUE),
            ACPE = sum(ACPE, na.rm=TRUE),
            POGR = sum(POGR, na.rm=TRUE))

Prism_Bin2 <- Prism_Bin %>%
  group_by(PlotTransBin)%>%
  summarize(TSCAper = (sum(TSCA, na.rm=TRUE))/sum(TSCA:POGR),
            total= sum(TSCA:POGR))
  

