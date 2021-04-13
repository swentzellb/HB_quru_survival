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

Cover_Sub <- read.csv("HB_VW_Oak_Trans_ShrubCover_2020_CLEAN.csv",header=TRUE,  fileEncoding = "UTF-8-BOM") 
#refers to canopy and shrub cover and substrate

PlotCoord <- read.csv("Plot_Table_Basic.csv")
#refers to plot coordinate - or distance dataset


#####################################################################################################################
#####################################################################################################################
###Data Cleaning - Cover Substrate

#separate out and create new identifying columns
Cover_Sub <- Cover_Sub %>% 
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


#######################################################################################################
#######################################################################################################
# Calculate Survival between 2018 and 2019
#######################################################################################################
# make unique column for seedling ID
seedIDs <- mutate(Seed, seedID = paste(VW,SdlgNum,TransDir,Distance,sep="-"))


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
Seed <- dplyr::select(seedIDs, -SdlgSpp, -SdlgNum, -PointType,-Along, -TransEW, -FromTrans, -NSline)





# create test variables equal to seedling status variables for each yr
Seed$Stat19test <- Seed$Stat19
Seed$Stat18test <- Seed$Stat18

# change PD and NAs in 2019 column to 0
Seed$Stat19test[Seed$Stat19=="PD"| is.na(Seed$Stat19test)] <- 0

# change PD and NAs in 2018 column to 0
Seed$Stat18test[is.na(Seed$Stat18test) | Seed$Stat18=="PD" |
                  Seed$Stat18=="NF"] <- 0 

#check all NAs have been removed
table(is.na(Seed$Stat18test))

table(Seed$Stat19test)

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

# removes seedlings that were not present in 2019 (set to NA)
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

## Find total survival between 2018 and 2019
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
Seed_18_19 <- Seed %>%
  filter(Stat18test==1) %>%
  group_by(PlotTrans, DistBin) %>%
  summarize(Surv_18_19 = sum(Stat19test==1, na.rm=TRUE)/sum(Stat18test==1, na.rm=TRUE),
            Abund18 = sum(Stat18test==1), 
            Abund19 = sum(Stat19test, na.rm=TRUE)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))


# summarize seedling survival between 2019 and 2020 to 20m bins
Seed_19_20 <- Seed %>%
  filter(Stat19test==1) %>%
  group_by(PlotTrans, DistBin) %>%
  summarize(Surv_19_20 = sum(Stat20test==1, na.rm=TRUE)/sum(Stat19test==1, na.rm=TRUE),
            Abund18 = sum(Stat19test==1), 
            Abund19 = sum(Stat20test, na.rm=TRUE)) %>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))



#summarize Shrub Cover to 20m bins
Shrub_Bin <- Cover_Sub %>%
  group_by(PlotTag, PlotTrans, DistBin)%>%
  summarize(ShrubCover_mean = mean(ShrubCover, na.rm=TRUE))%>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))


#################################################################
#######################################################################################
###Joining
#join Cover subset data with Seedling Survival dataset (binned at 20m)
Shrub_Surv <- left_join(Shrub_Bin, Seed_18_19, by="PlotTransBin")

Shrub_Surv2 <- left_join(Shrub_Bin, Seed_19_20, by="PlotTransBin")


#model w/  shrub cover 18/19
Mod1 <- glm(Surv_18_19 ~ ShrubCover_mean, data = Shrub_Surv, family=quasibinomial)#generalized linear model with quasibinomial family
summary(Mod1)# show model output
# YOU POSSIBLY DID THIS WRONG --- 18/19 is ONLY THE ENVIRO TRANSECTS


#model w/ shrub cover 19/20
Mod2 <- glm(Surv_19_20 ~ ShrubCover_mean, data = Shrub_Surv2, family=quasibinomial(logit))#generalized linear model with quasibinomial family
summary(Mod2)# show model output


#PLOT 1 
#plot of mean shrub cover at 20m bins by percent seedling survival 2018-2019
ggplot(Shrub_Surv, aes(x=ShrubCover_mean, y=Surv_18_19))+
  geom_point(size=3)+
  labs(x="Shrub Cover (%)",
       y="Seedling survival (proportion)")+
  theme_classic()+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))


plot(p1)
#ggsave("Shrub_Surv1819.jpeg", p1, width=6, height=5)

#plot of mean shrub cover at 20m bins by percent seedling survival 2018-2019
ggplot(Shrub_Surv2, aes(x=ShrubCover_mean, y=Surv_19_20))+
  geom_point(size=3)+
  labs(x="Shrub Cover (%)",
       y="Seedling survival (proportion)")+
  theme_classic()+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))


#################################################################
##################################################################
#correlation tests

#comparing shrub cover to survival
cor.test(Shrub_Surv$ShrubCover_mean, Shrub_Surv$Surv_18_19) #p-value 0.01, cor=0.22

#comparing shrub cover to survival
cor.test(Shrub_Surv2$ShrubCover_mean, Shrub_Surv2$Surv_19_20) #p-value 0.28, cor=0.10



