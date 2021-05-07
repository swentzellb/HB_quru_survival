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
Litter <- read.csv("HB_VW_Oak_Trans_Enviro_Litter_CLEAN.csv") 
#refers to leaf litter layer
Cover_Sub <- read.csv("HB_VW_Oak_Trans_ShrubCover_2020_CLEAN.csv", header=TRUE,  fileEncoding = "UTF-8-BOM") 
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

# View(Litter)

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
Cover_Sub <- Cover_Sub %>% 
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
  group_by(PlotTrans, DistBin)%>%
  summarize(ShrubCover_mean = mean(ShrubCover, na.rm=TRUE), 
            Sub_mean = mean(Litter, na.rm=TRUE))%>%
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

#comparing distance to survival
cor.test(Seed_Dist$HQdist, Seed_Dist$Surv_18_19)

cor.test(Cover_Surv$Surv_18_19, Cover_Surv$Abund18)

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


Mod10 <- glm(Surv_18_19 ~ ShrubCover_mean, data=Can_Cover_Surv,
            family = quasibinomial(logit))
summary(Mod10)

nullMod <- glm(Surv_18_19 ~ 1, data=Can_Cover_Surv,
               family = quasibinomial(logit))
summary(nullMod)
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

# Mtest <- glm(Stat18_19 ~ Damage18, data = Seed2, family ="binomial")
# summary(Mtest)

## find means & sd survival for each damage category
Dam_0 <- Seed2 %>%
  filter(Damage18=="0") %>% 
  summarize(meanDam18 = mean(Stat18_19), 
            sdDam18 = sd(Stat18_19))

Dam_1 <- Seed2 %>%
  filter(Damage18=="1") %>% 
  summarize(meanDam1 = mean(Stat18_19), 
            sdDam1 = sd(Stat18_19))

Dam_2 <- Seed2 %>%
  filter(Damage18=="2") %>% 
  summarize(meanDam2 = mean(Stat18_19), 
            sdDam2 = sd(Stat18_19))

Dam_3 <- Seed2 %>%
  filter(Damage18=="3") %>% 
  summarize(meanDam3 = mean(Stat18_19), 
            sdDam3 = sd(Stat18_19))

Dam_4 <- Seed2 %>%
  filter(Damage18=="4") %>% 
  summarize(meanDam4 = mean(Stat18_19), 
            sdDam4 = sd(Stat18_19))

low <- c("0","1","2")

Dam_low <- Seed2 %>%
  filter(Damage18==low) %>% 
  summarize(meanDamlow = mean(Stat18_19), 
            sdDamlow = sd(Stat18_19))


high <- c("3","4")

Dam_high <- Seed2 %>%
  filter(Damage18==high) %>% 
  summarize(meanDamhigh = mean(Stat18_19), 
            sdDamhigh = sd(Stat18_19))



test <- aov(Stat18_19 ~ as.factor(Damage18), data = Seed2)
summary(test)

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
        axis.text = element_text(size=13))
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
        axis.text = element_text(size=13))
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
p8 <- ggplot(Damage, aes(x=Damage18, y=Surv_18_19))+
  geom_bar(stat="Identity")+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Leaf damage rating", y="Seedling survival (proportion)")+
  theme_classic()+
  scale_x_discrete(breaks =c("0", "1", "2", "3", "4"), 
                   labels = c("0%", "1-25%", "26-50%","51-75%","76-100%"))+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=13))
plot(p8)

library(cowplot)
# combine all into one plot
plot_grid(p6, p7, p8, rel_heights = c(1, 1), labels = "auto")

aligned_plots <- align_plots(p6, p7, p8, align=c("h"), greedy=TRUE)



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





#############################################################################
#############################################################################
## Leaf litter data
# summarize Leaf litter to 20m bins
Litter_Bin <- Litter %>%
  group_by(PlotTag, PlotTrans, DistBin)%>%
  summarize(Plottotal= sum(Total, na.rm=TRUE), #total number of leaves for plot
            ACRUsum = sum(ACRU, na.rm=TRUE), # sum of leaves by species for plot
            ACRUper = ACRUsum/Plottotal,     # percentage of leaves by species for plot
            FAGRsum = sum(FAGR, na.rm=TRUE),
            FAGRper = FAGRsum/Plottotal,
            BEALsum = sum(BEAL, na.rm=TRUE),
            BEALper = BEALsum/Plottotal,
            ACSAsum = sum(ACSA, na.rm=TRUE),
            ACSAper = ACSAsum/Plottotal,
            ACPEsum = sum(ACPE, na.rm=TRUE),
            ACPEper = ACPEsum/Plottotal,
            BEPAsum = sum(BEPA, na.rm=TRUE),
            BEPAper = BEPAsum/Plottotal,
            FRAMsum = sum(FRAM, na.rm=TRUE),
            FRAMper = FRAMsum/Plottotal,
            POGRsum = sum(POGR, na.rm=TRUE),
            POGRper = POGRsum/Plottotal,
            VIALsum = sum(VIAL, na.rm=TRUE),
            VIALper = VIALsum/Plottotal)%>%
  mutate(PlotTransBin = paste(PlotTrans, DistBin, sep = "_"))  #add 20m bin identifier

#############################################################################
#############################################################################
#### Basal Area Prism data

# summarize Basal Area prism data: mature tree community to 20m bins
Prism_Bin <- Prism %>%
  group_by(PlotTransBin)%>%
  summarize(Plottotal = sum(Total),              # total number of trees per plot
            TSCAsum = sum(TSCA, na.rm=TRUE),     # number of trees by species per plot
            FAGRsum = sum(FAGR, na.rm=TRUE),
            ACRUsum = sum(ACRU, na.rm=TRUE),
            BEALsum = sum(BEAL, na.rm=TRUE),
            PIRUsum = sum(PIRU, na.rm=TRUE),
            BEPAsum = sum(BEPA, na.rm=TRUE),
            FRAMsum = sum(FRAM, na.rm=TRUE),
            ACSAsum = sum(ACSA, na.rm=TRUE),
            PISTsum = sum(PIST, na.rm=TRUE),
            ACSPsum = sum(ACSP, na.rm=TRUE),
            ACPEsum = sum(ACPE, na.rm=TRUE),
            POGRsum = sum(POGR, na.rm=TRUE), 
            TSCAper = TSCAsum/Plottotal,        # percentage of species per plot
            FAGRper = FAGRsum/Plottotal,
            ACRUper = ACRUsum/Plottotal,
            BEALper = BEALsum/Plottotal,
            PIRUper = PIRUsum/Plottotal,
            BEPAper = BEPAsum/Plottotal,
            FRAMper = FRAMsum/Plottotal,
            ACSAper = ACSAsum/Plottotal,
            PISTper = PISTsum/Plottotal,
            ACSPper = ACSPsum/Plottotal,
            ACPEper = ACPEsum/Plottotal,
            POGRper = POGR, na.rm=TRUE)


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
