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
library(broom.mixed)
library(lattice)


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
# Calculate Survival across year intervals
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
Seed <- dplyr::select(seedIDs, -SdlgSpp, -SdlgNum, -PointType, -Distance, 
                      -TransDir, -PlotDir,-Along, -TransEW, -FromTrans, -NSline)

# replace PD and NF seedling status to 0
Seed[Seed=="PD"] <- "0"
Seed[Seed=="NF"] <- "0"
Seed[Seed=="nf"] <- "0"
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


#####################################################
#####################################################
### Create figure of proportion sdlg survival for each year interval
#####################################################

# using seedInterval data frame which includes seedling data, 
# indiv characteristics
# & survival of seedlings across year intervals
seedSurvival <- seedIntervalAge %>%
  group_by(interval) %>%
  summarise(propSurv = sum(survival==1)/n())
      
test <- seedIntervalAge %>%
  filter(interval=="2019-2020")

ggplot(seedSurvival, aes(x=interval, y=propSurv))+
  geom_bar(stat="identity", fill="skyblue4")+
  ylim(0,1)+
  xlab("Survival Year Interval")+
  ylab("Seedlings Survived (proportion)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=290, size = 12),
        axis.text.y = element_text(size = 12),
          axis.title = element_text(size=14))


#########################################################
#########################################################
#Statistical Analysis
#generalized linear models

#change variable classes to factor
seedIntervalAge$seedDamage <- as.factor(seedIntervalAge$seedDamage)
seedInterval$interval <- as.factor(seedInterval$interval)

#############
#Create a model for survival by Leaf number, Damage, Live/Dead Branch ratio, and Age 
Mod_1 <- glm(survival ~ leafNumber + seedDamage + brchLvD + yearsAlive, data=seedIntervalAge, 
             family = quasibinomial(logit))
summary(Mod_1)

# Find model RMSE
RSS <- c(crossprod(Mod_1$residuals))
MSE <- RSS / length(Mod_1$residuals)
RMSE <- sqrt(MSE)
#RMSE = 1.12

############
#Model for survival by Leaf Number
Mod_2 <- glm(survival ~ leafNumber, data=seedIntervalAge, family=quasibinomial(logit))
summary(Mod_2)

# Find model RMSE
RSS <- c(crossprod(Mod_2$residuals))
MSE <- RSS / length(Mod_2$residuals)
RMSE <- sqrt(MSE)
#RMSE = 770201.7

############
#Model for survival by Leaf Damage
Mod_3 <- glm(survival ~ seedDamage + interval, data=seedIntervalAge, family=quasibinomial(logit))
summary(Mod_3)

# Find model RMSE
RSS <- c(crossprod(Mod_3$residuals))
MSE <- RSS / length(Mod_3$residuals)
RMSE <- sqrt(MSE)
#RMSE = 4.15

############
#Model for survival by Live/Dead Branch Ratio
Mod_4 <- glm(survival ~ brchLvD, data=seedIntervalAge, family=quasibinomial(logit))
summary(Mod_4)

# Find model RMSE
RSS <- c(crossprod(Mod_4$residuals))
MSE <- RSS / length(Mod_4$residuals)
RMSE <- sqrt(MSE)
#RMSE = 20.39

############
#Model for survival by Age
Mod_5 <- glm(survival ~ yearsAlive, data=seedIntervalAge, family=quasibinomial(logit))
summary(Mod_5)
#yearsAlive is significant p-value<0.001

# Find model RMSE
RSS <- c(crossprod(Mod_5$residuals))
MSE <- RSS / length(Mod_5$residuals)
RMSE <- sqrt(MSE)
#RMSE = 2.69

########################################################################
#Add year as a random effect

seedIntervalAge <- seedIntervalAge %>%
  mutate(log_yearsAlive = log(yearsAlive))

summary(seedIntervalAge$log_yearsAlive)
##########
#Model for survival by Leaf number and Age with Year interval as a random effect
M1 <- glmer(survival ~ leafNumber + yearsAlive + (1|interval), 
            data=seedIntervalAge, family="binomial")
summary(M1)
ranef(M1)
glance(M1) # from broom.mixed package, returns table with AIC, BIC, & df of residual
tidy(M1)
coef(M1)

#make random effects data into a data frame
random <- data.frame(ranef(M1))

#plot random effects by interval
require(lattice)
dotplot(ranef(M1, condVar=TRUE))





# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M1))) #residual sum of squares
MSE <- RSS / length(residuals(M1)) #mean squared error
RMSE <- sqrt(MSE)
#RMSE = 0.376

##########
#Model for survival by Leaf number and Live/Dead Branch ratio 
#with Year interval as a random effect
M2 <- glmer(survival ~ leafNumber + brchLvD + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M2)
ranef(M2)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M2))) #residual sum of squares
MSE <- RSS / length(residuals(M2)) #mean squared error
RMSE <- sqrt(MSE)
#RMSE = 0.199

##########
#Model for survival by Leaf number, yearsAlive, and Live/Dead Branch ratio 
#with Year interval as a random effect
M3 <- glmer(survival ~ leafNumber + yearsAlive + brchLvD + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M3)
# yearsAlive still not significant
ranef(M3)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M3))) #residual sum of squares
MSE <- RSS / length(residuals(M3)) #mean squared error
RMSE <- sqrt(MSE)
#RMSE = 0.199

########## M4
#Model for survival by Leaf number, yearsAlive, Live/Dead Branch ratio, and Seed Damage 
#with Year interval as a random effect
M4 <- glmer(survival ~ leafNumber + yearsAlive + brchLvD + seedDamage + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M4)
ranef(M4)

############
#check distribution of seed Damage variable by year
table(seedIntervalAge$seedDamage)

seedDamageInt <- seedIntervalAge %>%
  group_by(interval) %>%
  count(seedDamage)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M4))) #residual sum of squares
MSE <- RSS / length(residuals(M4)) #mean squared error
RMSE <- sqrt(MSE)
#RMSE = 0.124

########## M5
#Model for survival by Leaf Damage
#with Year interval as a random effect
M5 <- glmer(survival ~ seedDamage + (1|interval), 
            data=seedIntervalAge, family = "binomial", 
            control = glmerControl(optimizer ="Nelder_Mead"))
summary(M5)
ranef(M5)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M5))) #residual sum of squares
MSE <- RSS / length(residuals(M5)) #mean squared error
RMSE <- sqrt(MSE)
#RMSE = 0.248


########## M6
#Model for survival by Leaf Number, Age, and Leaf Damage
#with Year interval as a random effect
M6 <- glmer(survival ~ leafNumber + yearsAlive + seedDamage + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M6)
ranef(M6)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M6))) #residual sum of squares
MSE <- RSS / length(residuals(M6)) #mean squared error
RMSE <- sqrt(MSE)

########## M7
#Model for survival by Leaf Number
#with Year interval as a random effect
M7 <- glmer(survival ~ leafNumber + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M7)
ranef(M7)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M7))) #residual sum of squares
MSE <- RSS / length(residuals(M7)) #mean squared error
RMSE <- sqrt(MSE)
RMSE


########## M8
#Model for survival by Age
#with Year interval as a random effect
M8 <- glmer(survival ~ yearsAlive + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M8)
ranef(M8)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M8))) #residual sum of squares
MSE <- RSS / length(residuals(M8)) #mean squared error
RMSE <- sqrt(MSE)
RMSE

########## M9
#Model for survival by Live/Dead Branch Ratio
#with Year interval as a random effect
M9 <- glmer(survival ~ brchLvD + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M9)
ranef(M9)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M9))) #residual sum of squares
MSE <- RSS / length(residuals(M9)) #mean squared error
RMSE <- sqrt(MSE)
RMSE


########## M10
#Model for survival by Age and Leaf Damage
#with Year interval as a random effect
M10 <- glmer(survival ~ yearsAlive + seedDamage + (1|interval), 
            data=seedIntervalAge, family = "binomial")
summary(M10)
ranef(M10)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M10))) #residual sum of squares
MSE <- RSS / length(residuals(M10)) #mean squared error
RMSE <- sqrt(MSE)
RMSE

########## M11
#Model for survival by leaf Number and Leaf Damage
#with Year interval as a random effect
M11 <- glmer(survival ~ leafNumber + seedDamage + (1|interval), 
             data=seedIntervalAge, family = "binomial")
summary(M11)
ranef(M11)

# Calculate model Root Mean Squared Error (RMSE)
RSS <- c(crossprod(residuals(M11))) #residual sum of squares
MSE <- RSS / length(residuals(M11)) #mean squared error
RMSE <- sqrt(MSE)
RMSE


########## M11
#Model for survival by leaf Number and Leaf Damage
#with Year interval as a random effect
M12 <- glmer(survival ~ leafNumber + log(yearsAlive) + (1|interval), 
             data=seedIntervalAge, family = "binomial")
summary(M11)
ranef(M11)

###################################################
##Correlation tests

#correlation between leaf number and branch ratio
cor.test(seedIntervalAge$leafNumber, seedIntervalAge$brchLvD)
#correlated p-value <0.01 and cor = 0.17

#correlation between age and branch ratio
cor.test(seedIntervalAge$yearsAlive, seedIntervalAge$brchLvD)
#correlated p-value <0.01 and cor = -0.23

#correlation between age and leaf number
cor.test(seedIntervalAge$yearsAlive, seedIntervalAge$leafNumber)
#correlated p-value <0.01 and cor = 0.2

#correlation between branch ratio and survival
cor.test(seedIntervalAge$brchLvD, seedIntervalAge$survival)
#correlated p-value <0.01 and cor = 0.106

#correlation between leaf number and survival
cor.test(seedInterval$leafNumber, seedInterval$survival)
#correlated p-value < 2.2e-16 and cor = 0.499

#correlation between Age and survival
cor.test(seedIntervalAge$yearsAlive, seedIntervalAge$survival)
#correlated p-value < 0.01 and cor = 0.1

#correlation between Age and seedling damage
cor.test(seedIntervalAge$yearsAlive, as.numeric(seedIntervalAge$seedDamage))
#correlated p-value =.002 and cor = 0.06

#correlation between survival and seedling damage
cor.test(seedIntervalAge$survival, as.numeric(seedIntervalAge$seedDamage))
#correlated p-value < 0.001 and cor = 0.33


##########################################################################
### Plots
##########################################################################

#plot Distribution of Seedling Ages for each Year interval
boxplot(yearsAlive ~ interval, data = seedIntervalAge, xlab = "year",
        ylab = "age", main = "Distribution of seedling ages each year")

table(seedIntervalAge$yearsAlive)

test2 <- seedIntervalAge %>%
  subset(seedIntervalAge$yearsAlive < 11)

#plot Distribution of Seedling Ages for each Year interval with outliers removed
# outliers defined as seedlings age 11 years and older
boxplot(yearsAlive ~ interval, data = test2, xlab = "year",
        ylab = "age", main = "Distribution of seedling ages each year")

boxplot(log_yearsAlive ~ interval, data = seedIntervalAge, xlab = "year",
        ylab = "age", main = "Distribution of seedling ages each year")




###############################################################################
################################################################
### Model predictions for individual sdlg survival

##########visualizing: binomial plots with prediction line
# specify formula
log_form <- y ~ 1 / (1 + exp(b * (x - a)))
##### age plot
log_model <- nls(formula = log_form, data = list(x = seedIntervalAge$yearsAlive, y = seedIntervalAge$survival), start = list(a = 2, b = 1),control=list(maxiter=500))
cor(seedIntervalAge$survival, predict(log_model))
coef(log_model)
# add column with age predictions to dataset
seedIntervalAge$predAge <- 1 / (1 + exp(coef(log_model)["b"] * (seedIntervalAge$yearsAlive - coef(log_model)["a"])))


##### Lvs plot
log_model_2 <- nls(formula = log_form, data = list(x = seedIntervalAge$leafNumber, y = seedIntervalAge$survival), start = list(a = 1.5, b = .9),control=list(maxiter=500))
cor(seedIntervalAge$survival, predict(log_model_2))
coef(log_model_2)
# add column with age predictions to dataset
seedIntervalAge$predLvs <- 1 / (1 + exp(coef(log_model_2)["b"] * (seedIntervalAge$leafNumber - coef(log_model_2)["a"])))

seedLvsPlot <- seedIntervalAge 


#Plot 6
#plot yearsAlive (or age of seedling) by survival status 
#multi year model - includes 2011-2020
ggplot(seedIntervalAge, aes(x=yearsAlive, y=survival, na.rm=TRUE))+
  scale_y_continuous(breaks=c(0, 0.5, 1.0))+
  geom_jitter( height=0.1, color="turquoise4")+
  theme_classic()+
  geom_line(aes(y=predAge), size=1)+
  labs(y="Survival", x="Age")+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
#plot(p6)


#Plot
#plot yearsAlive (or age of seedling) by survival status 
#multi year model - includes 2011-2020
ggplot(seedIntervalAge, aes(x=leafNumber, y=survival, na.rm=TRUE))+
  scale_y_continuous(breaks=c(0, 0.5, 1.0))+
  geom_jitter( height=0.1, color="turquoise4")+
  theme_classic()+
  geom_line(aes(y=predLvs), size=1)+
  labs(y="Survival", x="Number of leaves")+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
#plot(p7)


