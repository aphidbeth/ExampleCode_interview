# HEAT TOLERANCE OF S.AVENAE CLONES - model building and stats
#-------------------------------------------------------------
# Created by BM
# Created March 2020
# Updated Feb 2022

# DATA: 
# Data collected by BM in a series of waterbath experiments in 2018. 
# individuals of different clones of Sitobion avenae were submerged 
# in a waterbath and the temperature ramped up from 18 - 45 degrees.
# Three points were observed for each individual: 
# temp.feed- The temperature at which aphids stopped feeding
# temp.run- The temperature at which aphids stopped running
# temp.twitch - The temperature at which aphids stopped twitching

# Heat tolerance with SA6 and SA1 merged (as same genotype)
# SA1 - SA3_A
# SA6 - SA3_B

#===============================================================
# HOUSEKEEPING
#==============

# Clear wd
rm(list=ls())

# Load libraries
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(lattice)
library(AICcmodavg)
library(broom)
library(jtools)


#Set working directory to either home or work
setwd("C:/Users/beth-/OneDrive - University of Aberdeen/Data/Physiology experiments")
#setwd("C:/Users/BethMoore/Dropbox/OneDrive/Data/Physiology experiments") # HOME #

#Load data
heat<- read.csv("HT_aug.csv")
str(heat)

#==============================================================
# DATA PREP 
#===========

heat$temp.feed.n<- as.numeric(as.character(heat$temp.feed))
heat$temp.feed.n
heat$temp.move.n<- as.numeric(as.character(heat$temp.move))
heat$temp.move.n

#Renaming SA1 and SA6 & IX so more logical
levels(heat$clone)[levels(heat$clone)=="SA1"] <- "SA3_A"
levels(heat$clone)[levels(heat$clone)=="SA6"] <- "SA3_B"
levels(heat$clone)[levels(heat$clone)=="IX"] <- "SA44"

#Combining SA1 and SA6 as they are the same genotype
levels(heat$genotype)[levels(heat$genotype)=="SA3_A"] <- "SA3"
levels(heat$genotype)[levels(heat$genotype)=="SA3_B"] <- "SA3"
# Adding the culture records

# ADDING LOCATION
heat$location <- ifelse(heat$clone == "SA7", "Kinnaber",
                        ifelse(heat$clone == "SA3_B", "Kinnaber",
                               ifelse(heat$clone == "SA9", "Kinnaber",
                                      ifelse(heat$clone == "SA3_A", "Kinnaber",
                                             ifelse(heat$clone == "SA11", "Old Meldrum",
                                                    ifelse(heat$clone == "SA12A", "Old Meldrum",
                                                           ifelse(heat$clone == "IX", "Ixworth", 
                                                                  NA  )))))))

# remove sample with unknown clone
heat<- heat[complete.cases(heat[,4]),]

# reorder factor levels so nice for plotting
heat$genotype <- factor(heat$genotype, levels=c("SA3", "SA44", "SA7", "SA9", "SA11", "SA12A"))
heat$clone <- factor(heat$clone, levels=c("SA3_A", "SA3_B", "SA44", "SA7", "SA9", "SA11", "SA12A"))


str(heat)
head(heat)
tail(heat)
View(heat)

# DATA EXPLORATION AVAIABLE IN glmmHTresistance script

#=====================================================================
# MODEL BUILDING
#==================

#FEEDING CESSATION
#------------------
# Fixed effects - Clone & Winged
# Random effects - Date

# model selection fixed effects
null.feed.clone <- lmer(temp.feed.n ~1 + (1|date_exp), data=heat)
feed.clone.2 <- lmer(temp.feed.n ~ clone + winged + (1|date_exp), data=heat)
feed.clone.3 <- lmer(temp.feed.n ~ winged + (1|date_exp), data=heat)
feed.clone.4 <-lmer(temp.feed.n ~ clone + (1|date_exp), data=heat)

# Model selection
AICc(null.feed.clone) #807.156
AICc(feed.clone.2) #794.2369
AICc(feed.clone.3) #804.226
AICc(feed.clone.4) #794.2859


# Model diagnostics
plot(feed.clone.4) #residuals
plot(feed.clone.4,sqrt(abs(residuals(.))) ~ fitted(.),type=c("p","smooth")) #scale-location
qqmath(feed.clone.4, auto.key=TRUE) # nae good
qqmath(feed.clone.4, col= heat$genotype, auto.key=TRUE) #nae good - 1 genotype

# aa <- augment(feed.clone.4)
# ggplot(aa,aes(clone,.resid))+
#   geom_boxplot()+coord_flip() # nae good # note this bit of code not working on 20/02/2023

summary(feed.clone.4)

#--- Transformation to limit left skew ---#

heat$trans1.temp.feed.n <- sqrt((max(heat$temp.feed.n, na.rm=TRUE)+1) - heat$temp.feed.n)
heat$trans2.temp.feed.n <-log10((max(heat$temp.feed.n, na.rm=TRUE)+1) - heat$temp.feed.n)
heat$trans3.temp.feed.n <-1/((max(heat$temp.feed.n, na.rm=TRUE)+1) - heat$temp.feed.n)
heat$trans4.temp.feed.n <-sqrt(heat$temp.feed.n)

hist(heat$temp.feed.n)
hist(heat$trans1.temp.feed.n)
hist(heat$trans2.temp.feed.n)
hist(heat$trans3.temp.feed.n)
hist(heat$trans4.temp.feed.n)

# TRANS 1
null.root.feed <- lmer(trans1.temp.feed.n ~1 + (1|date_exp), data=heat)
root.feed.2 <- lmer(trans1.temp.feed.n ~ clone + winged + (1|date_exp), data=heat)
root.feed.4 <-lmer(trans1.temp.feed.n ~ clone + (1|date_exp), data=heat)

AICc(null.root.feed) # 274.452
AICc(root.feed.2) # 288.0585
AICc(root.feed.4) #281.642

anova(null.root.feed, root.feed.2) # but 2 fits better than null here
anova(null.root.feed, root.feed.4) # As does 4
anova(root.feed.2, root.feed.4)

plot(root.feed.4) #residuals
plot(root.feed.4,sqrt(abs(residuals(.))) ~ fitted(.),type=c("p","smooth")) #scale-location
qqmath(root.feed.4, auto.key=TRUE) 
qqmath(root.feed.4, col= heat$genotype, auto.key=TRUE) 
# Just made the skew positive

# TRANS 2
null.root2.feed <- lmer(trans2.temp.feed.n ~1 + (1|date_exp), data=heat)
root2.feed.2 <- lmer(trans2.temp.feed.n ~ clone + winged + (1|date_exp), data=heat)
root2.feed.4 <-lmer(trans2.temp.feed.n ~ clone + (1|date_exp), data=heat)

AICc(null.root2.feed) # -56.4721
AICc(root2.feed.2) #-26.76104
AICc(root2.feed.4) # -36.29857

anova(null.root2.feed, root2.feed.2) # but 2 fits better than null here
anova(null.root2.feed, root2.feed.4) # As does 4
anova(root2.feed.2, root2.feed.4)

plot(root2.feed.4) #residuals
plot(root2.feed.4,sqrt(abs(residuals(.))) ~ fitted(.),type=c("p","smooth")) #scale-location
qqmath(root2.feed.4, auto.key=TRUE) 
qqmath(root2.feed.4, col= heat$genotype, auto.key=TRUE) 
# looks a bit better but very hard to interpret

summary(root2.feed.4)

# With clone as a fixed effect to look at the vairation 
fix.clone.feed <- lmer(temp.feed.n ~1 + (1|date_exp) + (1|clone), data=heat)
summary(fix.clone.feed)

#================================
# WITH INVERTED VARIABLES
#==============================

# Create inverted variables
max(heat$temp.feed.n, na.rm=TRUE) # 41
max(heat$temp.twitch, na.rm=TRUE) # 45
heat$invert.temp.feed <- 42 - heat$temp.feed.n
heat$invert.temp.twitch <- 46 - heat$temp.twitch

str(heat)

# Distributions
hist(heat$temp.feed.n)
hist(heat$invert.temp.feed.n)
hist(heat$temp.twitch)
hist(heat$invert.temp.twitch)

#=============================================================
# MODEL BUILDING - PART 2
#========================

#==================
# FEEDING CESSATION
#==================
# model building
# fixed effects: clone + winged
# random effects: date 

# Model fitting of fixed effects
invert.feed1 <- lmer(invert.temp.feed ~ clone + (1|date_exp), data=heat)

invert.feed2 <- lmer(invert.temp.feed ~ clone + winged + (1|date_exp), data=heat)

invert.feed3 <- lmer(invert.temp.feed ~ winged + (1|date_exp), data=heat)

invert.feed.null <-lmer(invert.temp.feed ~ 1 + (1|date_exp), data=heat)

AICc(invert.feed1) #clone -> 794.2859 
AICc(invert.feed2) # clone + winged -> 794.2369
AICc(invert.feed3) # winged -> 804.226
AICc(invert.feed.null) # Null -> 807.156



#model diagnostics
plot(invert.feed1) #residuals
plot(invert.feed1,sqrt(abs(residuals(.))) ~ fitted(.),type=c("p","smooth")) #scale-location
qqmath(invert.feed1, auto.key=TRUE) 
qqmath(invert.feed1, col= heat$genotype, auto.key=TRUE) 
# no better

# Can you use dharma?
library(DHARMa)
testDispersion(invert.feed1) 
simulationOutput2 <- simulateResiduals(fittedModel = e_mod2)
plot(simulationOutput2) # 

# Look at model summary
summary(invert.feed1)
summary(invert.feed2)

###--TEMP TWITCH --###

# model specification
# fixed effects: clone + winged
# random effects: date 

# model fitting of fixed effects:
null.twitch.clone <- lmer(temp.twitch ~1 + (1|date_exp), data=heat)
twitch.clone.2 <- lmer(temp.twitch ~ clone + winged + (1|date_exp), data=heat)
twitch.clone.3 <- lmer(temp.twitch ~  winged + (1|date_exp), data=heat)
twitch.clone.4 <-lmer(temp.twitch ~ clone + (1|date_exp), data=heat)

AICc(null.twitch.clone) #[1] 531.0599
AICc(twitch.clone.2) # [1] 538.822
AICc(twitch.clone.3) # 536.7133
AICc(twitch.clone.4) # [1] 537.0524 #null fits best with AIC

# model fitting of random effects: 
null.twitch.clone <- lmer(temp.twitch ~1 + (1|date_exp), data=heat)
null.twitch.clone2 <- glm(temp.twitch ~1 , data=heat)

AICc(null.twitch.clone) # better fit with random effect: 534.2164
AICc(null.twitch.clone2)# 656.454


# Model diagnostics
summary(null.twitch.clone)
plot(null.twitch.clone) #residuals
plot(null.twitch.clone,sqrt(abs(residuals(.))) ~ fitted(.),type=c("p","smooth")) #scale-location
qqmath(null.twitch.clone, auto.key=TRUE) # nae good
qqmath(null.twitch.clone, col= heat$genotype, auto.key=TRUE) #nae good - 1 genotype

aa <- augment(null.twitch.clone)
ggplot(aa,aes(clone,.resid))+
  geom_boxplot()+coord_flip() # nae good


# With clone as a fixed effect to look at the vairation 
fix.clone.twitch <- lmer(temp.twitch ~1 + (1|date_exp) + (1|clone), data=heat)
summary(fix.clone.twitch)



#############
### Plots ###
#############

### Plot of stopped feeding ###


heat %>%
  mutate(clone = factor(clone, levels=c("SA3_A", "SA3_B","SA44", "SA7", "SA9", "SA11", "SA12A"))) %>%
  ggplot(., aes(x=clone, y=temp.feed, fill=genotype)) + geom_boxplot()+
  labs(x= "Clone", y="Temperature stopped feeding")+
  ggtitle("Clonal differences in temp stopped feeding") + 
  theme(plot.title= element_text(hjust=0.5, face= "bold", size=17),
        axis.title.x=element_text(size=15), 
        axis.title.y=element_text(size=15, vjust=2),
        axis.text.x=element_text(size=10, vjust=.5), 
        axis.text.y=element_text(size=10)) 

### Plot of stopped moving ###
heat %>%
  mutate(clone = factor(clone, levels=c("SA3_A", "SA3_B","SA44", "SA7", "SA9", "SA11", "SA12A"))) %>%
  ggplot(., aes(x=clone, y=temp.twitch, fill=genotype)) + geom_boxplot()+
  labs(x= "Clone", y="Temperature stopped twitching")+
  ggtitle("Clonal differences in temp stopped twitching") + 
  theme(plot.title= element_text(hjust=0.5, face= "bold", size=17),
        axis.title.x=element_text(size=15), 
        axis.title.y=element_text(size=15, vjust=2),
        axis.text.x=element_text(size=10, vjust=.5), 
        axis.text.y=element_text(size=10))

### Violin plots ###
sample_size =heat %>% group_by(clone) %>% summarize(num=n())

# Plot feed.temp # FOR POWEPOINT (no labels or title)

png(file = "RawFeedingTemp.png", width = 1024, height = 768, units = 'px')
RFT<-heat %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(clone, "\n", "n=", num)) %>%
  ggplot(., aes(x=clone, y=temp.feed, fill=genotype)) + geom_violin(position= position_dodge(0.01), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.6) +
  geom_jitter(color="black", size=2, alpha=0.15, width=0.2) + 
  labs(x= "Clone", y="Temperature") + 
  theme(plot.title= element_text(hjust=0.5, face= "bold", size=17),
        axis.title.x=element_text(size=35), 
        axis.title.y=element_text(size=35, vjust=1.5),
        axis.text.x=element_text(size=20, vjust=2.0), 
        axis.text.y=element_text(size=20),
        legend.position="none")
print(RFT)
dev.off()

# Plot twitch.temp
png(file = "RawTwitchingTemp.png", width = 1024, height = 768, units = 'px')
RTT<-heat %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(clone, "\n", "n=", num)) %>%
  ggplot(., aes(x=clone, y=temp.twitch, fill=genotype)) + geom_violin(position= position_dodge(0.01), width=1) +
  geom_boxplot(width=0.1, color="black", alpha=0.6) +
  geom_jitter(color="black", size=2, alpha=0.15, width=0.2) + 
  labs(x= "Clone", y="Temperature") + 
  theme(plot.title= element_text(hjust=0.5, face= "bold", size=17),
        axis.title.x=element_text(size=35), 
        axis.title.y=element_text(size=35, vjust=1.5),
        axis.text.x=element_text(size=20, vjust=2.0), 
        axis.text.y=element_text(size=20),
        legend.position="none")
print(RTT)
dev.off()
