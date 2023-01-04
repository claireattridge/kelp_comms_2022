## Barkley Sound RLS photo quadrats - substrate processing
## Author: Claire Attridge
## Origin date: Dec 2022

# Load base packages
library(tidyverse)
library(MetBrewer)


#### Similarity analysis ----

### Cleaning up the data

# Reading in the raw data
raw <- read.csv("./MSc_data/Data_new/BenthicCover_TestingSheet_CMA.csv", header=F)
names(raw) <- raw[2,] # Making second row the header names
raw <- raw[-c(1:2),] # Removing first two rows

# Cleaning up column types
raw <- raw %>%
  mutate(across(.cols=9:16, .fns=as.numeric)) %>%
  mutate(SiteName = as.factor(SiteName), TotalPts = as.factor(TotalPts), Depth_m = as.numeric(Depth_m), Date_photos = as.POSIXct(Date_photos, format="%m/%d/%Y")) %>%
  dplyr::select(-c(Processor, ImageID, ImageWidth_cm, ImageLength_cm, AreaPerPt))

# Calculating proportion cover types
prop <- raw %>%
  rowwise() %>%
  mutate(Pcanopy = (Allcanopy/Total),
         Punder = (Allunderstory/Total),
         Phard = (Hardbottom/Total),
         Psand = (Sand/Total),
         Pbiol = (Biological/Total),
         Pother = (Other/Total)) %>%
  ungroup() 

# Calculating averages by pt number
ave <- prop %>%
  group_by(TotalPts) %>%
  summarise(Pcanopy = mean(Pcanopy, na.rm=T),
            Punder = mean(Punder, na.rm=T),
            Phard = mean(Phard, na.rm=T),
            Psand = mean(Psand, na.rm=T),
            Pbiol = mean(Pbiol, na.rm=T),
            Pother = mean(Pother, na.rm=T))

# 
# # Selecting for new dataframes
# canopy <- subset(prop, select=c(1,3,4,12)) %>%
#   pivot_wider(names_from = TotalPts, values_from = Pcanopy)
# under <- subset(prop, select=c(1,3,4,13)) %>%
#   pivot_wider(names_from = TotalPts, values_from = Punder)
# hard <- subset(prop, select=c(1,3,4,14)) %>%
#   pivot_wider(names_from = TotalPts, values_from = Phard)


wider <- t(ave) # Transposing columns to rows
colnames(wider) <- wider[1,] # Setting first row values as column names 
wider <- wider[-1, ]
wider <- as.data.frame(wider)
wider <- wider %>% # Converting the columns from character to numeric
  mutate(across(everything(), as.numeric))
  

### Autocovariance fxn for similarity



dat <- acf(ave$Phard, plot=FALSE, na.action=na.omit)
acf_hard <- data.frame(lag = dat$lag[,,1], acf =  dat$acf[,,1])
acf_hard

dat <- acf(ave$Psand, plot=FALSE, na.action=na.omit)
acf_hard <- data.frame(lag = dat$lag[,,1], acf =  dat$acf[,,1])
acf_hard

test <- ave %>%
  mutate(ac = list(acf(ave$Phard, na.action=na.pass)))

new <- raw %>%
  group_by(TotalPts) %>%
  summarise(ac = list(acf(raw$Allunderstory, type="correlation", na.action=na.pass)))


#### Percent cover processing ----


### WILL PROBABLY NEED THIS FOR LATER SUBSTRATE PROCESSING...

# # sorting data categories into broader groupings
# rawgrp <- raw %>%
#   rowwise() %>% 
#   mutate(TotalPoints = sum(raw[,c(10:41)]),
#          Allbrown = DesmerestiaSpp + SargassumMuticum + EgregiaMenziesii + LaminariaSetchellii + AlariaMarginata + AgarumSpp + SaccharinaLatissima + CostariaCostata + PterygophoraCalifornica + PleurophycusGardneri + FucusSpp + HaplogloiaAndersonii + AnalipusJaponicus + UnknownBrown,
#          Allred = ChondracanthusExasperatus + MazzaellaSplendens + UnknownRed,
#          Allgreen = UlvaLactuca + UnknownGreen,
#          Allunderstory = DesmerestiaSpp + SargassumMuticum + EgregiaMenziesii + LaminariaSetchellii + AlariaMarginata + AgarumSpp + SaccharinaLatissima + CostariaCostata + PterygophoraCalifornica + PleurophycusGardneri + FucusSpp + HaplogloiaAndersonii + AnalipusJaponicus + UnknownBrown + ChondracanthusExasperatus + MazzaellaSplendens + UnknownRed + UlvaLactuca + UnknownGreen,
#          Allcanopy = MacrocystisPyrifera + NereocystisLuetkeana,
#          Allcoralline = CorallineStanding + CorallineStanding,
#          Hardbottom = Bedrock + Boulder + Cobble + Sand + CorallineCrustose, 
#          AnalyzedPoints = Total - Unknown - Object) 










