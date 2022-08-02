## Barkley Sound kelp forest data processing
## Author: Claire Attridge
## Origin date: August 2022

setwd("C:/Users/Claire/Desktop/MSc/Thesis")

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(MetBrewer)

#### Data cleaning ----

## Our kelp density data ##
dens <- read.csv("./MSc_data/Data_new/kelp_density_2022.csv") %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>% # Changing units to /m2 area
  rowwise() %>%
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup()

# Grouping / averaging for site specific density 
densgrp <- dens %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(KelpM = mean(Kelp), KelpSD = sd(Kelp))


## Our kelp height & biomass data ##
kelp <- read.csv("./MSc_data/Data_new/kelp_morphology_2022.csv")

# Special addition for plotting the raw data by density
kelpbydens <- merge(kelp, densgrp, by="SiteName", all=TRUE)

# Grouping/averaging for site specific height & biomass
kelpgrp <- kelp %>%
  mutate(SiteName = as.factor(SiteName), Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(HeightM = mean(Height_m, na.rm=T), HeightSD = sd(Height_m, na.rm=T))


### Joining the data together ###
kelpdat <- merge(densgrp, kelpgrp, by = "SiteName", all=TRUE)

#### Plotting the data ----

## Site specific density (summed Macro & Nereo) # Ordering by increasing density
d1 <- ggplot() +
  geom_point(data=dens, size=3, alpha=0.1, aes(x=reorder(SiteName,Kelp), y=Kelp)) +
  geom_pointrange(data=densgrp, size=0.6, aes(x=reorder(SiteName,KelpM), y=KelpM, ymin=KelpM-KelpSD, ymax=KelpM+KelpSD)) +
  scale_y_continuous(limits=c(-0.8,25), breaks=c(0,6,12,18,24)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("") + ylab("Density (/m2)")

## Site specific canopy height # Ordering by increasing density
c1 <- ggplot() + 
  geom_point(data=kelpbydens, size=3, alpha=0.1, aes(x=reorder(SiteName,KelpM), y=Height_m)) +
  geom_pointrange(data=kelpdat, size=0.6, aes(x=reorder(SiteName,KelpM), y=HeightM, ymin=HeightM-HeightSD, ymax=HeightM+HeightSD)) +
  scale_y_continuous(limits=c(-0.2,9), breaks=c(0,2,4,6,8)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("") + ylab("Canopy height (m)") +
  annotate("text", size=4, x=1, y=0.5, label="NA")


## Grouped multi panel plot ##

tiff(file="./MSc_code/draftplots/KelpMetrics.tiff", height = 8, width = 5.5, units = "in", res=500)

ggarrange(d1, c1, ncol=1, align="v", heights=c(1,1.4)) # Generating the paneled plot

dev.off() 


