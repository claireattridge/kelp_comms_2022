## Barkley Sound kelp forest data processing
## Author: Claire Attridge
## Origin date: August 2022

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
kelp <- read.csv("./MSc_data/Data_new/kelp_morphology_2022.csv") %>%
  as.data.frame()

# Converting circumference to diameter
kelp <- kelp %>%
  rowwise() %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159)) %>%
  ungroup()
  
# Equation for converting from sub-bulb to biomass as per our Nereo sub-bulb model
formula <- function(x){
  (150.7966*(x)^2 -216.2721*(x) + 315.0124)
}

options(scipen=999) # Turning off scientific notation

# Statement for applying sub-bulb equation when biomass is absent
kelp$Biomass_g <- ifelse(is.na(kelp$Biomass_g), ifelse(is.na(kelp$Sub_diam_cm), NA, formula(kelp$Sub_diam_cm)), kelp$Biomass_g)

# Special addition for plotting the raw data by density
kelpbydens <- merge(kelp, densgrp, by="SiteName", all=TRUE)


# Grouping/averaging for site specific height & biomass
kelpgrp <- kelp %>%
  mutate(SiteName = as.factor(SiteName), Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(HeightM = mean(Height_m, na.rm=T), HeightSD = sd(Height_m, na.rm=T),
            BiomassM = mean(Biomass_g, na.rm=T), BiomassSD = sd(Biomass_g, na.rm=T))


### Joining the data together ###
kelpdat <- merge(densgrp, kelpgrp, by = "SiteName", all=TRUE) %>%
               merge(tempgrp) # Adding in the temp logger data


#### Plotting the data ----

## Site specific density (summed Macro & Nereo) # Ordered by increasing density
d1 <- ggplot() +
  geom_point(data=dens, size=3, alpha=0.1, aes(x=reorder(SiteName,Kelp), y=Kelp)) +
  geom_pointrange(data=kelpdat, size=0.7, aes(x=reorder(SiteName,KelpM), y=KelpM, ymin=KelpM-KelpSD, ymax=KelpM+KelpSD)) +
  scale_y_continuous(limits=c(-0.8,25), breaks=c(0,6,12,18,24)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("") + ylab("Density (/m2)")

## Site specific canopy height # Ordered by increasing density
c1 <- ggplot() + 
  geom_point(data=kelpbydens, size=3, alpha=0.1, aes(x=reorder(SiteName,KelpM), y=Height_m)) +
  geom_pointrange(data=kelpdat, size=0.7, aes(x=reorder(SiteName,KelpM), y=HeightM, ymin=HeightM-HeightSD, ymax=HeightM+HeightSD)) +
  # scale_y_continuous(limits=c(-0.4,9), breaks=c(0,2,4,6,8)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("") + ylab("Canopy height (m)") +
  annotate("text", size=4, x=1, y=0.5, label="NA")

## Site specific biomass # Ordered by increasing density
b1 <- ggplot() + 
  geom_point(data=kelpbydens, size=3, alpha=0.1, aes(x=reorder(SiteName,KelpM), y=Biomass_g)) +
  geom_pointrange(data=kelpdat, size=0.7, aes(x=reorder(SiteName,KelpM), y=BiomassM, ymin=BiomassM-BiomassSD, ymax=BiomassM+BiomassSD)) +
  # scale_y_continuous(limits=c(-0.4,9), breaks=c(0,2,4,6,8)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("") + ylab("Biomass (g)") +
  annotate("text", size=4, x=1, y=0.5, label="NA")

## Grouped multi panel plot ##

tiff(file="./MSc_code/draftplots/KelpMetrics.tiff", height = 10, width = 5.5, units = "in", res=500)

ggarrange(d1, c1, b1, ncol=1, align="v", heights=c(1,1,1.4)) # Generating the paneled plot

dev.off() 



#### Regressions of the data ----

kelpdat2 <- merge(kelpdat, rls_richness, by="SiteName", all=T)

ggplot(kelpdat2, aes(x=HeightM, y=species_richness)) +
  geom_point() +
  # geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T) + # Lm quadratic
  # geom_smooth(method="lm", formula = y ~ x, se=T) + # Lm
  geom_smooth(method="gam", formula = y ~ s(x, bs="cs", k=3), se=T) + # GAM
  theme_classic()
