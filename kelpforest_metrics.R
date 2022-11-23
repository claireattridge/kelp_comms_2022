## Barkley Sound kelp forest data processing
## Author: Claire Attridge
## Origin date: August 2022

# Loading base packages
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(MetBrewer)

#### Kelp density cleaning ----

## Our kelp density data
dens <- read.csv("./MSc_data/Data_new/kelp_density_2022.csv") %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>% # Changing units to /m2 area
  rowwise() %>%
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup()

# Grouping / averaging for site specific density 
densgrp <- dens %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(KelpM = mean(Kelp), KelpSD = sd(Kelp))


#### Kelp height & biomass cleaning ----

## Our kelp height and biomass data
kelp <- read.csv("./MSc_data/Data_new/kelp_morphology_2022.csv") %>%
  as.data.frame()

# Converting circumference to diameter
kelp <- kelp %>%
  rowwise() %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159)) %>%
  ungroup()
  
# Equation for converting from sub-bulb to biomass as per our Nereo sub-bulb model (all sites combined)
# See 'nereo_biomass.R'
formula <- function(x){
  (150.7966*(x)^2 -216.2721*(x) + 315.0124)
}

options(scipen=999) # Turning off scientific notation

# Statement for applying sub-bulb equation when biomass is absent (i.e., when it needs to happen)
kelp$Biomass_g <- ifelse(is.na(kelp$Biomass_g), ifelse(is.na(kelp$Sub_diam_cm), NA, formula(kelp$Sub_diam_cm)), kelp$Biomass_g)


#
# A 'sidenote' dataframe for plotting the raw data by density order
kelpbydens <- merge(kelp, densgrp, by="SiteName", all=TRUE)
#


# Grouping/averaging for site specific height & biomass
kelpgrp <- kelp %>%
  mutate(SiteName = as.factor(SiteName), Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(HeightM = mean(Height_m, na.rm=T), HeightSD = sd(Height_m, na.rm=T),
            BiomassM = mean(Biomass_g, na.rm=T), BiomassSD = sd(Biomass_g, na.rm=T))


#### Kelp area cleaning ----

library(scales)
library(sf)
library(ggsn)
library(MetBrewer)
library(maptools)
library(rgeos)

### Prepping the backing map

## setting projections at outset 
proj <- st_crs(3005) # ESPG for BC/Albers 
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")

## setting map extent for Barkley Sound 
ymax <- 48.922
ymin <- 48.80
xmax <- -125.05
xmin <- -125.26

## making corners for the area of interest 
corners <- st_multipoint(matrix(c(xmax,xmin,ymax,ymin),ncol=2)) %>% 
  st_sfc(crs=latlong) %>%
  st_sf() %>%
  st_transform(proj) 
plot(corners)


### Loading in the kelp forest shapefiles

# get all files with the .shp extension from the correct folder (i.e., the QGIS cleaned versions)
path <- "./MSc_data/Data_new/Kelp_area/CleanGPX" # set path to the folder

n <- list.files(path, pattern = "*shp", full.names=FALSE) # pull only the relevant names

shps <- list.files(path, pattern = "*shp", full.names=TRUE) # pull full file location names
lshps <- as.list(shps) # turn character vector into list

names(lshps) <- n # assign relevant names to elements of the list

f <- function(x){
  st_read(x, crs=latlong) # base crs of the GPS used to create the shapefiles is WGS84 
}

allshps <- lapply(lshps, f) # apply function to read in as sf objects with crs set as WGS84


list2env(allshps, envir=.GlobalEnv) # bring each unique sf object to the working environment for plotting             


### Calculating area of each kelp forest shape file

# I know this is terrible, still working on code to automate/loop st_area() through the .shp files
# But, at least it works! Keep eyes out for updates to this code. 
a <- st_area(allshps[[1]])
b <- st_area(allshps[[2]])
c <- st_area(allshps[[3]])
d <- st_area(allshps[[4]])
e <- st_area(allshps[[5]])
f <- st_area(allshps[[6]])
g <- st_area(allshps[[7]])
h <- st_area(allshps[[8]])
i <- st_area(allshps[[9]])
j <- st_area(allshps[[10]])
k <- st_area(allshps[[11]])
l <- st_area(allshps[[12]])
m <- st_area(allshps[[13]])
n <- st_area(allshps[[14]])
o <- st_area(allshps[[15]])
p <- st_area(allshps[[16]])
q <- st_area(allshps[[17]])
r <- st_area(allshps[[18]])
s <- st_area(allshps[[19]])
t <- st_area(allshps[[20]])
u <- st_area(allshps[[21]])
v <- st_area(allshps[[22]])
w <- st_area(allshps[[23]])

# Making dataframe of kelp forest areas
area.frame <- data.frame(SiteName = names(allshps), Area_m2 = cbind(c(a,d,b,c,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w))) %>%
  mutate(SiteName = as.factor(SiteName), Area_m2 = as.numeric(Area_m2))

############### TRYING TO FIGURE OUT HOW TO AUTOMATE THIS OVER THE LIST EGAUYGDUGHAIJIJIJXIAIEUHW

f2 <- function(polygon){
  st_area(polygon) #
}


test <- vector('list', length(allshps))
for (i in seq_along(allshps)) {
  x <- st_area(allshps[[i]])
  test[i] <- x
}
testdf <- do.call(rbind, test)


# Empty List for Centroids
test <- list()

# For Loop
for (i in seq_along(allshps)) {
    sf::st_area(allshps[[i]])
}

################ Anyways, moving on for now...

#### Joining the data together ----


### Merging density, canopy height, biomass, and area data 
kelpdat <- merge(densgrp, kelpgrp, area.frame, by = "SiteName", all=TRUE)


#### Plotting out the data ----

tiff(file="./MSc_plots/KelpMetrics3.tiff", height = 5, width = 7, units = "in", res=500)

## Site specific density (summed Macro & Nereo) # Ordered by increasing density
d1 <- ggplot() +
  geom_point(data=dens, size=2, alpha=0.1, aes(x=reorder(SiteName,Kelp), y=Kelp)) +
  geom_pointrange(data=kelpdat, size=0.5, aes(x=reorder(SiteName,KelpM), y=KelpM, ymin=KelpM-KelpSD, ymax=KelpM+KelpSD)) +
  scale_y_continuous(limits=c(-0.8,25), breaks=c(0,6,12,18,24)) +
  theme_classic() +
  theme(
        axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        # axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab("Density (/m2)")
d1

dev.off()

## Site specific canopy height # Ordered by increasing density
c1 <- ggplot() + 
  geom_point(data=kelpbydens, size=2, alpha=0.1, aes(x=reorder(SiteName,KelpM), y=Height_m)) +
  geom_pointrange(data=kelpdat, size=0.5, aes(x=reorder(SiteName,KelpM), y=HeightM, ymin=HeightM-HeightSD, ymax=HeightM+HeightSD)) +
  # scale_y_continuous(limits=c(-0.4,9), breaks=c(0,2,4,6,8)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab("Canopy height (m)") +
  annotate("text", size=4, x=1, y=0.5, label="NA")


# kelpbydens <- kelpbydens[-346,] # Removing outlier point at Second Beach South!!


## Site specific biomass # Ordered by increasing density
b1 <- ggplot() + 
  geom_point(data=kelpbydens, size=2, alpha=0.1, aes(x=reorder(SiteName,KelpM), y=Biomass_g)) +
  geom_pointrange(data=kelpdat, size=0.5, aes(x=reorder(SiteName,KelpM), y=BiomassM, ymin=BiomassM-BiomassSD, ymax=BiomassM+BiomassSD)) +
  # scale_y_continuous(limits=c(-0.4,9), breaks=c(0,2,4,6,8)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab("Biomass (g)") +
  annotate("text", size=4, x=1, y=0.5, label="NA")


## Grouped multi panel plot ##

tiff(file="./MSc_plots/KelpMetrics.tiff", height = 10, width = 7, units = "in", res=500)

ggarrange(d1, c1, b1, ncol=1, align="v", heights=c(1,1,1.4)) # Generating the paneled plot

dev.off() 


#### Regressions of the data ----

kelpdat2 <- merge(kelpdat, rls_richness, by="SiteName", all=T)

ggplot(kelpdat2, aes(x=HeightM, y=KelpM)) +
  geom_point() +
  # geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T) + # Lm quadratic
  # geom_smooth(method="lm", formula = y ~ x, se=T) + # Lm
  # geom_smooth(method="gam", formula = y ~ s(x, bs="cs", k=3), se=T) + # GAM
  theme_classic() +
  xlab("Kelp Density") + ylab("Kelp Biomass")
