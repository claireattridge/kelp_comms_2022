## Barkley Sound kelp forest data processing
## Author: Claire Attridge
## Origin date: August 2022

# Loading base packages
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(MetBrewer)

#### Kelp density, height, & biomass cleaning ----

####### Reading in our kelp density data
dens <- read.csv("./MSc_data/Data_new/kelp_density_2022.csv") %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>% # Changing units to /m2 area
  rowwise() %>% # To sum across rows
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup() # Stop rowwise 

# Removing the site 'redos' for current purposes
dens <- dens %>%
  filter(!str_detect(SiteName, 'REDO'))

# Dataframe for site-specific density averages 
densgrp <- dens %>%
  group_by(SiteName, Depth_datum_m) %>% # Averaging to site
  summarise(KelpM = mean(Kelp), KelpSD = sd(Kelp))


####### Reading in our kelp height and biomass data
kelp <- read.csv("./MSc_data/Data_new/kelp_morphology_2022.csv") %>%
  as.data.frame()

kelp <- kelp[-266,] # Removing the outlier point at Second Beach South!

# Removing the site 'redos' for current purposes
kelp <- kelp %>%
  filter(!str_detect(SiteName, 'REDO'))


# Converting sub-bulb circumference to diameter
kelp <- kelp %>%
  rowwise() %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159)) %>%
  ungroup()
  
# Equation for converting from sub-bulb diamater (cm) to biomass (g) as per our Nereo sub-bulb model (using equation for all 3 sample sites combined)
# See 'nereo_biomass.R' for the origin relationship
formula <- function(x){
  (150.7966*(x)^2 -216.2721*(x) + 315.0124)
}

options(scipen=999) # Turning off scientific notation

# Statement for applying sub-bulb equation when biomass is absent (i.e., when it needs to happen)
kelp$Biomass_g <- ifelse(is.na(kelp$Biomass_g), ifelse(is.na(kelp$Sub_diam_cm), NA, formula(kelp$Sub_diam_cm)), kelp$Biomass_g)


# Averaging to transect from individual samples (i.e., ave individual biomass / transect)
kelptrans <- kelp %>%
  mutate(SiteName = as.factor(SiteName), Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Transect) %>% # Averaging to transect 
  summarise(HeightT = mean(Height_m, na.rm=T),
            BiomassTind = mean(Biomass_g, na.rm=T)) # Ave individual biomass (g) / transect


# Bringing in the kelp density (transect level) data
kelptog <- merge(kelptrans, dens, by = c("SiteName", "Transect"), all=TRUE)


kelptog <- kelptog %>%
  rowwise() %>% # To sum across rows
  mutate(BiomassTkg = ((BiomassTind/1000)*Kelp)) %>% # Convert ave ind biomass (g to kg) and multiply by transect density
  ungroup() %>% # Stop rowwise 
  mutate(Biomassm2kg = (BiomassTkg/5)) # From ave transect area biomass (/5m2) to /m2 area
kelptog[sapply(kelptog, is.nan)] <- 0 # NaNs to 0s for working with


# Grouping/averaging from transect to site level
kelpgrp <- kelptog %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  group_by(SiteName, Depth_datum_m) %>% # Averaging to site & keeping depth (m)
  summarise(HeightM = mean(HeightT, na.rm=T), HeightSD = sd(HeightT, na.rm=T), # Ave height (m)
            BiomassM = mean(Biomassm2kg, na.rm=T), BiomassSD = sd(Biomassm2kg, na.rm=T), # Ave biomass (kg / m2)
            DensityM = mean(Kelp), DensitySD = sd(Kelp, na.rm=T), # Ave density any kelp (m2)
            MacroM = mean(Macro), MacroSD = sd(Macro, na.rm=T), # Ave density Macro (m2)
            NereoM = mean(Nereo), NereoSD = sd(Nereo, na.rm=T)) # Ave density Nereo (m2)


# Adding composition identity column for each site
kelpgrp <- kelpgrp %>%
  mutate(Composition = case_when(NereoM != 0 & MacroM != 0 ~ "Mixed",
                                 MacroM != 0 ~ "Macro",
                                 NereoM != 0 ~ "Nereo",
                                 TRUE ~ as.character("None")))

#### Kelp area cleaning ----

# Loading packages for working with sf and mapping objects
library(scales)
library(sf)
library(ggsn)
library(MetBrewer)
library(maptools)
library(rgeos)

######## Prepping the backing map (in case you want to plot out any of the .shp files)

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
  st_sfc(crs=latlong) %>% # Origin crs as WGS84
  st_sf() %>%
  st_transform(proj) # Projecting into working crs (BC/ALBERS)
plot(corners)


######## Loading in the kelp forest shape files

# Get all files with the .shp extension from the correct folder (i.e., the QGIS cleaned versions)
path <- "./MSc_data/Data_new/Kelp_area/CleanGPX" # set path to the folder

n <- list.files(path, pattern = "*shp", full.names=FALSE) # Pull only the relevant names

shps <- list.files(path, pattern = "*shp", full.names=TRUE) # Pull full file location names
lshps <- as.list(shps) # Turn character vector into a list

names(lshps) <- sub(" area.shp", "", n) # Assign relevant names (cleaned) to elements of the list

f <- function(x){
  st_read(x, crs=latlong) # Base crs of the GPS used to create the shp files is WGS84 
}

allshps <- lapply(lshps, f) # Apply function to read in shape files as sf objects with crs WGS84

f2 <- function(x){
  st_transform(x, crs=proj) # TO reproject to project crs with units in metres (BC/Albers) 
}

allshps <- lapply(allshps, f2) # Apply function to translate shape files to BC/Albers


list2env(allshps, envir=.GlobalEnv) # Bring each unique sf object to the working environment for plotting             


####### Calculating area of each kelp forest shape file

# I know this is terrible, working on code to automate/loop st_area() through the .shp files
# But, at least it works! Keep an eye out for updates to this code. 
a <- cbind(st_area(allshps[[1]]), as.character(names(allshps[1])))
b <- cbind(st_area(allshps[[2]]), as.character(names(allshps[2])))
c <- cbind(st_area(allshps[[3]]), as.character(names(allshps[3])))
d <- cbind(st_area(allshps[[4]]), as.character(names(allshps[4])))
e <- cbind(st_area(allshps[[5]]), as.character(names(allshps[5])))
f <- cbind(st_area(allshps[[6]]), as.character(names(allshps[6])))
g <- cbind(st_area(allshps[[7]]), as.character(names(allshps[7])))
h <- cbind(st_area(allshps[[8]]), as.character(names(allshps[8])))
i <- cbind(st_area(allshps[[9]]), as.character(names(allshps[9])))
j <- cbind(st_area(allshps[[10]]), as.character(names(allshps[10])))
k <- cbind(st_area(allshps[[11]]), as.character(names(allshps[11])))
l <- cbind(st_area(allshps[[12]]), as.character(names(allshps[12])))
m <- cbind(st_area(allshps[[13]]), as.character(names(allshps[13])))
n <- cbind(st_area(allshps[[14]]), as.character(names(allshps[14])))
o <- cbind(st_area(allshps[[15]]), as.character(names(allshps[15])))
p <- cbind(st_area(allshps[[16]]), as.character(names(allshps[16])))
q <- cbind(st_area(allshps[[17]]), as.character(names(allshps[17])))
r <- cbind(st_area(allshps[[18]]), as.character(names(allshps[18])))
s <- cbind(st_area(allshps[[19]]), as.character(names(allshps[19])))
t <- cbind(st_area(allshps[[20]]), as.character(names(allshps[20])))
u <- cbind(st_area(allshps[[21]]), as.character(names(allshps[21])))
v <- cbind(st_area(allshps[[22]]), as.character(names(allshps[22])))
w <- cbind(st_area(allshps[[23]]), as.character(names(allshps[23])))


# Generating data frame of kelp forest areas
areagrp <- data.frame(rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w)) %>%
  mutate(Area_m2 = as.numeric(X1), SiteName = as.factor(X2)) %>%
  dplyr::select(-c(X1, X2))


# Removing the site 'redos' for current purposes
areagrp <- areagrp %>%
  filter(!str_detect(SiteName, 'REDO'))


# ############### TRYING TO FIGURE OUT HOW TO AUTOMATE THIS OVER THE LIST GAHHHHHHH
# 
# f2 <- function(polygon){
#   st_area(polygon) #
# }
# 
# 
# test <- vector('list', length(allshps))
# for (i in seq_along(allshps)) {
#   x <- st_area(allshps[[i]])
#   test[i] <- x
# }
# testdf <- do.call(rbind, test)
# 
# 
# # Empty List for Centroids
# test <- list()
# 
# # For Loop
# for (i in seq_along(allshps)) {
#     sf::st_area(allshps[[i]])
# }

############### Anyways, moving on for now...

#### Final data joining ----

### Merging density, canopy height, biomass, and area data 
kelpdat <- merge(kelpgrp, areagrp, by = "SiteName", all=TRUE)

# Filtering out the outlier site (Second Beach South) and no kelp sites (Less Dangerous Bay & Wizard I North)
kelpdat <- kelpdat %>%
  filter(SiteName != "Second Beach South") %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North")


# saving a .csv file of the kelp metrics by site
write.csv(kelpdat, "./MSc_data/Data_new/kelp_metrics_2022.csv", row.names=F)

#

#### Plotting out the data ----

# Filtering out the outlier site (Second Beach South) and no kelp sites (Less Dangerous Bay & Wizard I North)
# from the dataframe retaining all ungrouped data points ('kelptog')
# the raw data points are needed for the following plots, but these sites need to be excluded
kelptog_clean <- kelptog %>%
  filter(SiteName != "Second Beach South") %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North") %>%
  droplevels()


## Site specific density (summed Macro & Nereo) # Ordered by increasing density

d1 <- ggplot() +
  geom_pointrange(data=kelpdat, size=0.5, aes(x=reorder(SiteName,DensityM), y=DensityM, ymin=DensityM-DensitySD, ymax=DensityM+DensitySD)) +
  geom_point(data=kelptog_clean, size=2, alpha=0.2, aes(x=SiteName, y=Kelp)) +
  scale_y_continuous(limits=c(-0.8,25), breaks=c(0,6,12,18,24)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.89),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab(expression("Density (stipes /m"^2*")"))
d1



## Site specific canopy height # Ordered by increasing density

c1 <- ggplot() + 
  geom_pointrange(data=kelpdat, size=0.5, aes(x=reorder(SiteName,DensityM), y=HeightM, ymin=HeightM-HeightSD, ymax=HeightM+HeightSD)) +
  geom_point(data=kelptog_clean, size=2, alpha=0.2, aes(x=SiteName, y=HeightT)) +
  scale_y_continuous(limits=c(-0.4,7), breaks=c(0,2,4,6)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.89),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab("Canopy height (m)") 
c1


## Site specific biomass # Ordered by increasing density

b1 <- ggplot() + 
  geom_pointrange(data=kelpdat, size=0.5, aes(x=reorder(SiteName,DensityM), y=BiomassM, ymin=BiomassM-BiomassSD, ymax=BiomassM+BiomassSD)) +
  geom_point(data=kelptog_clean, size=2, alpha=0.2, aes(x=reorder(SiteName,Kelp), y=Biomassm2kg)) +
  scale_y_continuous(limits=c(-0.4,5.5), breaks=c(0,1.5,3,4.5)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.89),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab(expression("Biomass (kg wet weight /m"^2*")"))
b1


## Site specific area # Ordered by increasing density

a1 <- ggplot() + 
  geom_point(data=kelpdat, size=2, aes(x=reorder(SiteName,DensityM), y=Area_m2)) +
  theme_classic() +
  scale_y_continuous(limits=c(-500,15000)) +
  theme(axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.89),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab(expression("Forest area (/m"^2*")")) 
a1



## Grouped multi-panel plot ##

tiff(file="./MSc_plots/KelpMetricsAll.tiff", height = 12, width = 6, units = "in", res=400)

ggarrange(d1, c1, b1, a1, ncol=1, align="v", heights=c(1,1,1,1.4)) # Generating the paneled plot

dev.off() 


#### Silly plot for a presentation ----
# 
# # fxn to generate sequential exponential values
# future_value = function(years, x = 1, increase = 0.2) {
#   x * (1 + increase) ^ (1:length(years))
# }
# # test line
# future_value(1:23)
# 
# # adding fake exponential column
# kelpsilly <- kelpdat %>%
#   arrange(DensityM, desc=T) %>%
#   mutate(Frustration = future_value(1:23))
# 
# tiff(file="C:/Users/Claire/Desktop/sillykelp.tiff", height = 5, width = 7, units = "in", res=600)
# 
# # plotting out exponential curve to density
# silly <- ggplot(data=kelpsilly, size=2, aes(x=as.numeric(DensityM), y=as.numeric(Frustration))) +
#   geom_point(data=kelpsilly, size=2, aes(x=as.numeric(DensityM), y=as.numeric(Frustration))) +
#   geom_smooth(method="lm", formula = (y ~ x + I(x^2)), se=T) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(color="black", size="9.5"),
#     axis.text.y = element_text(color="black", size="10"),
#     axis.title.y = element_text(color="black", size="12", vjust=1)) +
#   xlab(expression("Kelp density (stipes /m"^2*")")) + ylab("Frustration (swears /min)")
# silly
# 
# dev.off()

#### Test regressions of the data ----

# Pulling from the 'RLS_species_code.R' sheet
kelpdat2 <- merge(kelpdat, rls_richness, by="SiteName", all=T)

ggplot(kelpdat2, aes(x=DensityM, y=species_richness)) +
  geom_point() +
  # geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T) + # Lm quadratic
  # geom_smooth(method="lm", formula = y ~ x, se=T) + # Lm
  geom_smooth(method="gam", formula = y ~ s(x, bs="cs", k=3), se=T) + # GAM
  theme_classic() +
  xlab("Kelp Density") + ylab("Spp richness")
