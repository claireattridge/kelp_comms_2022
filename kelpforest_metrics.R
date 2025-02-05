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

kelp <- kelp[-266,] # Removing the outlier point at Second Beach South

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


# Averaging to transect with respect to kelp species
# Average of n=5 samples per species
kelptrans <- kelp %>%
  mutate(SiteName = as.factor(SiteName), Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Transect, Species) %>% # Averaging to transect 
  summarise(HeightT = mean(Height_m, na.rm=T),
            BiomassTind_g = mean(Biomass_g, na.rm=T)) # Ave individual biomass (g) / transect


# Bringing in the kelp density data
# Everything exists here at the transect level average
kelpjoin <- merge(kelptrans, dens, by = c("SiteName", "Transect"), all=TRUE)


# i) Converting transect level ave individual kelp biomass from g to kg
# ii) Multiplying each ave individual kelp biomass by the species' density / transect
# Leaves the total biomass / species / transect
kelptog <- kelpjoin %>%
  rowwise() %>% 
  mutate(BiomassTind_kg = (BiomassTind_g/1000)) %>%
  mutate(BiomassTkg = case_when(Species == "Macro" ~ BiomassTind_kg*Macro_5m2,
                                Species == "Nereo" ~ BiomassTind_kg*Nereo_5m2))
kelptog[sapply(kelptog, is.nan)] <- 0 # NaNs to 0s for working with


# Adding the kelp species biomass / transect together 
# Averaging other kelp metrics (height, density) across both species
kelptog_joined <- kelptog %>%
  group_by(SiteName, Transect, Depth_datum_m) %>%
  summarise(BiomassTkg_all = sum(BiomassTkg),
            HeightM = mean(HeightT),
            HeightSD = sd(HeightT),
            DensityM = mean(Kelp),
            DensitySD = sd(Kelp)) 

###############LEAVING OFF HERE REMEMBER TO CARRY THE SD FORWARD FOR AVE OF AVES!!!

# Averaging across transects to the forest level
# Converting to / m2 from the transect / 5m2
kelpgrp <- kelptog_joined %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  group_by(SiteName, Depth_datum_m) %>% 
  summarise(HeightM = weighted.mean(HeightT, na.rm=T), HeightSD = sd(HeightT, na.rm=T), # Ave height (m)
            BiomassM = mean(BiomassTkg, na.rm=T), BiomassSD = sd(BiomassTkg, na.rm=T), # Ave biomass (kg / m2)
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

# I know this is not the most efficient code solution,
# but it does work to produce the areas from each .shp file
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


#### Final data joining ----

### Merging density, canopy height, biomass, and area data 
kelpdat <- merge(kelpgrp, areagrp, by = "SiteName", all=TRUE)

# Filtering out the no kelp sites (Less Dangerous Bay & Wizard I North)
kelpdat <- kelpdat %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North")

# Filtering out the 'redos' of site sampling in the fall
kelpdat <- kelpdat %>%
  filter(SiteName != "Ross Islet 2 REDO") %>%
  filter(SiteName != "Second Beach South REDO") %>%
  filter(SiteName != "Less Dangerous Bay REDO")


# saving a .csv file of the kelp metrics by site
write.csv(kelpdat, "./MSc_data/Data_new/kelp_metrics_2022.csv", row.names=F)

#

#### Plot for kelp forest metrics by sites ----

## Calling the averaged and cleaned kelp data sheet
kelpdat <- read_csv("./MSc_data/Data_new/kelp_metrics_2022.csv")


# Filtering out the no kelp sites (Less Dangerous Bay & Wizard I North)
# from the raw dataframe retaining all ungrouped data points ('kelptog').
# The other raw data points are needed for the following plots, but 
# these 'no-kelp' sites need to be excluded from the datasets
kelptog_clean <- kelptog %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North") %>%
  droplevels()

# Filtering out the 'redos' of site sampling in the fall
kelptog_clean <- kelptog_clean %>%
  filter(SiteName != "Ross Islet 2 REDO") %>%
  filter(SiteName != "Second Beach South REDO") %>%
  filter(SiteName != "Less Dangerous Bay REDO")


# saving a .csv file of the clean raw kelp data points
write.csv(kelptog_clean, "./MSc_data/Data_new/kelp_rawdata_clean_2022.csv", row.names=F)



### SUPPLEMENTAL FIGURE OF KELP FOREST ATTRIBUTES

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
  geom_point(data=kelptog_clean, size=2, alpha=0.2, aes(x=reorder(SiteName,Kelp), y=BiomassTkg)) +
  scale_y_continuous(limits=c(-1,30), breaks=c(0,10,20,30)) +
  theme_classic() +
  theme(
        # axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.89),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size="10"),
        axis.title.y = element_text(color="black", size="12", vjust=1)) +
  xlab("") + ylab(expression("Biomass (kg /m"^2*")"))
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

tiff(file="./MSc_plots/SuppFigs/KelpMetricsAll.tiff", height = 12, width = 6, units = "in", res=400)

ggarrange(d1, c1, b1, a1, ncol=1, align="v", heights=c(1,1,1,1.4)) # Generating the paneled plot

dev.off() 



#### Seasonality plot (2022-2023) ----

## Loading the 2022 kelp metrics
kelp2022 <- read_csv("./MSc_data/Data_new/kelp_metrics_2022.csv")

## Loading the 2023 raw data sheets
dens2023 <- read_csv("./MSc_data/Data_new/kelp_density_2023.csv")
kelp2023 <- read_csv("./MSc_data/Data_new/kelp_morphology_2023.csv")
# Has already been processed in 'forest_area_processing_2023.R'
area2023 <- read_csv("./MSc_data/Data_new/forest_areas_2023.csv")




### Cleaning up the kelp density data
dens2023_select <- dens2023 %>%
  select(SiteName, Transect, Depth_datum_m, Macro_5m2, Nereo_5m2)

## New columns for Macro & Nereo density / m2
dens2023_m2 <- dens2023_select %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>%
  rowwise() %>% # To sum across rows
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup() # Stop rowwise 

## Grouping to site-specific density averages 
dens2023_group <- dens2023_m2 %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(Depth_datum_m = mean(Depth_datum_m), KelpM = mean(Kelp), KelpSD = sd(Kelp))





### Cleaning up the kelp morphology data

## Converting sub-bulb circumference to diameter for Nereo
kelp_subbulb <- kelp2023 %>%
  rowwise() %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159)) %>%
  ungroup()

# Equation for converting from sub-bulb diamater (cm) to biomass (g) as per our Nereo sub-bulb model (using equation for all 3 sample sites combined)
# See 'nereo_biomass.R' for the origin relationship
formula <- function(x){
  (150.7966*(x)^2 -216.2721*(x) + 315.0124)
}

options(scipen=999) # Turning off scientific notation

## Statement for applying sub-bulb equation when biomass is absent (i.e., when it needs to happen)
kelp_subbulb$Biomass_g <- ifelse(is.na(kelp_subbulb$Biomass_g), 
                                 ifelse(is.na(kelp_subbulb$Sub_diam_cm), NA, 
                                        formula(kelp_subbulb$Sub_diam_cm)), kelp_subbulb$Biomass_g)


## Averaging to transect from individual samples (i.e., ave individual biomass / transect)
kelp_transect_biomass <- kelp_subbulb %>%
  mutate(SiteName = as.factor(SiteName), Height_m = as.numeric(Height_m)) %>%
  group_by(SiteName, Transect) %>% # Averaging to transect 
  summarise(Depth_datum_m = mean(Depth_datum_m, na.rm=TRUE),
            HeightT = mean(Height_m, na.rm=TRUE),
            BiomassTind = mean(Biomass_g, na.rm=TRUE)) # Ave individual biomass (g) / transect



### Adding together the kelp morphology and density data

## Removing the depth column from density to merge with the 
## kelp morphology data 
dens2023_m2 <- dens2023_m2 %>%
  select(-Depth_datum_m)

## Bringing in the kelp density (transect level) data
kelptog <- merge(kelp_transect_biomass, dens2023_m2, by = c("SiteName", "Transect"), all=TRUE)


## Using kelp densities to estimate biomass at transect level
kelptog_cleaned <- kelptog %>%
  rowwise() %>% # To sum across rows
  mutate(BiomassTkg = ((BiomassTind/1000)*Kelp)) %>% # Convert ave ind biomass (g to kg) and multiply by transect density
  ungroup() %>% # Stop rowwise 
  mutate(Biomassm2kg = (BiomassTkg/5)) # From ave transect area biomass (/5m2) to /m2 area
kelptog_cleaned[sapply(kelptog_cleaned, is.nan)] <- NA # NaNs to NAs for working with


# Grouping/averaging from transect to site level
kelptog_grouped <- kelptog_cleaned %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(Depth_datum_m = mean (Depth_datum_m, na.rm=T),
            HeightM = mean(HeightT, na.rm=T), HeightSD = sd(HeightT, na.rm=T), # Ave height (m)
            BiomassM = mean(Biomassm2kg, na.rm=T), BiomassSD = sd(Biomassm2kg, na.rm=T), # Ave biomass (kg / m2)
            DensityM = mean(Kelp), DensitySD = sd(Kelp, na.rm=T), # Ave density any kelp (m2)
            MacroM = mean(Macro), MacroSD = sd(Macro, na.rm=T), # Ave density Macro (m2)
            NereoM = mean(Nereo), NereoSD = sd(Nereo, na.rm=T)) # Ave density Nereo (m2)
kelptog_grouped[sapply(kelptog_grouped, is.nan)] <- NA # NaNs to NAs for working with

# Adding composition identity column for each site
kelptog_grouped <- kelptog_grouped %>%
  mutate(Composition = case_when(NereoM != 0 & MacroM != 0 ~ "Mixed",
                                 MacroM != 0 ~ "Macro",
                                 NereoM != 0 ~ "Nereo",
                                 TRUE ~ as.character("None")))


### Joining with the kelp area data

## Bringing in the kelp density (transect level) data
kelptog_grouped_all <- merge(kelptog_grouped, area2023, by = c("SiteName"), all=TRUE)

## saving a .csv file of the kelp metrics by site
write.csv(kelptog_grouped_all, "./MSc_data/Data_new/kelp_metrics_2023.csv", row.names=F)




### MAKING THE KELP SEASONALITY PLOT (2022 vs. 2023)

## Filtering out sites without all metrics examined
## Adding column to identify year
kelpdat2023 <- read_csv("./MSc_data/Data_new/kelp_metrics_2023.csv") %>%
  filter(SiteName != "Less Dangerous Bay" &
           SiteName != "Wizard Islet North" &
           SiteName != "Flemming 112" &
           SiteName != "Ross Islet Slug Island") %>%
  mutate(Year = "2023") %>%
  mutate(SiteName = as.factor(SiteName), Composition = as.factor(Composition))
  

# extracting the site names from those surveyed in 2023
sites2023 <- levels(as.factor(kelpdat2023$SiteName))


## Adding in the 2022 kelp metric data
kelpdat2022 <- read_csv("./MSc_data/Data_new/kelp_metrics_2022.csv") %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  filter(SiteName %in% sites2023) %>%
  droplevels() %>%
  mutate(Year = as.factor("2022"))

# extracting the site names from those surveyed in both 2022 & 2023
sites_shared <- levels(as.factor(kelpdat2022$SiteName))

## Back-filtering the 2023 frame again for only the shared sites with 2022
kelpdat2023 <- kelpdat2023 %>%
  filter(SiteName %in% sites_shared) %>%
  droplevels()


## joining the 2022 & 2023 kelp data together
## Site averages
kelpdat_biannual <- rbind(kelpdat2022, kelpdat2023)




# Kelp density plot
densplot2023 <- ggplot() +
  geom_pointrange(data=kelpdat_biannual, size=0.3, 
                  aes(x=SiteName, y=DensityM, 
                      ymin=(DensityM-DensitySD), ymax=(DensityM+DensitySD),
                      group=Year, shape=Year), position=position_dodge(width=0.4)) +
  # geom_point(data=dens2023_m2, size=2, alpha=0.2, aes(x=SiteName, y=Kelp, group=Year)) +
  scale_y_continuous(breaks=c(0,5,10,15,20, 25), limits=c(-2,25)) +
  scale_shape_manual(values=c(15,16)) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
    # axis.text.x = element_text(color="white", size="9.5", angle=45, vjust=0.89, hjust=0.89),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size="8.5"),
    axis.title.y = element_text(color="black", size="10", vjust=1)) +
  xlab("") + ylab(expression("Density (stipes /m"^2*")"))


# Forest biomass plot
biomassplot2023 <- ggplot() +
  geom_pointrange(data=kelpdat_biannual, size=0.3, 
                  aes(x=SiteName, y=BiomassM, 
                      ymin=(BiomassM-BiomassSD), ymax=(BiomassM+BiomassSD),
                      group=Year, shape=Year), position=position_dodge(width=0.4)) +
  # geom_point(data=kelptog_clean, size=2, alpha=0.2, aes(x=SiteName, y=Kelp)) +
  # scale_y_continuous(limits=c(-2,25)) +
  scale_shape_manual(values=c(15,16)) +
  theme_classic() +
  theme(legend.position = "none",
    # axis.text.x = element_text(color="white", size="9.5", angle=45, vjust=0.89, hjust=0.89),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size="8.5"),
    axis.title.y = element_text(color="black", size="10", vjust=1)) +
  xlab("") + ylab(expression("Biomass (kg wet weight /m"^2*")"))


# Forest biomass plot
areaplot2023 <- ggplot() +
  geom_point(data=kelpdat_biannual, size=1.7, 
                  aes(x=SiteName, y=Area_m2, 
                      group=Year, shape=Year), position=position_dodge(width=0.4)) +
  # geom_point(data=kelptog_clean, size=2, alpha=0.2, aes(x=SiteName, y=Kelp)) +
  scale_y_continuous(breaks=c(0,3000,6000,9000,12000)) +
  scale_shape_manual(values=c(15,16)) +
  theme_classic() +
  theme(legend.position = "none",
    axis.text.x = element_text(color="black", size="8.5", angle=45, vjust=0.89, hjust=0.89),
    axis.text.y = element_text(color="black", size="8.5"),
    axis.title.y = element_text(color="black", size="10", vjust=1)) +
  xlab("") + ylab(expression("Forest area (/m"^2*")"))


## Grouped multi-panel plot for 2022 vs. 2023 kelp forest metrics

tiff(file="./MSc_plots/SuppFigs/Kelp_2022_2023.tiff", height = 7.5, width = 5, units = "in", res=400)

ggarrange(densplot2023, biomassplot2023, areaplot2023, 
          ncol=1, align="v", heights=c(1.2,1,1.4)) 

dev.off() 
