## Barkley Sound Electric Blue temperature logger processing
## Author: Claire Attridge
## Origin date: October 2022

# Loading base packages
library(tidyverse)
library(scales)
library(data.table)
library(sf)
library(ggsn)
library(MetBrewer)
library(wesanderson)

#### 2022 temperature records ----


### Loading & cleaning the envlogger data

# Setting pathway to 'Temp_loggers' folder as a variable
mypath <- "./MSc_data/Data_new/Temp_loggers/Temps2022"

# File list with path in name
files <- list.files(path = mypath, pattern="*.csv", full.names = TRUE)
# File list without path in name
files_names <- list.files(path = mypath, pattern="*.csv", full.names = FALSE)

# Reading all temp csvs into one long dataframe
temp_h <- NULL
#if you rerun the loop, make sure to rerun the line above as well
for (i in 1:length(files)){
  tmp_h <- read_csv(files[[i]], col_names = FALSE) %>% 
    slice(-c(1:21)) %>% 
    mutate(SiteName = files_names[[i]]) %>% 
    separate(SiteName, sep = "[.]", into = c("SiteName", NA))
  temp_h <- temp_h %>% 
    bind_rows(tmp_h)
}

# Setting time as POSIXct
temp <- temp_h %>%
  mutate(Date = as.POSIXct(X1, format="%Y-%m-%d %H:%M:%S"), 
         Temp = as.numeric(X2), 
         SiteName = as.factor(SiteName)) %>%
  ungroup() %>%
  dplyr::select(-X1, -X2)
temp <- as.data.frame(temp) # Making into a data frame


# Cleaning up the lead & butt ends based on deployment timelines
tempf <- temp %>%
  filter(ifelse(SiteName == "Ohiat", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-21"),
          ifelse(SiteName == "Tzartus 116", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
            ifelse(SiteName == "Nanat Bay", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-17"),
              ifelse(SiteName == "Between Scotts and Bradys", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-20"),
                ifelse(SiteName == "Bordelais Island", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-20"),
                  ifelse(SiteName == "Ross Islet 2", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-12"),
                    ifelse(SiteName == "Wizard Islet South", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-19"),
                      ifelse(SiteName == "North Helby Rock", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-21"),
                        ifelse(SiteName == "Turf Island 2", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-20"),
                          ifelse(SiteName == "Robbers Passage 2", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                            ifelse(SiteName == "Flemming 114", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-19"),
                              ifelse(SiteName == "Flemming 112", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                                ifelse(SiteName == "Dodger Channel 1", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-13"),
                                  ifelse(SiteName == "Cable Beach (Blow Hole)", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-14"),
                                    ifelse(SiteName == "Ed King East Inside", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-13"),
                                      ifelse(SiteName == "Ross Islet Slug Island", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-20"),
                                        ifelse(SiteName == "Taylor Rock", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-04"),
                                          ifelse(SiteName == "Dodger Channel 2", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-13"),
                                            ifelse(SiteName == "Second Beach South", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-13"),
                                              ifelse(SiteName == "3 Tree Island 2", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-21"),
                                                ifelse(SiteName == "Second Beach", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-14"),
                                                  ifelse(SiteName == "Wizard Islet North", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-19"),
                                                    ifelse(SiteName == "Dixon Island Back (Bay)", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-19"),
                                                      ifelse(SiteName == "Past Roquefoil Bay", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-19"),
                                                        ifelse(SiteName == "Danvers Danger Rock", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-17"),
                                                          ifelse(SiteName == "Swiss Boy", Date >= as.POSIXct("2022-07-22") & Date <= as.POSIXct("2022-09-18"),
                                                            ifelse(SiteName == "Less Dangerous Bay", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-17"),
                                                                    Date == Date))))))))))))))))))))))))))))



### Averaging/grouping the envlogger data & combining with depths

# Filtering out the sites that were not included in the community analyses
# Filtering out the no kelp sites (Less Dangerous Bay & Wizard I North)
tempf <- tempf %>%
  filter(SiteName != "3 Tree Island 2") %>%
  filter(SiteName != "Ohiat") %>%
  filter(SiteName != "Past Roquefoil Bay") %>%
  filter(SiteName != "Robbers Passage 2") %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North") %>%
  droplevels()


# Finding average, max, min by day, month and site
tempgrp_day <- tempf %>%
  group_by(SiteName, lubridate::day(Date)) %>%
  summarise(Tempave_m = mean(Temp, na.rm=T), Tempmax_m = max(Temp, na.rm=T), Tempmin_m = min(Temp, na.rm=T)) %>%
  ungroup() %>%
  as.data.frame()

# Finding site average, max, min
tempgrp <- tempgrp_day %>%
  group_by(SiteName) %>%
  summarise(Tempave = mean(Tempave_m, na.rm=T), SD_Tempave = sd(Tempave_m, na.rm=T), Tempmin = mean(Tempmin_m, na.rm=T), Tempmax = mean(Tempmax_m, na.rm=T)) %>%
  arrange(desc(Tempmax)) %>% 
  ungroup() %>%
  as.data.frame()


tempgrp <- as.data.frame(tempgrp) # Making into a data frame


# saving a .csv file of the site temp groupings
write.csv(tempgrp, "./MSc_data/Data_new/temps_2022.csv", row.names=F)


## Loading temp logger info sheet
info <- read_csv("./MSc_data/Data_new/temp_logger_info.csv") %>%
  dplyr::select(c(ID, SiteName, Depth_logger_m, Depth_logger_datum_m)) %>%
  mutate(SiteName = as.factor(SiteName), ID = as.factor(ID)) %>%
  filter(ID != "L31" & ID != "L32" & ID != "L33") %>% # Removing the site deployment loggers
  as.data.frame()

tempgrp <- merge(tempgrp, info, by= "SiteName", all = FALSE) # Adding info the the logger data frame


# saving a .csv file of the site temp groupings
write.csv(tempgrp, "./MSc_data/Data_new/temps_2022.csv", row.names=F)


#

                        
#### 2023 temperature records ----

### Loading & cleaning the envlogger data (SUMMER RETRIEVALS)

# Setting pathway to 'Temp_loggers' folder as a variable
mypath <- "./MSc_data/Data_new/Temp_loggers/Temps2023/Summer"

# File list with path in name
files <- list.files(path = mypath, pattern="*.csv", full.names = TRUE)
# File list without path in name
files_names <- list.files(path = mypath, pattern="*.csv", full.names = FALSE)

# Reading all temp csvs into one long dataframe
temp_h <- NULL
#if you rerun the loop, make sure to rerun the line above as well
for (i in 1:length(files)){
  tmp_h <- read_csv(files[[i]], col_names = FALSE) %>% 
    slice(-c(1:21)) %>% 
    mutate(SiteName = files_names[[i]]) %>% 
    separate(SiteName, sep = "[.]", into = c("SiteName", NA))
  temp_h <- temp_h %>% 
    bind_rows(tmp_h)
}

# Setting time as POSIXct
temp <- temp_h %>%
  mutate(Date = as.POSIXct(X1, format="%Y-%m-%d %H:%M:%S"), 
         Temp = as.numeric(X2), 
         SiteName = as.factor(SiteName)) %>%
  ungroup() %>%
  dplyr::select(-X1, -X2)
temp <- as.data.frame(temp) # Making into a data frame


# Cleaning up the lead & butt ends based on deployment timelines
tempf <- temp %>%
  filter(ifelse(SiteName == "Ohiat", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-21"),
                ifelse(SiteName == "Tzartus 116", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                       ifelse(SiteName == "Nanat Bay", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-17"),
                              ifelse(SiteName == "Between Scotts and Bradys", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-20"),
                                     ifelse(SiteName == "Bordelais Island", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-20"),
                                            ifelse(SiteName == "Ross Islet 2", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-12"),
                                                   ifelse(SiteName == "Wizard Islet South", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-19"),
                                                          ifelse(SiteName == "North Helby Rock", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-21"),
                                                                 ifelse(SiteName == "Turf Island 2", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-20"),
                                                                        ifelse(SiteName == "Robbers Passage 2", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                                                                               ifelse(SiteName == "Flemming 114", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-19"),
                                                                                      ifelse(SiteName == "Flemming 112", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                                                                                             ifelse(SiteName == "Dodger Channel 1", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-13"),
                                                                                                    ifelse(SiteName == "Cable Beach (Blow Hole)", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-14"),
                                                                                                           ifelse(SiteName == "Ed King East Inside", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-13"),
                                                                                                                  ifelse(SiteName == "Ross Islet Slug Island", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-20"),
                                                                                                                         ifelse(SiteName == "Taylor Rock", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-04"),
                                                                                                                                ifelse(SiteName == "Dodger Channel 2", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-13"),
                                                                                                                                       ifelse(SiteName == "Second Beach South", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-13"),
                                                                                                                                              ifelse(SiteName == "3 Tree Island 2", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-21"),
                                                                                                                                                     ifelse(SiteName == "Second Beach", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-14"),
                                                                                                                                                            ifelse(SiteName == "Wizard Islet North", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-19"),
                                                                                                                                                                   ifelse(SiteName == "Dixon Island Back (Bay)", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-19"),
                                                                                                                                                                          ifelse(SiteName == "Past Roquefoil Bay", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-19"),
                                                                                                                                                                                 ifelse(SiteName == "Danvers Danger Rock", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-17"),
                                                                                                                                                                                        ifelse(SiteName == "Swiss Boy", Date >= as.POSIXct("2022-07-22") & Date <= as.POSIXct("2022-09-18"),
                                                                                                                                                                                               ifelse(SiteName == "Less Dangerous Bay", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-17"),
                                                                                                                                                                                                      Date == Date))))))))))))))))))))))))))))

#### Combined temperature records ----


#### Plots of data ----


# Calling the temps datasheet
tempgrp <- read_csv("./MSc_data/Data_new/temps_2022.csv")

# Calling the site nums datasheet
sitenums <- read_csv("./MSc_data/Data_new/Site_list_numbers.csv")

# Making a named vector of labels for facet
facet_labels <- as.character(sitenums$SiteNum)
names(facet_labels) <- c(seq(1,21))




## Site logger depth series

tiff(file="./MSc_plots/SuppFigs/SiteDepthsLogger.tiff", height = 6, width = 9, units = "in", res=600)

ggplot() +
  geom_point(data=tempgrp, size=3, alpha=0.9, aes(x=reorder(SiteName,Depth_logger_datum_m), y=Depth_logger_datum_m)) +
  theme_classic() +
  scale_y_continuous(breaks=c(1,2,3,4,5,6)) +
  theme(
    axis.text.x = element_text(color="black", size="9", angle=55, vjust=0.89, hjust=0.85),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.y = element_text(color="black", size="11", vjust=3)) +
  ylab("Depth (m) at datum") + xlab("")

dev.off()





# Extending existing color scheme from Met Brewer
colourCount = length(unique(tempf$SiteName))
getPalette = colorRampPalette(met.brewer(name="Redon", n=12))


## All temperature series over time (multi-panel for Supp. figs)

tiff(file="./MSc_plots/SuppFigs/TempSeriesAll.tiff", height = 9, width = 8.5, units = "in", res=500)

ggplot(tempf, aes(x=Date, y=Temp, group=SiteName, color=SiteName)) +
  geom_line(linewidth=0.5) +
  scale_x_datetime(labels = date_format("%b %d"), breaks = date_breaks("weeks")) + 
  xlab("") + ylab("Temperature (°C)") +
  scale_color_manual(values=getPalette(colourCount)) +
  theme_minimal() +
  theme(legend.position="none",
        strip.text.x=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=12),
        axis.text.y=element_text(color="black", size=9),
        axis.text.x=element_text(color="black", size=8, angle=55, vjust=0.89, hjust=0.85),
        panel.grid.major=element_line(color="grey66")) +
  facet_wrap(~SiteName, ncol=3) +
  theme(strip.background = element_rect(fill = "white", color = "white"),
        strip.text.x = element_text(hjust = 0, size = 12, color="black"))



dev.off()



## Temperature against deployment depths
ggplot(tempgrp, aes(x=Depth_logger_datum_m, y=Tempave)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x, se=T) + # Lm
  # scale_x_reverse() +
  theme_classic()


#

#### Interpolated temp maps ----

###               ###
### MAP PREP WORK ###
###               ###

library(gstat)

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



## reading in full size BC sea map (eez) *Should only have to do this the first time
sea <- sf::st_read("./Maps_BC/eez_bc/eez.shp", crs=4326) %>%
  st_sf() 

# projecting to UTM
seap <- sea %>%
  st_transform(proj)

# sea <- st_make_valid(sea) # going over the map to check for 'invalid' points of error (i.e., data cleaning)

seacrop <- seap %>% # cropping to the Bamfield corners
  st_crop(st_bbox(corners))

## saving the cropped size BC sea map (eez) *Should be able to use this going forward now!
write_sf(seacrop, "./Maps_BC/eez_bc/eez_crop_BarkleyS.shp", overwrite=T)




###                ###
### DATA PREP WORK ###
###                ###

# reading in out site data
sites <- read.csv("./MSc_data/Data_new/Final_site_list_2022.csv") %>%
  mutate_all(na_if, "") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon), Estimated_dens=factor(Estimated_dens, levels=c("High", "Med", "Low", "None"))) %>%
  mutate(Deploy = as.factor(ifelse(!is.na(Temp_logger), "YES", "NO")))

# converting location data to sf object
sites <- sites %>%
  drop_na(Lat) %>%
  dplyr::select(c(SiteName, Lat, Lon, Estimated_dens, Composition, Deploy))
sitessf <- st_as_sf(sites, coords = c(3,2), crs=latlong) # set points as original crs (WGS84)

sitessf <- st_transform(sitessf,proj) # now assign the new system to project into (is still the original in this case)

# cropping to only points within map limits
sitessf <- sitessf %>%
  st_crop(st_bbox(corners))


# joining to the temp data
plotdata <- merge(sitessf, tempgrp, by = "SiteName", all = TRUE)
# Removing the NA temp sites
plotdata <- plotdata %>%
  filter(!is.na(Tempave))



###                  ###
### RASTER PREP WORK ###
###                  ###

library(ncdf4)
library(raster)
library(stars)

# loading in cropped version of the eez ocean shape file (Barkley Sound)
seacrop <- sf::st_read("./Maps_BC/eez_bc/eez_crop_BarkleyS.shp") %>%
  st_set_crs(proj)

# loading in the NOAA netcdf file for the Barkley Sound, Canada 1 arc-second Coastal Digital Elevation Model (DEM)
dem <- raster("./MSc_data/Data_sourced/DEM/barkley_sound_1_navd88_2016.nc")

# projecting to the same UTM crs as our sea layer
dem <- projectRaster(dem, crs=crs(seacrop))


# setting spatial extent of the DEM raster to extent of plot domain (rectangular output based on the limits) 
extent(dem) <- extent(seacrop)


# # going to a more course cell resolution (taking too long to process at 1 arc (i.e., 30m))
# dem150 <- aggregate(dem, fact = 5, fun = "mean") # Up to 150 x 150 m


# rasterizing the sea multipolygon to the DEM cell size *Should only have to do this the first time
searas <- rasterize(seacrop, dem, fun="last")

## saving the raster version of BC sea map (eez) *Can use this going forward now!
writeRaster(searas, "./Maps_BC/eez_bc/eez_crop_ras_BarkleyS.tif", format="GTiff", overwrite=T)



# converting the sf multipolygon of sea to a spatial polygon for plotting & limits
seacropsp <- seacrop %>%
  as("Spatial")



###               ###
### INTERPOLATION ###
###               ###

library(geoR)
library(automap)

# loading in the raster version of BC sea map cropped (eez) *saved in code above
searas <- raster("./Maps_BC/eez_bc/eez_crop_ras_BarkleyS.tif")

# # going to a more course cell resolution (taking too long to process when at 1 arc (i.e., 30m))
# searas150 <- aggregate(searas, fact = 5, fun = "mean") # Up to 150 x 150 m


########### Checking for assumptions prior to interpolating

## NORMALITY ##
hist(plotdata$Tempave) # looks ok
hist(plotdata$Tempmax) # looks good
hist(plotdata$Tempmin) # looks ok

## SPATIAL CONTINUITY (Variogram) *See helpful links below

#https://gis.stackexchange.com/questions/287988/kriging-example-with-sf-object
#https://vsp.pnnl.gov/help/Vsample/Kriging_Variogram_Model.htm#:~:text=Nugget%3A,and%20not%20enough%20spatial%20correlation.
#https://www.youtube.com/watch?v=iKTut0JfvRg
  

####### Making the variogram fit

# converting the sf dataframe to spatial points
plotsp <- as_Spatial(plotdata)

# fitting variogram of data points to the spatial map data
vario <- variogram(Tempave~1, data=plotsp)
# eyefit(vario)

# fitting model to the variogram
fit <- fit.variogram(vario, model = vgm(model="Exp")) # fit model
## Did not specify all parameters so that fit.variogram() will optimize
# Model type: Exponential (Exp) and Spherical (Sph) are most widely applied. 
# Nugget: Estimate of non-spatial variance
# PSill: Estimate of total variance (max value of the variogram) minus the nugget
# Range: Distance at which the variogram attains the sill (max) (i.e., Dist beyond which data are no longer correlated)

plot(vario, fit)


####### Interpolating using ordinary Kriging (AUTOKRIGE) 
# Helpful link: https://www.neonscience.org/resources/learning-hub/tutorials/raster-data-r

# Double checking that the land is set to NA values for my raster layer
plot(searas, colNA="red")


## Converting the raster to spatial points for interpolating over
rassp <- as(searas,"SpatialPoints")

## Converting the sf dataframe of my data points to spatial 
plotsp <- as_Spatial(plotdata)


## Converting the raster of dem to spatial 
demfr <- as.data.frame(dem,xy=TRUE)
demsf <- st_as_sf(x = demfr, 
                        coords = c("x", "y"),
                        crs = proj)
demsfclean <- na.omit(demsf)
demsp <- as(demsfclean, "Spatial")


### running the ordinary kriging

k <- autoKrige(Tempave~1, input_data=plotsp, new_data=rassp)

plot(k$krige_output)
plot(k$krige_output[,"var1.pred"])
plot(k$krige_output[,"var1.var"])
plot(k$krige_output[,"var1.stdev"])


### saving the predicted outputs from the ordinary kriging

kras <- rasterFromXYZ(k$krige_output) # building a raster from the spatial pts df 

crs(kras) <- crs(searas) # setting the crs of new raster to the BC/Albers layers

writeRaster(kras, "./MSc_data/Data_new/Maps/AveTempRas.tif", format="GTiff", overwrite=T)

templayer <- raster("./MSc_data/Data_new/Maps/AveTempRas.tif") # Can load this instead going forward now!



############ Interpolating using co-Kriging (STILL NOT FUNCTIONAL - TBA)
# Helpful link: https://gis.stackexchange.com/questions/280850/preparing-dataset-to-perform-co-kriging-in-r-gstat?fbclid=IwAR1-facpiy2c-pZc2F4EkA0bSNCojYcPLxlx2Dt8oi5GjfhTyCUfqd4KaX0
# Intent is to use co-kriging to predict the spatial variation of subtidal temps with respect to influence of depth

g <- gstat(NULL, id="Temps", form=Tempave~1, data=plotsp) # gstat object for temp data
g <- gstat(g, id="bathy", form=GDAL.Band.Number.1~1, data=demsp) # gstat object for dem data
v <- variogram(g, cutoff=2000) # variograms for both & the cross-variogram

g <- gstat(g, id="Temps", model=v, fill.all=T) # adding the variogram models to the gstat object
pred <- predict.gstat(g, rassp) # predicting over the empty grid of my study area



###                   ###
### PLOTTING THE MAPS ###
###                   ###


## kelp forest (logger) temp by point

tiff(file="./MSc_plots/MapTempPt.tiff", height = 6, width = 9, units = "in", res=600)

data_map <- ggplot(land)+
  geom_sf(fill = "grey53", color = NA) + #color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=4, shape=21, color="black", aes(fill=Tempave)) +
  # geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) +
  scale_fill_gradientn(colors=rev(met.brewer("VanGogh2")), name="Ave. Temp (°C)") +
  theme(panel.grid.major = element_line(colour = "transparent"), #hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  north(land, symbol=12, location="topleft") +
  ggsn::scalebar(land,
                 location = "bottomright",
                 transform = F,
                 st.bottom = F,
                 st.size = 4,
                 height = 0.01,
                 dist = 2,
                 dist_unit = "km",
                 model = 'WGS84')
data_map

dev.off()



## kelp forest (logger) temp by kriging prediction gradient

# quick edit of the point data by column for plotting
plotdata <- plotdata %>%
  mutate(Composition=ifelse(is.na(Composition), "None", Composition)) %>% # changing NAs to a factor level for no kelp
  mutate(Composition = as.factor(Composition)) # making into factor column


tiff(file="./MSc_plots/MapTempGrd.tiff", height = 6.5, width = 8, units = "in", res=600)

col <- wes_palette("Zissou1", 20, type="continuous") # setting gradient color palette
image(templayer, zlim=c(12,15), col=col, xlab="", ylab="", axes=F) # plotting sea gradient
box(col="black") # adding outer box lines
plot(seacropsp, add=T) # adding land borders 

fill <- met.brewer(name="Kandinsky", n=4, type="continuous") # setting point color palette
palette(fill)
plot(plotdata, add=T, pch=21, cex=1.5, col="black", bg=plotdata$Composition) # plotting points
legend("topleft", pch=21, legend=c(expression(italic("M. pyrifera")), "Mixed", expression(italic("N. luetkeana")), "None"), title = "", pt.bg = fill, bty="n", cex=1)

dev.off()

