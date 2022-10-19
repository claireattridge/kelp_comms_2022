## Barkley Sound site mapping
## Author: Claire Attridge
## Origin date: May 2022 

library(tidyverse)
library(sf)
library(ggsn)
library(cowplot)
library(raster)
library(rgeos)
library(MetBrewer)

#### Map prep work ####

## setting projections at outset
proj <- st_crs(3005) # ESPG code for BC Albers projection
latlong <- "+proj=longlat +datum=WGS84 +no_defs" # identifying the coordinate system being used for the corner values

## seting map extent for Bamfield 
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

## reading in BC sea map (eez)
sea <- sf::st_read("./Maps_BC/eez_bc/eez.shp") %>%
  st_sf() %>%
  st_transform(proj)
sea <- sea %>% # cropping to the Bamfield corners
  st_crop(st_bbox(corners))

## converting from sea area to land area (because this map is the water not the land)
mask <- st_bbox(sea) %>% # generating a box that covers the cropped multipolygon area
  st_as_sfc() # making a simple geometry feature

sea <- st_union(sea) # unionizing the sea layer

land <- st_difference(mask, sea) # removing the intersecting area (i.e., removing sea and leaving land)
land <- land %>%
  st_sf() %>% # making into sf object 
  st_transform(proj) # projecting back into ESPG 3005


#### Data prep work ####

##                             ##
## Joanne Lessard survey sites ##
##                             ##

library(foreign) # Package for reading in .dbf files

#reading in site location data
jl <- read.dbf("./MSc_data/Data_sourced/Lessard_sites/2021SITE.DBF")

#converting location data to sf object
jlsf <- jl %>%
  dplyr::select(c(TRANSECT, LATEND, LONEND)) %>%
  st_as_sf(coords = c(3,2))

st_crs(jlsf) = latlong
jlsf <- st_transform(jlsf,proj) #make sure you assign this, otherwise it does not save the as the transformed version

#cropping to only points within map limits
jlsf <- jlsf %>%
  st_crop(st_bbox(corners))


##                          ##
## CN/SS & JMS survey sites ##
##                          ## 

rec <- read.csv("./MSc_data/Data_sourced/Fieldedits_Kelp_site_suggestions_2022.csv") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon))

#converting location data to sf object
recsf <- rec %>%
  drop_na(Lat) %>%
  dplyr::select(c(SiteName, Lat, Lon)) %>%
  st_as_sf(coords = c(3,2))

st_crs(recsf) = latlong
recsf <- st_transform(recsf,proj) #make sure you assign this, otherwise it does not save the as the transformed version

#cropping to only points within map limits
recsf <- recsf %>%
  st_crop(st_bbox(corners))


##                                ##
## Our finalized survey site list ##
##                                ##

sites <- read.csv("./MSc_data/Data_new/Final_site_list_2022.csv") %>%
  mutate_all(na_if, "") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon), Estimated_dens=factor(Estimated_dens, levels=c("High", "Med", "Low", "None"))) %>%
  mutate(Deploy = as.factor(ifelse(!is.na(Temp_logger), "YES", "NO")))

#converting location data to sf object
sitessf <- sites %>%
  drop_na(Lat) %>%
  dplyr::select(c(SiteName, Lat, Lon, Estimated_dens, Composition, Deploy)) %>%
  st_as_sf(coords = c(3,2))

st_crs(sitessf) = latlong
sitessf <- st_transform(sitessf,proj) #make sure you assign this, otherwise it does not save the as the transformed version

#cropping to only points within map limits
sitessf <- sitessf %>%
  st_crop(st_bbox(corners))



##                       ##
## Our kelp density data ##
##                       ##

# Will add this as a value column to the site mapping data above
dens <- read.csv("./MSc_data/Data_new/kelp_density_2022.csv") %>%
  mutate(Macro = (Macro_5m2/5), Nereo=(Nereo_5m2/5)) %>% # Changing units to /m2 area
  rowwise() %>%
  mutate(Kelp = sum(Macro,Nereo)) %>% # Sum macro and nereo to get total kelp dens / transect
  ungroup()

densgrp <- dens %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  group_by(SiteName) %>% # Averaging to site
  summarise(KelpM = mean(Kelp), KelpSD = sd(Kelp))

densgrp <- as.data.frame(densgrp)


# Adding the additional kelp data metrics to the mapping data
plotdata <- merge(sitessf, densgrp, by = "SiteName", all = TRUE) %>% #Adding on the kelp density data
            merge(kelpgrp) # Adding in the canopy height and biomass data from 'kelpforest_metrics.R'
# Should leave just the 23 sites that we have all kelp data for


#### Making the maps ####

tiff(file="./MSc_plots/MapDens.tiff", height = 5.5, width = 8.5, units = "in", res=600)

## plotting the kelp composition map (discrete fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey53", color = NA) + #color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=4, shape=21, color="black", aes(fill=Composition)) + 
  geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) +
  scale_fill_manual(values=met.brewer("Degas", 4), name="Canopy composition") +
  theme(panel.grid.major = element_line(colour = "transparent"), #hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.title = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=11),
        panel.border = element_rect(colour="black", fill=NA, size=1)) + # Adding in outer border
  north(land, symbol=12, location="topleft") +
  ggsn::scalebar(land,
               location = "bottomright",
               transform = F,
               st.bottom = F,
               st.size = 3.5,
               height = 0.01,
               dist = 1,
               dist_unit = "km",
               model = 'NAD83')
data_map

dev.off()


## plotting the kelp density map (continuous fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey53", color = NA) + #color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=4, shape=21, color="black", aes(fill=KelpM)) +
  # geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3"), name="Kelp density (/m2)", limits=c(0,16)) +
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
                 model = 'NAD83')
data_map

