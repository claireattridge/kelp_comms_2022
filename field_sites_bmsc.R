## Barkley Sound site mapping
## Author: Claire Attridge
## Origin date: May 2022 

# Loading base packages
library(tidyverse)
library(sf)
library(ggsn)
library(cowplot)
library(raster)
library(rgeos)
library(MetBrewer)

#### Map prep work ----

## setting projections at outset
proj <- st_crs(3005) # ESPG code for BC Albers projection
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")

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


# ## saving the cropped size map of land *Should be able to use this going forward now!
# write_sf(land, "./Maps_BC/eez_bc/land_crop_BarkleyS.shp", overwrite=T)


#### Data prep work ----

# ##                             ##
# ## Joanne Lessard survey sites ##
# ##                             ##
# 
# library(foreign) # Package for reading in .dbf files
# 
# #reading in site location data
# jl <- read.dbf("./MSc_data/Data_sourced/Lessard_sites/2021SITE.DBF")
# 
# #converting location data to sf object
# jlsf <- jl %>%
#   dplyr::select(c(TRANSECT, LATEND, LONEND)) %>%
#   st_as_sf(coords = c(3,2))
# 
# st_crs(jlsf) = latlong
# jlsf <- st_transform(jlsf,proj) #make sure you assign this, otherwise it does not save the as the transformed version
# 
# #cropping to only points within map limits
# jlsf <- jlsf %>%
#   st_crop(st_bbox(corners))


# ##                                         ##
# ## Chris Neufeld & Sam Starko survey sites ##
# ##                                         ## 
# 
# rec <- read.csv("./MSc_data/Data_sourced/Fieldedits_Kelp_site_suggestions_2022.csv") %>%
#   mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon))
# 
# #converting location data to sf object
# recsf <- rec %>%
#   drop_na(Lat) %>%
#   dplyr::select(c(SiteName, Lat, Lon)) %>%
#   st_as_sf(coords = c(3,2))
# 
# st_crs(recsf) = latlong
# recsf <- st_transform(recsf,proj) #make sure you assign this, otherwise it does not save the as the transformed version
# 
# #cropping to only points within map limits
# recsf <- recsf %>%
#   st_crop(st_bbox(corners))


##                            ##
## Our finalized survey sites ##
##                            ##

sites <- read.csv("./MSc_data/Data_new/Final_site_list_2022.csv") %>%
  mutate_all(na_if, "") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon), Estimated_dens=factor(Estimated_dens, levels=c("High", "Med", "Low", "None"))) %>%
  mutate(Deploy = as.factor(ifelse(!is.na(Temp_logger), "YES", "NO")))

# Adding on a column to identify the cluster groups (see 'comm_analyses.R' for details)
sites <- sites %>%
mutate(SiteName = as.character(SiteName)) %>% 
  mutate(Cluster = case_when(SiteName == "Between Scotts and Bradys" | SiteName == "Danvers Danger Rock" | SiteName == "Dixon Island Back (Bay)" | SiteName == "Dodger Channel 1" | SiteName == "Flemming 112" | SiteName == "Flemming 114" | SiteName == "North Helby Rock" | SiteName == "Tzartus 116" ~ "C3",
                             SiteName == "Dodger Channel 2" | SiteName == "Nanat Bay" ~ "C2",
                             SiteName == "Ed King East Inside" | SiteName == "Ross Islet 2" | SiteName == "Ross Islet Slug Island" | SiteName == "Taylor Rock" | SiteName == "Turf Island 2" | SiteName == "Wizard Islet South" ~ "C4",
                             SiteName == "Bordelais Island" | SiteName == "Second Beach" ~ "C5",
                             SiteName == "Cable Beach (Blow Hole)" | SiteName == "Less Dangerous Bay" | SiteName == "Swiss Boy" | SiteName == "Wizard Islet North" ~ "C1",
                             TRUE ~ "Aux"))


#converting location data to sf object
sitessf <- sites %>%
  drop_na(Lat) %>%
  dplyr::select(c(SiteName, Lat, Lon, Estimated_dens, Composition, Deploy, Cluster)) %>%
  st_as_sf(coords = c(3,2))

st_crs(sitessf) <- latlong # setting 
sitessf <- st_transform(sitessf,proj) # make sure you assign this, otherwise it won't use the transformed project crs

#cropping to only points within map limits
sitessf <- sitessf %>%
  st_crop(st_bbox(corners))



##                       ##
## Our kelp metric data  ##
##                       ##

# Adding the kelp data metrics to the mapping data!
# 'kelpforest_metrics.R' (load the objects in this sheet first)
plotdata <- merge(sitessf, kelpdat, by = "SiteName", all = TRUE) %>% # Adding on the kelp data
  mutate(Cluster = factor(Cluster, levels=c("C1", "C2", "C3", "C4", "C5", "Aux")))


#### Plotting the maps ----

tiff(file="./MSc_plots/MapForests.tiff", height = 5.5, width = 8.5, units = "in", res=600)

## plotting the kelp composition map (discrete fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey53", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=4, shape=21, color="black", aes(fill=Composition)) + 
  geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) + # provides site names at points
  scale_fill_manual(values=met.brewer("Degas", 4), name="Canopy composition") +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
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



tiff(file="./MSc_plots/MapDens.tiff", height = 5, width = 10, units = "in", res=600)

## plotting the kelp density map (continuous fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey53", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=3, shape=21, color="black", aes(fill=KelpM)) +
  # geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3"), name=expression("Kelp density (stipes /m"^2*")"), limits=c(0,16)) +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(color="black", size=11),
        legend.text = element_text(color="black", size=10),
        panel.border = element_rect(colour="black", fill=NA, size=1)) + # Adding in outer border
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

dev.off()



# List of relevant sites with 2022 data collection 
studysites <- as.vector(c("Between Scotts and Bradys","Bordelais Island","Cable Beach (Blow Hole)","Danvers Danger Rock","Dixon Island Back (Bay)","Dodger Channel 1","Dodger Channel 2","Ed King East Inside","Flemming 112","Flemming 114","Less Dangerous Bay","Nanat Bay","North Helby Rock","Ross Islet 2","Ross Islet Slug Island","Second Beach","Second Beach South","Swiss Boy","Taylor Rock","Turf Island 2","Tzartus 116","Wizard Islet North","Wizard Islet South"))
  
plotdataclust <- plotdata %>%
  filter(SiteName %in% studysites)

cols <- c("#e4632d", "#994455", "#015f60", "#4477aa", "#997700")

tiff(file="./MSc_plots/MapClust.tiff", height = 6, width = 8, units = "in", res=400)

## plotting the kelp cluster groups map (discrete fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey70", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdataclust, size=4, shape=21, color="black", aes(fill=Cluster)) + 
  scale_fill_manual(values=c(cols, "grey45"), name="") +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.position = "top",
        legend.title = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=14, face="bold"),
        panel.border = element_rect(colour="black", fill=NA, linewidth=1), # border
        plot.margin = unit(c(0,0.5,0.5,-4), "cm")) + # shrinking margins
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

#### Paper Fig. 1 ----

# Calling cluster plot as magick image from 'comm_analyses.R'
mclust <- image_read("./MSc_plots/Clusters_kmeans.tiff") 
# mclustscl <- image_scale(mclust, "2000") # Scaling it to smaller size if needed

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/MapClust.tiff") 
# mplotscl <- image_scale(mplot, "2000") # Scaling it to smaller size if needed

# Adding in stack of kelp forest images from 'comm_analyses.R'
mmap <- image_composite(mplot, imgstk, offset="+2565+370")

# Stacking the map & cluster plots together
mfinal <- image_append(c(mmap, mclust), stack=TRUE)

# Saving the final fig
image_write(mfinal, path = "./MSc_plots/PaperFigs/Fig1a.tiff", format = "tiff", density=600)

