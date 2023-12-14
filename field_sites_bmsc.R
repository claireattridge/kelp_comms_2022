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


#

#### Spatial package prep work ----

## Have to do this section because maptools and gssn were retired/archived
## There currently (2023) isn't a better replacement package for the tools
## Used for map features such as scale bars and north arrows

## Installing map tools & ggsn & rgeos (from archived versions bc they've been retired as of 2023)
## Maptools is a dependency that needs to be loaded prior to ggsn
url1 <- "https://cran.r-project.org/src/contrib/Archive/maptools/maptools_1.1-8.tar.gz"
url2 <- "https://cran.r-project.org/src/contrib/Archive/ggsn/ggsn_0.5.0.tar.gz"
url3 <- "https://cran.r-project.org/src/contrib/Archive/rgeos/rgeos_0.6-4.tar.gz"

# repos=NULL indicates that you're not using a repository for installation
# type=source specifies that it is a source package
install.packages(url1, type="source", repos=NULL)  
install.packages(url2, type="source", repos=NULL)  
install.packages(url3, type="source", repos=NULL)  

#


#### Map prep work ----

## setting projections at outset
proj <- st_crs(3005) # ESPG code for BC Albers projection
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")

# seting map extent for Bamfield
ymax <- 48.922
ymin <- 48.80
xmax <- -125.05
xmin <- -125.26

# # seting map extent for vancouver Island
# ymax <- 48.0
# ymin <- 51.0
# xmax <- -120.00
# xmin <- -129.00

## making corners for the area of interest 
corners <- st_multipoint(matrix(c(xmax,xmin,ymax,ymin),ncol=2)) %>% 
  st_sfc(crs=latlong) %>%
  st_sf() %>%
  st_transform(proj) 
plot(corners)

## reading in high res map of BC coast from Hakai
land <- sf::st_read("./Maps_BC/Hakai_BC/COAST_TEST2.shp") %>%
  st_sf() %>%
  st_transform(proj)
land_bs <- land %>% # cropping to the Bamfield corners
  st_crop(st_bbox(corners))

## converting from land area to sea area (because this map is the land not the water)
mask <- st_bbox(land_bs) %>% # generating a box that covers the cropped multipolygon area
  st_as_sfc() # making a simple geometry feature

land_bs <- st_union(land_bs) # unionizing the land layer

sea_bs <- st_difference(mask, land_bs) # removing the intersecting area (i.e., removing sea and leaving land)
sea_bs <- sea_bs %>%
  st_sf() %>% # making into sf object 
  st_transform(proj) # projecting back into ESPG 3005

## making the land back into an sf object
land_bs <- land_bs %>%
  st_sf() %>% # making into sf object 
  st_transform(proj) # projecting back into ESPG 3005


## saving the cropped size map of land 
## saving the cropped size map of sea
write_sf(land_bs, "./Maps_BC/land_Hakai_BarkleyS.shp", overwrite=T)
write_sf(sea_bs, "./Maps_BC/sea_Hakai_BarkleyS.shp", overwrite=T)


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

# Loading site list and filtering for only those surveyed
sites <- read.csv("./MSc_data/Data_new/Final_site_list_2022.csv") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon), Estimated_dens=factor(Estimated_dens, levels=c("High", "Med", "Low", "None")),
         Kelpdom = as.factor(ifelse(Composition == "Macro" | Composition == "Mixed", "Macro", "Nereo"))) %>%
  filter(Surveyed == "YES")

# Removing the outlier site (Second Beach South) and no kelp sites (Less Dangerous Bay & Wizard I North)
sites <- sites %>%
  filter(SiteName != "Second Beach South") %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North")
# Should leave the 20 study sites remaining


#converting site location data to sf object
sitessf <- sites %>%
  drop_na(Lat) %>%
  dplyr::select(c(SiteName, Lat, Lon, Estimated_dens, Surveyed, Kelpdom)) %>%
  st_as_sf(coords = c(3,2))

# Setting the initial crs, then reprojecting to BC/Albers crs
st_crs(sitessf) <- latlong
sitessf <- st_transform(sitessf,proj) # make sure you assign this, otherwise it won't use the transformed project crs

#cropping to only points within map limits
sitessf <- sitessf %>%
  st_crop(st_bbox(corners))



##                       ##
## Our kelp metric data  ##
##                       ##

# Loading the kelp data file first
kelpdata <- read_csv("./MSc_data/Data_new/kelp_metrics_2022.csv")


# Adding the kelp data to the site location data
# Reversing level order for the Kandinsky colors to fill in correctly (gold to Macro teal to Nereo)
plotdata <- merge(sitessf, kelpdata, by = "SiteName", all = TRUE) %>%
  mutate(Kelpdom = factor(Kelpdom, levels=c("Nereo", "Macro"))) 


#### Plotting the Barkley Sound map ----

tiff(file="./MSc_plots/MapForests.tiff", height = 5.5, width = 8.5, units = "in", res=400)

## plotting the kelp species map (discrete fill)
data_map <- ggplot(land_bs)+
  geom_sf(fill = "grey65", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=4, shape=21, color="black", aes(fill=Kelpdom)) + 
  # geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) + # provides site names at points
  scale_fill_manual(values=met.brewer("Kandinsky", 2), name="") +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text = element_text(color="black", size=12),
        legend.title = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=11),
        plot.margin = unit(c(0.5,0.5,0.2,0), "cm"), # Selecting margins
        panel.border = element_rect(colour="black", fill=NA, linewidth=1)) + # Adding in outer border
  north(land, symbol=12, location="topleft") +
  ggsn::scalebar(land_bs,
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



tiff(file="./MSc_plots/MapDepth.tiff", height = 7, width = 10, units = "in", res=600)

## plotting the kelp density map (continuous fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey53", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, size=3, shape=21, color="black", aes(fill=Depth_datum_m)) +
  # geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3"), name=expression("Depth (m below datum)"), limits=c(0,5)) +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.position = "top",
        legend.title = element_text(color="black", size=11),
        legend.text = element_text(color="black", size=10),
        plot.margin = unit(c(0,0,0.5,0), "cm"), # Selecting margins
        panel.border = element_rect(colour="black", fill=NA, linewidth=1)) + # Adding in outer border
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



# List of relevant sites with 2022 data collection 
studysites <- as.vector(c("Between Scotts and Bradys","Bordelais Island","Cable Beach (Blow Hole)","Danvers Danger Rock","Dixon Island Back (Bay)","Dodger Channel 1","Dodger Channel 2","Ed King East Inside","Flemming 112","Flemming 114","Less Dangerous Bay","Nanat Bay","North Helby Rock","Ross Islet 2","Ross Islet Slug Island","Second Beach","Second Beach South","Swiss Boy","Taylor Rock","Turf Island 2","Tzartus 116","Wizard Islet North","Wizard Islet South"))
  
plotdataclust <- plotdata %>%
  filter(SiteName %in% studysites)

colsmap <- c("#225555", "#4477aa", "#997700", "#e4632d")

tiff(file="./MSc_plots/MapClust4.tiff", height = 6, width = 8, units = "in", res=400)

## plotting the kelp cluster groups map (discrete fill)
data_map <- ggplot(land)+
  geom_sf(fill = "grey70", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdataclust, size=4, shape=21, color="black", aes(fill=Cluster)) + 
  scale_fill_manual(values=c(colsmap, "grey45"), name="") +
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




## Vancouver Island map

tiff(file="./MSc_plots/MapVanIsland.tiff", height = 5.5, width = 8.5, units = "in", res=600)

## plotting Barkley Sound area (Vancouver Island)
data_map <- ggplot(land)+
  geom_sf(fill = "grey70", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.position = "top",
        legend.title = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=14, face="bold"),
        panel.border = element_rect(colour="black", fill=NA, linewidth=1), # border
        plot.margin = unit(c(0,0.5,0.5,0), "cm")) + # shrinking margins
  north(land, symbol=12, location="topleft") +
  ggsn::scalebar(land,
                 location = "bottomright",
                 transform = F,
                 st.bottom = F,
                 st.size = 3,
                 height = 0.02,
                 dist = 50,
                 dist_unit = "km",
                 model = 'NAD83')
data_map

dev.off()




#### Plotting the broader map inset ----


# Loading the additional required packages
library(rnaturalearth)
library(rnaturalearthdata)


# Loading and defining the world map as sf of countries
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## setting projections at outset
proj <- st_crs(4269) # ESPG code for NAD1983
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")

## setting map extent for BC -> Cali
ymax <- 51.00
ymin <- 20.00
xmax <- -127.00
xmin <- -116.00

## making corners for the area of interest 
corners <- st_multipoint(matrix(c(xmax,xmin,ymax,ymin),ncol=2)) %>% 
  st_sfc(crs=latlong) %>%
  st_sf() %>%
  st_transform(proj) 
plot(corners)


# projecting the world map
world_proj <- world %>%
  st_as_sfc(crs=latlong) %>%
  st_sf() %>%
  st_transform(proj)

# Rebuilding to a valid geometry
world_proj <- world_proj %>%
  st_make_valid()

# cropping to the corners
world_crop <- world_proj %>% 
  st_crop(st_bbox(corners))



# Plotting the map inset 
tiff(file="./MSc_plots/MapInset.tiff", height = 10, width = 4, units = "in", res=400)

data_map <- ggplot(world_crop)+
  geom_sf(fill = "grey65", color = "black") + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text = element_text(color="black", size=12),
        legend.position = "top",
        legend.title = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=14, face="bold"),
        # panel.border = element_rect(colour="black", fill=NA, linewidth=1), # border
        plot.margin = unit(c(0.2,0,0,0), "cm"),
        axis.text = element_blank()) # shrinking margins
# north(world_crop, symbol=12, location="topright") +
# ggsn::scalebar(world_crop,
#                location = "bottomleft",
#                dist=100,
#                dist_unit="km",
#                transform=TRUE,
#                height=0.01,
#                model = 'WGS84')
data_map

dev.off()

#### Putting together Fig. 1 map figure ----

# Calling cluster plot as magick image from 'comm_analyses.R'
mclust <- image_read("./MSc_plots/Clusters_kmeans4.tiff") 
# mclustscl <- image_scale(mclust, "2000") # Scaling it to smaller size if needed

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/MapClust4.tiff") 
# mplotscl <- image_scale(mplot, "2000") # Scaling it to smaller size if needed

# Adding in stack of kelp forest images from 'comm_analyses.R'
mmap <- image_composite(mplot, imgstk, offset="+2590+480")

# Stacking the map & cluster plots together
mfinal <- image_append(c(mmap, mclust), stack=TRUE)

# Saving the final fig
image_write(mfinal, path = "./MSc_plots/PaperFigs/Fig1b.tiff", format = "tiff", density=600)

