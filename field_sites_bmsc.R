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
library(elementalist)

library(imager)
library(magick)


#

#### Spatial package prep work *Only need to do this the 1st time ----

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


#### Map prep work *Only need to do this the 1st time ----


## setting projections at outset
proj <- st_crs(3005) # ESPG code for BC Albers projection
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")



### Barkley Sound processing

# seting map extent for Bamfield / Barkley Sound
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
sea_bs <- sea_bs %>%
  st_sf() %>% # making into sf object 
  st_transform(proj) # projecting back into ESPG 3005


## saving the cropped size map of land 
## saving the cropped size map of sea
write_sf(land_bs, "./Maps_BC/Hakai_BC/land_Hakai_BarkleyS.shp", overwrite=T)
write_sf(sea_bs, "./Maps_BC/Hakai_BC/sea_Hakai_BarkleyS.shp", overwrite=T)




### Vancouver Island processing

# seting map extent for Vancouver Island
ymax <- 52.0
ymin <- 48.5
xmax <- -120.00
xmin <- -129.00

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
land_vi <- land %>% # cropping to the Van Island corners
  st_crop(st_bbox(corners))

## converting from land area to sea area (because this map is the land not the water)
mask <- st_bbox(land_vi) %>% # generating a box that covers the cropped multipolygon area
  st_as_sfc() # making a simple geometry feature

land_vi <- st_union(land_vi) # unionizing the land layer

sea_vi <- st_difference(mask, land_vi) # removing the intersecting area (i.e., removing sea and leaving land)
sea_vi <- sea_vi %>%
  st_sf() %>% # making into sf object 
  st_transform(proj) # projecting back into ESPG 3005

## making the land back into an sf object
land_vi <- land_vi %>%
  st_sf() %>% # making into sf object 
  st_transform(proj) # projecting back into ESPG 3005


## saving the cropped size map of land 
## saving the cropped size map of sea
write_sf(land_vi, "./Maps_BC/Hakai_BC/land_Hakai_VanIsland.shp", overwrite=T)
write_sf(sea_vi, "./Maps_BC/Hakai_BC/sea_Hakai_VanIsland.shp", overwrite=T)





### BC (Van Island) + USA coast processing

# seting map extent for Vancouver Island
ymax <- 52.0
ymin <- 47.0
xmax <- -120.00
xmin <- -129.00

## making corners for the area of interest 
corners <- st_multipoint(matrix(c(xmax,xmin,ymax,ymin),ncol=2)) %>% 
  st_sfc(crs=latlong) %>%
  st_sf() %>%
  st_transform(proj) 
plot(corners)

## reading in high res map of NE coast 
land <- sf::st_read("./Maps_BC/ne_10m_land/ne_10m_land.shp") %>%
  st_sf() %>%
  st_transform(proj)
land_usa <- land %>% # cropping to the Van Island corners
  st_crop(st_bbox(corners))

# Removing extra features
land_usa <- land_usa %>%
  dplyr::select(-c("scalerank", "min_zoom"))

## Saving the cropped map of land
write_sf(land_usa, "./Maps_BC/ne_10m_land/ne_10m_land_VanIsland.shp", overwrite=T)


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

# Removing the no kelp sites (Less Dangerous Bay & Wizard I North)
sites <- sites %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North")
# Should leave 21 study sites remaining


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


#### Plotting the Barkley Sound & Van Island inset maps ----

# Calling the new cropped land & sea maps for Barkley Sound
land_bs <- sf::st_read("./Maps_BC/Hakai_BC/land_Hakai_BarkleyS.shp") 
sea_bs <- sf::st_read("./Maps_BC/Hakai_BC/sea_Hakai_BarkleyS.shp")


# Adding new column for the legend levels (common names of kelps)
# Adding new columns for just Longitude and for numbering sites
plotdata <- plotdata %>%
  mutate(Kelpdom_comm = as.factor(ifelse(Kelpdom == "Macro", "Giant", "Bull"))) %>%
  merge(sites[,c("SiteName","Lon")], by="SiteName") %>%
  mutate(SiteNum = row_number(Lon))
  

# Second Beach South is plotting on top of Second Beach and blocking it
# Need to swap order for the plot
topframe <- plotdata[15,]


## Barkley Sound map

tiff(file="./MSc_plots/MapForests.tiff", height = 5.5, width = 8.5, units = "in", res=400)

## plotting the kelp species map (discrete fill)
BSdata_map <- ggplot(land_bs)+
  geom_sf(fill = "grey65", color = NA) + # 'color = NA' removes country borders
  theme_minimal(base_size = 16) +
  geom_sf(data = plotdata, shape=21, color="black", aes(size=Area_m2, fill=Kelpdom_comm)) + 
  geom_sf(data = topframe, shape=21, color="black", aes(size=Area_m2, fill=Kelpdom_comm)) + # Making Second Beach point visible
  geom_sf_text(mapping=aes(), data = plotdata, label = plotdata$SiteNum, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=50, nudge_x=-700) + # provides site nums at points
  guides(fill = "none", size = guide_legend(title = expression("Forest area (m"^2*")"))) +
  scale_size_continuous(name="area", breaks = c(250,500,1000,2000,3000)) + # Adding more levels to size legend
  scale_fill_manual(values=met.brewer("Kandinsky", 2), name="") + # Designating colors by kelp spp
  theme(panel.grid.major = element_line(colour = "transparent"), # Hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1.215,0.18),
        axis.text = element_text(color="black", size=10),
        legend.title = element_text(color="black", size=11),
        legend.text = element_text(color="black", size=10),
        # plot.margin = unit(c(0.5,0.5,0.2,0), "cm"), # margins for stand alone
        plot.margin = unit(c(0.2, 5, 0.2, 0.2), "cm"), # margins for adding side images
        panel.border = element_rect_round(colour="black", fill=NA)) + # Adding in outer black border
  north(land_bs, 
        location="bottomright", 
        anchor = c(x = 1069400, y = 421508),
        symbol=12, 
        scale=0.08) +
  ggsn::scalebar(land_bs,
               location = "bottomright",
               anchor = c(x = 1069300, y = 420908),
               transform = F,
               st.bottom = F,
               st.size = 3,
               height = 0.01,
               dist = 1.5,
               dist_unit = "km",
               model = 'NAD83') +
  annotate("text", x = 1069400, y = 433958, label = "(b)", size=3, fontface=2)
  
BSdata_map

dev.off()



# tiff(file="./MSc_plots/MapDepth.tiff", height = 7, width = 10, units = "in", res=600)
# 
# ## plotting the kelp density map (continuous fill)
# data_map <- ggplot(land)+
#   geom_sf(fill = "grey53", color = NA) + # color = NA will remove country border
#   theme_minimal(base_size = 16) +
#   geom_sf(data = plotdata, size=3, shape=21, color="black", aes(fill=Depth_datum_m)) +
#   # geom_sf_text(mapping=aes(), data=sitessf, label=sitessf$SiteName, stat="sf_coordinates", inherit.aes=T, size=2.5, nudge_y=100, nudge_x=-1200) +
#   scale_fill_gradientn(colors=met.brewer("VanGogh3"), name=expression("Depth (m below datum)"), limits=c(0,5)) +
#   theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text = element_text(color="black", size=12),
#         legend.position = "top",
#         legend.title = element_text(color="black", size=11),
#         legend.text = element_text(color="black", size=10),
#         plot.margin = unit(c(0,0,0.5,0), "cm"), # Selecting margins
#         panel.border = element_rect(colour="black", fill=NA, linewidth=1)) + # Adding in outer border
#   north(land, symbol=12, location="topleft") +
#   ggsn::scalebar(land,
#                  location = "bottomright",
#                  transform = F,
#                  st.bottom = F,
#                  st.size = 3.5,
#                  height = 0.01,
#                  dist = 1,
#                  dist_unit = "km",
#                  model = 'NAD83')
# data_map
# 
# dev.off()
# 
# 
# 
# # List of relevant sites with 2022 data collection 
# studysites <- as.vector(c("Between Scotts and Bradys","Bordelais Island","Cable Beach (Blow Hole)","Danvers Danger Rock","Dixon Island Back (Bay)","Dodger Channel 1","Dodger Channel 2","Ed King East Inside","Flemming 112","Flemming 114","Less Dangerous Bay","Nanat Bay","North Helby Rock","Ross Islet 2","Ross Islet Slug Island","Second Beach","Second Beach South","Swiss Boy","Taylor Rock","Turf Island 2","Tzartus 116","Wizard Islet North","Wizard Islet South"))
#   
# plotdataclust <- plotdata %>%
#   filter(SiteName %in% studysites)
# 
# colsmap <- c("#225555", "#4477aa", "#997700", "#e4632d")
# 
# tiff(file="./MSc_plots/MapClust4.tiff", height = 6, width = 8, units = "in", res=400)
# 
# ## plotting the kelp cluster groups map (discrete fill)
# data_map <- ggplot(land)+
#   geom_sf(fill = "grey70", color = NA) + # color = NA will remove country border
#   theme_minimal(base_size = 16) +
#   geom_sf(data = plotdataclust, size=4, shape=21, color="black", aes(fill=Cluster)) + 
#   scale_fill_manual(values=c(colsmap, "grey45"), name="") +
#   theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text = element_text(color="black", size=12),
#         legend.position = "top",
#         legend.title = element_text(color="black", size=12),
#         legend.text = element_text(color="black", size=14, face="bold"),
#         panel.border = element_rect(colour="black", fill=NA, linewidth=1), # border
#         plot.margin = unit(c(0,0.5,0.5,-4), "cm")) + # shrinking margins
#   north(land, symbol=12, location="topleft") +
#   ggsn::scalebar(land,
#                  location = "bottomright",
#                  transform = F,
#                  st.bottom = F,
#                  st.size = 3.5,
#                  height = 0.01,
#                  dist = 1,
#                  dist_unit = "km",
#                  model = 'NAD83')
# data_map
# 
# dev.off()





# Calling the new cropped land & sea maps for Vancouver Island
land_vi <- sf::st_read("./Maps_BC/Hakai_BC/land_Hakai_VanIsland.shp") 
sea_vi <- sf::st_read("./Maps_BC/Hakai_BC/sea_Hakai_VanIsland.shp")

# Calling the new cropped land map for Vancouver Island + USA coast
land_vi_usa <- sf::st_read("./Maps_BC/ne_10m_land/ne_10m_land_VanIsland.shp") 


## Vancouver Island map

## plotting Barkley Sound area on Vancouver Island
VIdata_map <- ggplot(land_vi_usa)+
  geom_sf(fill = "grey70", color = NA) + # color = NA will remove country border
  theme_minimal(base_size = 16) +
  theme(panel.grid.major = element_line(colour = "transparent"), # hiding graticules
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect_round(fill = 'white', colour = 'white'), # forcing white background
        legend.position = "top",
        legend.title = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=14, face="bold"),
        panel.border = element_rect_round(colour="black", fill=NA), # border
        plot.margin = unit(c(0,0.5,0.5,0), "cm"),) +  # shrinking margins
  # north(land, symbol=12, location="topleft") +
  # ggsn::scalebar(land,
  #                location = "bottomright",
  #                transform = F,
  #                st.bottom = F,
  #                st.size = 3,
  #                height = 0.02,
  #                dist = 50,
  #                dist_unit = "km",
  #                model = 'NAD83') +
  annotate("text", x = 1372000 , y = 760000, label = "(a)", size=3, fontface=2)
  
VIdata_map



tiff(file="./MSc_plots/MapVanIsland.tiff", height = 5.5, width = 8.5, units = "in", res=600)

## Adding bounding box to highlight Barkley Sound region
VIdata_map <- 
  VIdata_map +
  geom_rect(
    xmin = 1054306.2934839332,
    ymin = 420893.0657236207,
    xmax = 1069898.5088127996,
    ymax = 434255.0779045987,
    fill = NA, 
    colour = "black",
    linewidth = 0.3
  )
VIdata_map

dev.off()





## Adding the Vancouver Island inset to the Barkley Sound kelp spp map

tiff(file="./MSc_plots/MapBarkleyS_Inset.tiff", height = 5, width = 8, units = "in", res=300)

ggdraw(BSdata_map) +
  draw_plot(
    {
      VIdata_map + 
        coord_sf(
          xlim = c(808809.6, 1410805),
          ylim = c(225734.9, 793926.2),
          expand = TRUE) +
        theme(legend.position = "none")
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    # x = 0.064, # For stand alone
    x = 0.0112, # For adding side images
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    # y = 0.55, # For stand alone
    y = 0.55, # For adding side images
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.45, 
    height = 0.4)

dev.off()

#

#### Putting together the Fig. 1 map (+ kelp forest photos and kelp density plot) ----

# https://cran.r-project.org/web/packages/magick/vignettes/intro.html


### MAKING THE MAPS + PHOTOS UPPER PANEL

## Calling map + inset plot as magick image
mmap <- image_read("./MSc_plots/MapBarkleyS_Inset.tiff") 
# mclustscl <- image_scale(mclust, "2000") # Scaling it to smaller size if needed


## Forming image stack to go beside map

# Giant
img1 <- image_read('./MSc_plots/Giantforest_crop.JPG')
img1bd <- img1 %>%
  magick::image_border(color="#ce9642", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)
img1bd <- img1bd %>%
  magick::image_border(color="black", "8x8") %>% # Adding thinner black border
  magick::image_border(color="white", "40x40") # Adding white border to create space

# Bull
img2 <- image_read('./MSc_plots/Bullforest_crop.JPG')
img2bd <- img2 %>%
  magick::image_border(color="#3b7c70", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)
img2bd <- img2bd %>%
  magick::image_border(color="black", "8x8") %>% # Adding thinner black border
  magick::image_border(color="white", "40x40") # Adding white border to create space
  

## Combining the images into a single vector
imgs <- c(img1bd,img2bd)


## Stacking the images vertically
imgstk <- magick::image_append(image_scale(imgs, "500"), stack=TRUE)


## Adding the stack of kelp forest images to map margin
mplot <- image_composite(mmap, imgstk, offset="+1850+72")


## Saving the final fig
image_write(mplot, path = "./MSc_plots/PaperFigs/uppermap.tiff", format = "tiff", density=600)




### MAKING THE KELP DENSITY LOWER PANEL


## Calling the kelp data sheet
kelpdat <- read_csv("./MSc_data/Data_new/kelp_metrics_2022.csv")

# Adding the longitude column to the kelpdat sheet
# Adding the site number column
kelpdat <- kelpdat %>%
  merge(sites[,c("SiteName","Lat", "Lon", "Kelpdom")], by="SiteName") %>%
  mutate(SiteNum = as.factor(row_number(Lon))) %>%
  mutate(SiteName = as.factor(SiteName), Composition = as.factor(Composition), Lon = as.factor(Lon), Kelpdom=factor(Kelpdom, levels=c("Nereo", "Macro")))



# Pulling out & saving the column of sites #s by longitude
sitenums <- kelpdat %>%
  dplyr::select(SiteName,SiteNum,Lat,Lon)

write.csv(sitenums, "./MSc_data/Data_new/Site_list_numbers.csv", row.names=F)



## Calling the kelp data sheet
kelptog_clean <- read_csv("./MSc_data/Data_new/kelp_rawdata_clean_2022.csv")

# Adding the site number & kelpdom columns from the kelpdat sheet
kelptog_clean <- kelptog_clean %>%
  merge(kelpdat[,c("SiteName","SiteNum", "Kelpdom")], by="SiteName") %>%
  mutate(SiteName = as.factor(SiteName), SiteNum = as.factor(SiteNum), Kelpdom=factor(Kelpdom, levels=c("Nereo", "Macro"))) %>%
  dplyr::select(-c(Surveyor, Date, Time_start))



tiff(file="./MSc_plots/Kelpdens_bylon.tiff", height = 2, width = 6.1, units = "in", res=400)

# The density plot
sitedens <- ggplot() +
  geom_pointrange(data=kelpdat, size=0.6, shape=21, color="grey20", stroke=0.3, linewidth=0.3, aes(x=SiteNum, y=DensityM, ymin=DensityM-DensitySD, ymax=DensityM+DensitySD, fill=Kelpdom)) +
  geom_point(data=kelptog_clean, size=2.6, shape=21, color="grey40", stroke=0.3, alpha=0.2, aes(x=SiteNum, y=Kelp, fill=Kelpdom)) +
  scale_y_continuous(limits=c(-0.8,25), breaks=c(0,6,12,18,24)) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c("Bull", "Giant")) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"), # margins for matching map alignment
    panel.grid.major = element_line(colour = "transparent"), # Hiding the grid
    panel.grid.minor = element_line(colour = "transparent"), # Hiding the grid
    legend.position = "none",
    axis.text.x = element_text(color="black", size="9"),
    axis.text.y = element_text(color="black", size="8"),
    axis.title.y = element_text(color="black", size="9", vjust=3),
    panel.border = element_rect_round(colour="black", fill=NA)) + # Adding in the outer border
  xlab("") + ylab(expression("Kelp density (stipes /m"^2*")"))
sitedens

dev.off()

#



# Calling dens plot as magick image
mdens <- image_read("./MSc_plots/Kelpdens_bylon.tiff") 

# Stacking the map & cluster plots together
mfinal <- image_append(c(mplot, mdens), stack=TRUE)

# Saving the final fig
image_write(mfinal, path = "./MSc_plots/PaperFigs/Fullmap.tiff", format = "tiff", density=300)


#


#### Plotting with the broader map inset (Pacific Northeast) ----


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

#### Putting together Fig. 1 map figure V1 (clusters) ----

# Calling cluster plot as magick image from 'comm_analyses.R'
mclust <- image_read("./MSc_plots/Clusters_kmeans.tiff") 
# mclustscl <- image_scale(mclust, "2000") # Scaling it to smaller size if needed

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/MapClust.tiff") 
# mplotscl <- image_scale(mplot, "2000") # Scaling it to smaller size if needed

# Adding in stack of kelp forest images from 'comm_analyses.R'
mmap <- image_composite(mplot, imgstk, offset="+2590+480")

# Stacking the map & cluster plots together
mfinal <- image_append(c(mmap, mclust), stack=TRUE)

# Saving the final fig
image_write(mfinal, path = "./MSc_plots/PaperFigs/Fig1b.tiff", format = "tiff", density=600)



#
