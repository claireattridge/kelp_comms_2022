## Barkley Sound site fetch & wave exposure
## Author: Claire Attridge
## Origin date: Jan 2023

# Load base packages
library(tidyverse)
library(MetBrewer)
library(reader)
library(lubridate)
library(ggpubr)
library(ggspatial)
library(sf)
library(rgeos)
library(windfetch)

#

#### Map prep ----

## setting projections at outset
proj <- st_crs(3005) # ESPG code for BC Albers projection
latlong <- st_crs(4326) # for baseline/existing WGS84 ("+proj=longlat +datum=WGS84 +no_defs")

## setting map extent for Bamfield (setting a wider range than in previous plots - capturing the Broken Group islands for these calculations)
ymax <- 48.98
ymin <- 48.78
xmax <- -125.00
xmin <- -125.45

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


#### Site data prep ----

# Reading in general site data
sites <- read.csv("./MSc_data/Data_new/Final_site_list_2022.csv") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon), 
         Estimated_dens=factor(Estimated_dens, levels=c("High", "Med", "Low", "None")),
         Surveyed = as.factor(Surveyed))


# selecting for only the surveyed study sites
sites <- sites %>%
  filter(Surveyed == "YES")


# nudging the lat & lon values for sites that register as 'on land' in the high res coastal map
sites <- sites %>%
  mutate(Lat = ifelse(sites$SiteName == "Swiss Boy", (sites$Lat - 0.0001), 
                        ifelse(sites$SiteName == "Taylor Rock", (sites$Lat + 0.0001), sites$Lat)))

# shifting the lat & lon values for sites that were inaccurately placed gps points
sites <- sites %>%
  mutate(Lat = ifelse(sites$SiteName == "Between Scotts and Bradys", as.numeric(48.83257), sites$Lat)) %>%
  mutate(Lon = ifelse(sites$SiteName == "Between Scotts and Bradys", as.numeric(-125.1490), sites$Lon))

sites <- sites %>%  
  mutate(Lat = ifelse(sites$SiteName == "Less Dangerous Bay", as.numeric(48.87545), sites$Lat)) %>%
  mutate(Lon = ifelse(sites$SiteName == "Less Dangerous Bay", as.numeric(-125.0908), sites$Lon))



# converting location data to sf object
sitessf <- sites %>%
  drop_na(Lat) %>%
  dplyr::select(c(SiteName, Lat, Lon, Estimated_dens, Composition, Surveyed)) %>%
  st_as_sf(coords = c(3,2))

st_crs(sitessf) <- latlong # setting 
sitessf <- st_transform(sitessf,proj) # make sure you assign this, otherwise it won't use the transformed project crs

# cropping to only points within map limits
sitessf <- sitessf %>%
  st_crop(st_bbox(corners))



### Filtering out the 'no kelp' sites (Less Dangerous Bay, Wizard Islet N)
### Should go from 23 to 20 sites
sitessf <- sitessf %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  filter(SiteName != "Less Dangerous Bay") %>%
  filter(SiteName != "Wizard Islet North")



# locking in factor level order
# renaming 'SiteName' to 'Names' so that the windfetch package registers the names
sitessf <- sitessf %>%
  mutate(Names = as.factor(SiteName), 
         Composition = as.factor(Composition)) 


#### Fetch calculation ----

### Running the fetch function (windfetch)
fetch <- windfetch(polygon_layer=land_bs, site_layer=sitessf, max_dist=60, n_directions=8)
# 60 km max dist (the offshore wind buoy we will use is ~ 60 km offshore -> 'La Perouse')
# 32 directions (8 fetch measures / quadrant of 90 deg. - N, E, S, W)
# *Fetch output is in meters though


## visualizing the fetch
plot(land_bs, col="darkgrey")
plot(fetch, col="black", add=TRUE)


fetch_sf <- as_sf(fetch) # turning into sf object to capture all 36 directional measures

fetchtab <- as.data.frame(fetch_sf) %>% # into a standard dataframe
  dplyr::select(-geometry) %>% # will add the geometry back in later by site name
  mutate(site_name = as.factor(site_name))

# Saving fetch info as a .csv file
write.csv(fetchtab, "./MSc_data/Data_new/fetch_2022.csv", row.names=F)


#### Wave & wind processing ----
 
## uploading wind and wave data from offshore buoy (La Perouse Bank)
weather <- read.csv("./MSc_data/Data_new/Wind_wave/c46206.csv") %>%
                mutate(DATE = as.POSIXct(DATE, format="%m/%d/%Y %H:%M")) %>% # makes date format
  dplyr::select(-c("WSS.","VWH.","VTP.","WSS..1","ATMS.1", "X")) %>% # removing unnecessary columns 
  dplyr::rename(swh = 'VCAR', # changing to column name that's more intuitive 
                Date = DATE) %>%
  filter(Date >= "2012-04-01") # selecting for data from 2012-2022 (records go up to Apr 2022)

# removing NAs from relevant data columns
weather <- weather %>%
  dplyr::filter(!is.na(swh)) %>% 
  dplyr::filter(!is.na(WDIR)) %>% 
  dplyr::filter(!is.na(WSPD)) 


## isolating wave data
wavedata <- weather %>%
  dplyr::select(LATITUDE, LONGITUDE, Date, swh) %>% #selecting for only necessary columns
  dplyr::mutate(id = row.names(.)) %>% # creating id column for later manipulations
  dplyr::arrange(Date) # sorting data by time
wavedata <- as.data.frame(wavedata) # from tibble to dataframe


## isolating wind data
winddata <- weather %>%
  dplyr::select(LATITUDE, LONGITUDE, Date, WDIR, GSPD) %>% # for wind direction & gust speeds
  dplyr::mutate(id = row.names(.)) %>% #creating id column for later manipulations
  dplyr::mutate(winddir_10s = (WDIR/10), windspd = GSPD) # converting to 10s of degs
winddata <- winddata %>%
  dplyr::mutate(winddir_r = round(winddir_10s, digits=0), windspd = GSPD) %>% # to nearest 10th deg
  dplyr::select(-c(WDIR,GSPD, winddir_10s)) %>% # removing unnecessary columns 
  dplyr::arrange(Date) # sorting data by time
winddata <- winddata %>% 
  mutate(id = as.numeric(id)) # make id column numeric


# cleaning up the joined wind & wave data
windwave <- merge(winddata, wavedata) %>%
  mutate(winddir_r = dplyr::case_when(
    winddir_r == 36 ~ 0, # making the 360 deg equal to the 0 deg (are the same direction)
    TRUE ~ winddir_r
  ))
wavey <- windwave %>% 
  dplyr::filter(swh > 2) %>% # significant wave height as waves greater than 2 m
  dplyr::select(-id)
j <- nrow(wavey)


# data manipulations to prepare for REI calcs
wind_36 <- wavey %>%
  dplyr::group_by(winddir_r) %>%
  dplyr::summarize(avg_wind = mean(quantile(windspd, .9)), # top 10% for wind speeds
                   frequency = dplyr::n()) 
wind_36$frequency <- as.numeric(wind_36$frequency)

wind_36 <- wind_36 %>% 
  dplyr::mutate(per_freq = wind_36$frequency/j) %>% # j is the number of rows of 'wavey' (data frame prior to grouping into wind direction bins)
  dplyr::select(-frequency) %>% 
  dplyr::rename("theta" = "winddir_r") %>% 
  mutate(theta = theta*10)

wind_36 <- wind_36 %>%
  mutate(theta = as.factor(theta)) %>%
  as.data.frame()


# data manipulations to prepare for wind rose plots
wavey_deg <- wavey %>% 
  mutate(winddir_r = winddir_r*10) #to full degree units

## making the wind rose
library(openair)

openair::windRose(wavey_deg, 
         ws = "windspd", 
         wd = "winddir_r", 
         angle = 10, # number of degrees for each petal
         paddle = F, 
         cols = "default",
         ws.int = 3,
         breaks = 5,
         grid.line = 5)



#### REI calculations ----

## loading fetch data and wind data for the study sites!

# loading & mutating fetch data frame to have coordinates associated
fetch_df <- read_csv("./MSc_data/Data_new/fetch_2022.csv") %>% # loading the fetch by site
  mutate(SiteName = as.factor(site_name)) # renaming to match for merging
sites_loc <- sites %>% # selecting only relevant location columns for sites
  dplyr::select(SiteName, Lat, Lon)
fetch_loc <- fetch_df %>% 
  left_join(sites_loc, by="SiteName") %>% # joining by the common column
  dplyr::select(-site_name) %>% # removing unnecessary naming column
  st_as_sf(coords = c(6,5)) # specifying columns for geometry (lat & lon)
st_crs(fetch_loc) <- latlong # setting crs of the coordinates as WGS84
fetch_loc <- fetch_loc %>%
  st_transform(proj) %>% # transforming again to NAD83 BC Albers
  mutate(theta = as.factor(directions)) %>% # name to match the wind data frame
  dplyr::select(-directions)


### REI calculations
exp_df <- merge(x=wind_36, y=fetch_loc, by = "theta") %>%
  dplyr::select(avg_wind, per_freq, fetch, geometry, SiteName) %>% 
  dplyr::mutate(exposure_bin = (avg_wind*per_freq*fetch)) # calcuating REI for each wind bin
st_geometry(exp_df) <- exp_df$geometry

## final averaging for REI / site
exp_summary <- exp_df %>% 
  dplyr::group_by(SiteName) %>% 
  st_as_sf() %>%
  dplyr::summarise(exp_36 = sum(exposure_bin)) %>% #sums REI per wind bin for 1 REI per site
  dplyr::rowwise() %>%
  mutate(x = st_coordinates(geometry)[1],
         y = st_coordinates(geometry)[2]) %>%
  dplyr::ungroup()

exp_summary_df <- data.frame(exp_summary) %>% # returning to a dataframe from sf
  dplyr::select(-geometry)

# saving a .csv file of the final REI site values
write.csv(exp_summary_df, "./MSc_data/Data_new/REI_2022.csv", row.names=F)



tiff(file="./MSc_plots/SuppFigs/SiteREI.tiff", height = 6, width = 9, units = "in", res=400)

## plotting out the trend in site specific REI
ggplot() +
  geom_point(data=exp_summary, size=3, alpha=0.9, aes(x=reorder(SiteName,exp_36), y=exp_36)) +
  theme_classic() +
  scale_y_continuous(breaks=c(0, 5000,10000, 15000, 20000, 25000)) +
  theme(
    axis.text.x = element_text(color="black", size="9", angle=45, vjust=0.89, hjust=0.85),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.y = element_text(color="black", size="11", vjust=3)) +
  ylab("REI") + xlab("")

dev.off()
