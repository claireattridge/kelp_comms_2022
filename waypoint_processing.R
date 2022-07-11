## SIMRAD waypoint processing
## Author: Claire Attridge
## Origin date: May 2022 

setwd("C:/Users/Claire/Desktop/MSc/Thesis")

library(tidyverse)
library(plotKML)
library(tidyr)

r <- readGPX("./MSc_data/Data_new/BMSCSITESV3.gpx", waypoints=T, tracks=F, routes=F)

df <- r %>%
  map_df(as_tibble) # Turning list into tibble






