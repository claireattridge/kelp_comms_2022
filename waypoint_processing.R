## 
## Author: Claire Attridge
## Origin date: May 2022 


library()
library(plotKML)
library(tidyr)

r <- readGPX("./MSc_data/Data_new/BMSCLOGGERSV2.gpx", waypoints=T, tracks=F, routes=F)

df <- data.frame(matrix(unlist(r), nrow=length(r), byrow=TRUE)) %>% # Turning list into df
  add_column(ID = c("Lon", "Lat", "Date", "Waypoint", "Shape")) %>% # Setting names for my rows
  mutate(ID = as.factor(ID))

wp <- df %>%
  rownames_to_column() %>% # Necessary when the original df has meaningful row names
  gather(variable, value, -ID) %>% # Gathering all values by all variables 
  spread(ID, value) # Spreading out based on the meaningful row names


