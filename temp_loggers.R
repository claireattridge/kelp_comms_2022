## Barkley Sound HOBO temperature logger processing
## Author: Claire Attridge
## Origin date: October 2022

library(tidyverse)
library(scales)
library(data.table)

#### Data cleaning ----

# Setting pathway to Temp_loggers folder as a variable
mypath <- "C:/Users/Claire/Desktop/MSc/Thesis/kelp_comms_2022/MSc_data/Data_new/Temp_loggers"

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
temp <- as.data.frame(temp)

# Cleaning up the lead & butt ends based on deployment timelines
tempf <- temp %>%
  filter(fifelse(SiteName == "Ohiat", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-21"),
          fifelse(SiteName == "Tzartus 116", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
            fifelse(SiteName == "Nanat Bay", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-17"),
              fifelse(SiteName == "Between Scotts and Bradys", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-20"),
                fifelse(SiteName == "Bordelais Island", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-20"),
                  fifelse(SiteName == "Ross Islet 2", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-12"),
                    fifelse(SiteName == "Wizard Islet South", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-19"),
                      fifelse(SiteName == "North Helby Rock", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-21"),
                        fifelse(SiteName == "Turf Island 2", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-20"),
                          fifelse(SiteName == "Robbers Passage 2", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                            fifelse(SiteName == "Flemming 114", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-19"),
                              fifelse(SiteName == "Flemming 112", Date >= as.POSIXct("2022-06-30") & Date <= as.POSIXct("2022-09-18"),
                                fifelse(SiteName == "Dodger Channel 1", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-13"),
                                  fifelse(SiteName == "Cable Beach (Blow Hole)", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-14"),
                                    fifelse(SiteName == "Ed King East Inside", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-13"),
                                      fifelse(SiteName == "Ross Islet Slug Island", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-20"),
                                        fifelse(SiteName == "Taylor Rock", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-20"),
                                          fifelse(SiteName == "Dodger Channel 2", Date >= as.POSIXct("2022-06-23") & Date <= as.POSIXct("2022-09-13"),
                                            fifelse(SiteName == "Second Beach South", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-13"),
                                              fifelse(SiteName == "3 Tree Island 2", Date >= as.POSIXct("2022-06-25") & Date <= as.POSIXct("2022-09-21"),
                                                fifelse(SiteName == "Second Beach", Date >= as.POSIXct("2022-06-24") & Date <= as.POSIXct("2022-09-14"),
                                                  fifelse(SiteName == "Wizard Islet North", Date >= as.POSIXct("2022-06-22") & Date <= as.POSIXct("2022-09-19"),
                                                    fifelse(SiteName == "Dixon Island Back (Bay)", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-19"),
                                                      fifelse(SiteName == "Past Roquefoil Bay", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-19"),
                                                        fifelse(SiteName == "Danvers Danger Rock", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-17"),
                                                          fifelse(SiteName == "Swiss Boy", Date >= as.POSIXct("2022-07-22") & Date <= as.POSIXct("2022-09-18"),
                                                            fifelse(SiteName == "Less Dangerous Bay", Date >= as.POSIXct("2022-07-01") & Date <= as.POSIXct("2022-09-17"),
                                                                    Date == Date))))))))))))))))))))))))))))
                         
#### Grouping & adding depths ----

# Finding site average, max, min
tempgrp <- tempf %>%
  group_by(SiteName) %>%
  summarise(Tempave = mean(Temp, na.rm=T), SD_Tempave = sd(Temp, na.rm=T), Tempmax = max(Temp, na.rm=T), Tempmin = min(Temp, na.rm=T)) %>%
  arrange(desc(Tempmax)) # Arranging descending order

tempgrp <- as.data.frame(tempgrp) # Making into a dataframe

## Loading temp logger info sheet
info <- read_csv("./MSc_data/Data_new/temp_logger_info.csv") %>%
  dplyr::select(c(ID, SiteName, Depth_logger_m, Depth_logger_datum_m)) %>%
  mutate(SiteName = as.factor(SiteName), ID = as.factor(ID)) %>%
  filter(ID != "L31" & ID != "L32" & ID != "L33") %>% # Removing the site deployment loggers
  as.data.frame()

tempgrp <- merge(tempgrp, info, by= "SiteName", all = TRUE)


#### Regressions of data ----


## Site depth series
ggplot() +
  geom_point(data=tempgrp, size=3, alpha=0.1, aes(x=reorder(SiteName,Depth_logger_m), y=Depth_logger_m)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size="9.5", angle=45, vjust=0.89, hjust=0.85),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.y = element_text(color="black", size="11", vjust=3))


## All temperature series over time
ggplot(tempf, aes(x=Date, y=Temp, group=SiteName, color=SiteName)) +
  geom_line() +
  scale_x_datetime(labels = date_format("%b %d"), breaks = date_breaks("weeks")) + 
  xlab("Date") +
  theme_classic()

## Temperature against deployment depths
ggplot(tempgrp, aes(x=Depth_logger_datum_m, y=Tempmax)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x, se=T) + # Lm
  scale_x_reverse() +
  theme_classic()





#### Maps of temperature ----

# TBA
