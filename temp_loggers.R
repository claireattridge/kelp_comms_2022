## Barkley Sound Electric Blue temperature logger processing
## Author: Claire Attridge
## Origin date: October 2022

# Loading base packages
library(tidyverse)
library(scales)
library(data.table)
library(MetBrewer)
library(wesanderson)

#

#### 2022 temperature dataset ----

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
  theme_classic() +
  theme(legend.position="none",
        axis.line=element_line(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text.x=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=12),
        axis.text.y=element_text(color="black", size=9),
        axis.text.x=element_text(color="black", size=8, angle=55, vjust=0.89, hjust=0.85),
        panel.grid.major.y=element_line(color="grey66"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
        ) +
  facet_wrap(~SiteName, ncol=3)


dev.off()



## Temperature against deployment depths
ggplot(tempgrp, aes(x=Depth_logger_datum_m, y=Tempave)) +
  geom_point() +
  geom_smooth(method="lm", formula = y ~ x, se=T) + # Lm
  # scale_x_reverse() +
  theme_classic()


#

