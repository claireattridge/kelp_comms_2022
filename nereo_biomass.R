## Barkley Sound N. luetkeana biomass data
## Author: Claire Attridge
## Origin date: September 2022

setwd("C:/Users/Claire/Desktop/MSc/Thesis")

library(tidyverse)
library(gam)

nereo <- read.csv("./MSc_data/Data_new/nereo_biomass_2022.csv") 

# converting from circumference to diameter
nereo <- nereo %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159))

# plotting sub-bulb (cm) to biomass (g) relationship
ggplot(data=nereo, aes(x=Sub_circ_cm, y=Biomass_g)) +
  geom_point(size=3, shape=21) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=3), se=F) +
  scale_fill_manual(values=c("green", "orange")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size="9.5"),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("Sub bulb (cm)") + ylab("Biomass (g)")

# generating gam model of sub-bulb (cm) to biomass (g) relationship
model <- gam(Biomass_g~s(Sub_diam_cm, bs = "cs", k=2), family="poisson", data=nereo)
summary(model)
