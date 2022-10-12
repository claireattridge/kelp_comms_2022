## Barkley Sound N. luetkeana biomass data
## Author: Claire Attridge
## Origin date: September 2022

library(tidyverse)
library(mgcv) # For gam models

nereo <- read.csv("./MSc_data/Data_new/nereo_biomass_2022.csv") 

# converting from circumference to diameter
nereo <- nereo %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159))

# plotting sub-bulb (cm) to biomass (g) relationship
ggplot(data=nereo, aes(x=Sub_diam_cm, y=Biomass_g, color=SiteName, group=SiteName)) +
  geom_point(size=3, shape=21) +
  # geom_smooth(method="gam", formula = y ~ s(x, bs="cs", k=3), se=F) + # GAM
  geom_smooth(method="lm", formula = y ~ x + I(x^2), se=F) + # Lm quadratic
  # geom_smooth(method="lm", formula = y ~ x, se=F) + # Lm
  scale_fill_manual(values=c("green", "orange")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size="9.5"),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("Sub-bulb diameter (cm)") + ylab("Biomass (g)")

# # generating gam model of sub-bulb (cm) to biomass (g) relationship
# model <- mgcv::gam(Biomass_g~s(Sub_diam_cm, bs="cs", k=3), family="poisson", data=nereo)
# summary(model)

# generating linear model with quadratic term of sub-bulb (cm) to biomass (g) relationship
model <- lm(Biomass_g ~ Sub_diam_cm + I(Sub_diam_cm^2), data=nereo)
summary(model)
coef(model)

## Model equation: Biomass = 150.7966(Sub-bulb)^2 -216.2721(Sub-bulb) + 315.0124 ##

# Testing out model predictions
predict(model, data.frame(Sub_diam_cm=0.509))
