## Barkley Sound N. luetkeana biomass data
## Author: Claire Attridge
## Origin date: September 2022

library(tidyverse)
library(MetBrewer)
library(mgcv) # For gam models
library(cowplot)
library(gridExtra)

#### Data cleaning ----

nereo <- read.csv("./MSc_data/Data_new/nereo_biomass_2022.csv") 

# converting from circumference to diameter
nereo <- nereo %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159))

#### Plots of relationships ----

tiff(file="./MSc_plots/Nereo_relation.tiff", height = 5, width = 7, units = "in", res=600)

## Plot of sub-bulb (cm) to biomass (g) relationship by all sites
n1 <- ggplot(data=nereo, aes(x=Sub_diam_cm, y=Biomass_g, color=SiteName, group=SiteName, linetype=SiteName)) +
  geom_point(size=3, shape=21) +
  # geom_smooth(method="gam", formula = y ~ s(x, bs="cs", k=3), se=F) + # GAM
  geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T) + # Lm quadratic
  # geom_smooth(method="lm", formula = y ~ x, se=F) + # Lm
  scale_color_manual(values=met.brewer("Lakota", 3)) +
  theme_classic() +
  scale_x_continuous(breaks=c(0,2,4,6), limits=c(0,6.5)) +
  scale_y_continuous(breaks=c(0,2000,4000,6000), limits=c(-500,6300)) +
  theme(
    legend.position = c(0.2,0.8),
    legend.title = element_blank(),
    legend.text = element_text(color="black", size=11),
    axis.text.x = element_text(color="black", size="9.5"),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.x = element_blank(), # For plotting purposes
    axis.title.y = element_text(color="black", size="11", vjust=3)) +
    guides(color=guide_legend(override.aes=list(fill=NA))) + # Removing the grey background on legend from the 'se'
  xlab("") + ylab("")
n1

dev.off()



### So it looks like Second Beach South relationship is unique from the other two sites
### This suggests that the sub-bulb relationship is likely depth or temp rather than spatially dependent
### We will have to split! 

# Making dataframe w only Second Beach South
nereo2 <- nereo %>%
  filter(SiteName == "Second Beach South")

#Making dataframe w not Second Beach South
nereo3 <- nereo %>%
  filter(SiteName != "Second Beach South")



tiff(file="./MSc_plots/Nereo_relation.tiff", height = 5, width = 7, units = "in", res=600)

## Plot of sub-bulb (cm) to biomass (g) relationship for Second Beach South
n2 <- ggplot(data=nereo2, aes(x=Sub_diam_cm, y=Biomass_g)) +
  geom_point(size=3, shape=21) +
  geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T, color="black") + # Lm quadratic
  theme_classic() +
  scale_y_continuous(breaks=c(0,2000,4000,6000), limits=c(-500,6300)) +
  scale_x_continuous(breaks=c(0,2,4,6), limits=c(0,6.5)) +
  theme(
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=11, vjust=3)) +
  xlab("") + ylab("")
n2

dev.off()


tiff(file="./MSc_plots/Nereo_relation.tiff", height = 5, width = 7, units = "in", res=600)

## Plot of sub-bulb (cm) to biomass (g) relationship for Ed King E Inside & Swiss Boy
n3 <- ggplot(data=nereo3, aes(x=Sub_diam_cm, y=Biomass_g)) +
  geom_point(size=3, shape=21) +
  geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T, color="black") + # Lm quadratic
  theme_classic() +
  scale_y_continuous(breaks=c(0,2000,4000,6000), limits=c(-500,6300)) +
  scale_x_continuous(breaks=c(0,2,4,6), limits=c(0,6.5)) +
  theme(
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_blank(), # For plotting purposes
    # axis.text.y = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=11, vjust=3)) +
  xlab("") + ylab("")
n3

dev.off()


### Going to arrange all of these plots into one now
### Label them uniquely too

# Creating the bottom row arrangement
bottom <- plot_grid(n2, n3, axis="b", align="hv", ncol=2, labels=c("B", "C"))
# Putting together with top row plot
nall <- plot_grid(n1, bottom, labels=c("A", ""), label_size=12, ncol=1, rel_widths=c(2,2))
nall

# Creating the axis labels
y.grob <- textGrob("Biomass (g)", 
                   gp=gpar(col="black", fontsize=12), rot=90)

x.grob <- textGrob("Sub-bulb diameter (cm)", 
                   gp=gpar(col="black", fontsize=12))

# Adding the axis titles onto the plot arrangement
grid.arrange(arrangeGrob(nall, left = y.grob, bottom = x.grob))


#### Models of relationships ----

# generating linear model with quadratic term of sub-bulb (cm) to biomass (g) relationship
model <- lm(Biomass_g ~ Sub_diam_cm + I(Sub_diam_cm^2), data=nereo)
summary(model)
coef(model)

## Model equation: Biomass = 150.7966(Sub-bulb)^2 -216.2721(Sub-bulb) + 315.0124 ##

# Testing out model predictions
predict(model, data.frame(Sub_diam_cm=0.509))






