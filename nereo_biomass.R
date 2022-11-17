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


tiff(file="./MSc_plots/Nereo_relation_tog.tiff", height = 5, width = 7, units = "in", res=600)

## Plot of sub-bulb (cm) to biomass (g) relationship by all sites
n0 <- ggplot(data=nereo, aes(x=Sub_diam_cm, y=Biomass_g)) +
  geom_point(size=3, shape=21) +
  geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T) + # Lm quadratic
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
    axis.title.x = element_text(color="black", size="11"), # For plotting purposes
    axis.title.y = element_text(color="black", size="11")) +
  guides(color=guide_legend(override.aes=list(fill=NA))) + # Removing the grey background on legend from the 'se'
  xlab("Sub-bulb diameter (cm)") + ylab("Biomass (g)") +
  geom_text(aes(1.5,5750, label=(paste(expression("y = 150.7966 x"^2*" - 216.2721 x + 315.0124")))),parse = TRUE, size=3.5) +
  geom_text(aes(1.35,5200, label=(paste(expression("R"^2*" = 0.81")))),parse = TRUE, size=3.5) # Adding in the R squared value
n0

dev.off()




tiff(file="./MSc_plots/Nereo_relation_p1.tiff", height = 5, width = 7, units = "in", res=600)

## Plot of sub-bulb (cm) to biomass (g) relationship split by sites
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
### We may have to split! 

# Making dataframe w only Second Beach South
nereo2 <- nereo %>%
  filter(SiteName == "Second Beach South")

#Making dataframe w not Second Beach South
nereo3 <- nereo %>%
  filter(SiteName != "Second Beach South")




## Plot of sub-bulb (cm) to biomass (g) relationship for Second Beach South
n2 <- ggplot(data=nereo2, aes(x=Sub_diam_cm, y=Biomass_g)) +
  geom_point(size=3, shape=21) +
  geom_smooth(method="lm", formula = y ~ x + I(x^2), se=T, color="black") + # Lm quadratic
  theme_classic() +
  scale_y_continuous(breaks=c(0,2000,4000,6000), limits=c(-500,6300)) +
  scale_x_continuous(breaks=c(0,2,4,6), limits=c(0,6.5)) +
  theme(
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10)
    ) +
  xlab("") + ylab("") +
  geom_text(aes(2.8,5750, label=(paste(expression("y = 99.24119 x"^2*" + 238.81054 x - 26.33073")))),parse = TRUE, size=3.5) + # Adding in model equation
  geom_text(aes(0.55,5200, label=(paste(expression("R"^2*" = 0.93")))),parse = TRUE, size=3.5) # Adding in the R squared value
n2


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
    ) +
  xlab("") + ylab("") +
  geom_text(aes(2.85,5750, label=(paste(expression("y = 43.21482 x"^2*" + 169.49876 x + 60.37447")))),parse = TRUE, size=3.5) + # Adding in model equation
  geom_text(aes(0.55,5200, label=(paste(expression("R"^2*" = 0.81")))),parse = TRUE, size=3.5) # Adding in the R squared value
n3



### Going to arrange all of these plots into one now
### Label them uniquely too

# Creating the bottom row arrangement
bottom <- plot_grid(n2, n3, axis="b", align="hv", ncol=2, labels=c("B", "C"))
# Putting together with top row plot
nall <- plot_grid(n1, bottom, labels=c("A", ""), label_size=12, ncol=1, rel_widths=c(2,2))
nall

# Creating the axis labels
y.grob <- textGrob("Biomass (g)", 
                   gp=gpar(col="black", fontsize=13), rot=90)

x.grob <- textGrob("Sub-bulb diameter (cm)", 
                   gp=gpar(col="black", fontsize=13))


tiff(file="./MSc_plots/Nereo_relation_pALL.tiff", height = 7, width = 8, units = "in", res=600)

# Adding the axis titles onto the plot arrangement
grid.arrange(arrangeGrob(nall, left = y.grob, bottom = x.grob))

dev.off()


#### Models of relationships ----

###
# Quadratic model of sub-bulb (cm) to biomass (g) relationship for all sites 
model <- lm(Biomass_g ~ Sub_diam_cm + I(Sub_diam_cm^2), data=nereo)
summary(model)
coef(model)

# Model equation: Biomass = 150.7966(Sub-bulb)^2 -216.2721(Sub-bulb) + 315.0124 
# R2 of 0.82

# Testing out model predictions
predict(model, data.frame(Sub_diam_cm=0.509))
###


###
# Quadratic model of sub-bulb (cm) to biomass (g) relationship for Second Beach South (Panel B)
model2 <- lm(Biomass_g ~ Sub_diam_cm + I(Sub_diam_cm^2), data=nereo2)
summary(model2)
coef(model2)

# Model equation: Biomass = 99.24119(Sub-bulb)^2 238.81054(Sub-bulb) - 26.33073
# R2 of 0.93

# Testing out model predictions
predict(model2, data.frame(Sub_diam_cm=0.509))
###


###
# Quadratic model of sub-bulb (cm) to biomass (g) relationship for Ed King E Inside & Swiss Boy (Panel C)
model3 <- lm(Biomass_g ~ Sub_diam_cm + I(Sub_diam_cm^2), data=nereo3)
summary(model3)
coef(model3)
# R2 of 0.81

# Model equation: Biomass = 43.21482(Sub-bulb)^2 + 169.49876(Sub-bulb) + 60.37447

# Testing out model predictions
predict(model3, data.frame(Sub_diam_cm=0.509))
###









