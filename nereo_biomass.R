## Barkley Sound N. luetkeana biomass data
## Author: Claire Attridge
## Origin date: September 2022

nereo <- read.csv("./MSc_data/Data_new/nereo_biomass_2022.csv") 

nereo <- nereo %>%
  mutate(Sub_diam_cm = (Sub_circ_cm/3.14159))

ggplot(data=nereo, aes(x=Sub_diam_cm, y=Biomass_g)) +
  geom_point(size=3, shape=21) +
  geom_smooth(method="lm", formula = y~x, se=F) +
  scale_fill_manual(values=c("green", "orange")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size="9.5"),
    axis.text.y = element_text(color="black", size="10"),
    axis.title.y = element_text(color="black", size="11", vjust=3)) +
  xlab("Sub bulb (cm)") + ylab("Biomass (g)")

model <- lm(data=nereo, Biomass_g~Sub_diam_cm)
summary(model)
