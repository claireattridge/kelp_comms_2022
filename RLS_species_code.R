# Code to play with RLS data
# Author: Em Lim
# Editor: Claire Attridge
# Updated CMA - Feb 2023

# Load packages
library(tidyverse)
library(visreg)
library(ggplot2)
library(RColorBrewer)
library(ggeffects)
theme_set(theme_bw())

# Load data -------------

rls <- read_csv("./MSc_data/Data_new/RLS_2022_KDC_CMA.csv") %>%
  dplyr::select(-1) %>% # cuts the first column which is blank
  dplyr::select(-Inverts) %>% # cuts the inverts column which is just NAs
  filter(Method != 0) %>% # get rid of all method 0's
  slice(2:n()) %>% # cuts the first blank row
  rename(
    site_ID = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`
  )  %>% # Rename columns with spaces
  filter(Species != "Debris - Metal") %>%
  filter(Species != "Debris - Other") %>%
  filter(Species != "Debris - Wood") %>%
  filter(Species != "Debris - Glass") %>%
  filter(Species != "Debris - Fishing Gear") %>%
  filter(Species != "Debris - Zero") %>%
  filter(Species != "Survey not completed") %>%
  filter(Species != "No species found") %>%
  filter(Species != "No species found") %>%
  filter(Species != "NOT PRESENT") # Cut the non-animal species

# Neat data manipulations ------

# Cleaning 
rls <- rls %>%
  mutate(Species = as.factor(Species)) 

# What's the most abundant species?
rls_abundant_species <- rls %>% group_by(Species) %>%
  summarise(sum = sum(Total))  %>% arrange(desc(sum))

# Which site has the most species?
rls_richness <- rls %>%
  group_by(site_name) %>%
  summarize(species_richness = n_distinct(Species))%>%
  arrange(desc(species_richness))
rls_richness <- rls_richness %>%
  rename(SiteName = "site_name")

length(unique(rls$Species))

# where are the most urchins???
# Can swap for any species
rls_urchins <- rls %>%
  filter(Species == "Mesocentrotus franciscanus") %>%
  group_by(site_name) %>%
  summarize(urchins = sum(Total)) %>%
  arrange(desc(urchins))


## Total counts for all given species at all sites
rls_all <- rls %>%
  as.data.frame() %>%
  ungroup() %>%
  group_by(site_name, Species, common_name) %>%
  summarise(Total = sum(Total)) %>%
  rename(SiteName = "site_name")


# saving a .csv file of the total species abundances by site
write.csv(rls_all, "./MSc_data/Data_new/RLS_2022.csv", row.names=F)


## Total counts for any given species at all sites
rls_sp <- rls %>%
  as.data.frame() %>%
  filter(Species == "Haliotis kamtschatkana") %>%
  group_by(site_name) %>%
  summarize(sp = sum(Total)) %>%
  arrange(desc(sp))
rls_sp <- rls_sp %>%
  rename(SiteName = "site_name")
  


# So basically you can manipulate the data to summarize whatever thing you're interested in, and then you can join that summarized data with another df of interest by the common site code/site name :)

# Claire you probably want to stop here what follows is a crime against humanity
# Basically I tried to guess how much the more abundant species would weigh and see where biomass was highest

# Very heinous mass estimates ------
rls_hell <- rls %>%
  mutate(
    mass = ifelse(Species == "Mesocentrotus franciscanus", 300,
                  ifelse(Species == "Pomaulax gibberosus", 60,
                         ifelse(Species == "Patiria miniata", 60,
                                ifelse(Species == "Rhinogobiops nicholsii", 6,
                                       ifelse(Species == "Apostichopus californicus", 695,
                                              ifelse(Species == "Haliotis kamtschatkana", 200,
                                                     ifelse(Species == "Acmaea mitra", 8,
                                                            ifelse(Species == "Ceratostoma foliatum", 14,
                                                                   ifelse(Species == "Orthasterias koehleri", 25,
                                                                          ifelse(Species == "Dermasterias imbricata", 60,
                                                                                 ifelse(Species == "Paguroidea spp.", 6,
                                                                                        ifelse(Species == "Strongylocentrotus purpuratus", 70,
                                                                                               ifelse(Species == "Pisaster ochraceus", 100, # top 16 inverts
                                                                                                      ifelse(Species == "Hexagrammos decagrammus", 600,
                                                                                                             ifelse(Species == "Sebastes flavidus", 1500,
                                                                                                                    ifelse(Species == "Sebastes melanops", 1500,
                                                                                                                           ifelse(Species == "Sebastes caurinus", 1500,
                                                                                                                                  ifelse(Species == "Sebastes maliger", 1500,
                                                                                                                                         ifelse(Species == "Sebastes nebulosus", 1500,
                                                                                                                                                ifelse(Species == "Sebastes spp. juv", 100,
                                                                                                                                                       ifelse(Species == "Jordania zonope", 10,
                                                                                                                                                              ifelse(Species == "Artedius harringtoni", 10,
                                                                                                                                                                     ifelse(Species == "Embiotoca lateralis", 750,
                                                                                                                                                                            ifelse(Species == "Oxylebius pictus", 500,
                                                                                                                                                                                   ifelse(Species == "Rhacochilus vacca", 500,
                                                                                                                                                                                          ifelse(Species == "Ophiodon elongatus", 1000,
                                                                                                                                                                                                 10)))))
                                                                                                                                                       )))))))))))))))))))))
  )
##### NEED TO FIGURE OUT BETTER MASS ESTIMATES!!!!!!!
# Jasmin says to use dry or ash free
# Rick has published bioass stuff

# weight calc for kelp greenling = log W = log a + b·log L
10^ ((-5.138) + 3.274*log10(200))

# weight calc for rockfish = log W = log a + b·log L
10^ ((-5.168) + 3.375*log10(690))

# blackeye goby
10^ ((-4.356) + 2.720*log10(75))

# Calculate how much mass each species contributes
rls_heaven <- rls_hell %>%
  mutate(total_mass = Total*mass)

# summarize total mass for each site
rls_summary <- rls_heaven %>%
  group_by(site_ID) %>%
  summarize(total_critter_mass = sum(total_mass)) %>%
  arrange(desc(total_critter_mass))

# what's biggest by weight? urchin then cukes
rls_biomass_ranked <- rls_heaven %>% 
  group_by(Species) %>%
  summarise(sum = sum(total_mass),
            total_animals = n())  %>%
  arrange(desc(total_animals))

#write_csv(rls_species, "Output/Output_data/RSL_abundant_species.csv")

# join mass and abundance data with ammonium data  
bottles_f2 <- bottles_f %>%
  left_join(rls_summary, by = "site_ID") %>%
  mutate(
    total_mass_kg = total_critter_mass/1000
  )

cuke_rank <- rls_heaven %>%
  filter(Species == "Apostichopus californicus") %>%
  group_by(site_name) %>%
  summarise(sum = sum(Total))  %>% arrange(desc(sum))


# more stats
model <- lmer(nh4_conc ~ depth + total_critter_mass + (1|site_ID), bottles_f2)
summary(model)
visreg(model)

# Graph mass vs pee????

ggplot(bottles_f2, aes(total_mass_kg, nh4_conc)) +
  geom_point(aes(colour = site), size = 2.5) +
  geom_smooth(method = lm) +
  labs(x = "Total animal biomass (kg)", y = "Ammonium concentration (umol)") +
  theme_black() 

#ggsave("Output/Figures/RLS_sites_pee_mass.png", device = "png",
#       height = 9, width = 16, dpi = 400)