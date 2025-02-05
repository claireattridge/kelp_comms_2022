## Barkley Sound community analyses
## Author: Claire Attridge
## Origin date: Feb 2023

# Loading necessary base packages
library(tidyverse)
library(ggpubr)
library(cowplot)
library(vegan)
library(MetBrewer)
library(ggforce)
library(grid)
library(gridExtra)
library(BiodiversityR) 
library(png)
library(ggrepel)
library(ecodist)
library(magick)

#

#### Loading & saving str & env data sheets ----

# Wave exp
rei <- read_csv("./MSc_data/Data_new/REI_2022.csv") %>%
  dplyr::select(-c(x,y)) %>%
  as.data.frame()

studysites <- as.vector(rei$SiteName) # Pulling a vector of study site names

# Temperature
temps <- read_csv("./MSc_data/Data_new/temps_2022.csv") %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  ungroup() %>%
  filter(SiteName %in% studysites) %>%
  droplevels() %>%
  dplyr::select(-SD_Tempave)

# Substrate
subs <- read_csv("./MSc_data/Data_new/Substrates_2022.csv") %>%
  dplyr::select(SiteName, Punderstory, Pcanopy, Phardbottom, Psoftbottom, Panimal, Pturf) %>%
  as.data.frame()

# Kelp metrics 
kelps <- read_csv("./MSc_data/Data_new/kelp_metrics_2022.csv") %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  dplyr::select(-c(HeightSD, BiomassSD, DensitySD, MacroSD, NereoSD))
kelps <- kelps %>%
  mutate(MacroP = MacroM/DensityM) # Creating a column for proportion Macro at a given site


# Listing the variables of interest together
dfs1 <- list(kelps, temps, rei, subs)

# Joining the predictor variables into one dataframe
preds <- reduce(dfs1, dplyr::left_join, by="SiteName")
preds[is.na(preds)] <- 0 # Translating any NAs to zeros



# Adding dominant kelp column (category) to the frame
preds_ord <- preds %>%
  mutate(Kelpdom = case_when(Composition == "Macro" | Composition == "Mixed" ~ "Giant", # Macro is dominant at mixed sites (too few mixed sites to statistically compare)
                             Composition == "Nereo" ~ "Bull",
                             Composition == "None" ~ "Giant")) %>% # Re. Starko (2022) - used to have M. pyrifera present around the area of Less Dangerous bay
  mutate(Kelpdom = as.factor(Kelpdom)) %>%
  droplevels()


###
## Removing second beach south as an outlier site 
preds_ord <- preds_ord %>%
  filter(SiteName != "Second Beach South") %>%
  droplevels()
## Removing the 'no kelp sites'
preds_ord <- preds_ord %>%
  filter(SiteName != "Wizard Islet North" & SiteName != "Less Dangerous Bay") %>%
  droplevels()
###


## saving a .csv file of the ecological and environmental predictors
write.csv(preds_ord, "./MSc_data/Data_new/AllPredictors_2022.csv", row.names=F)



# Removing the less important predictors
# Scaling all remaining variables
preds_ord_scale <- preds_ord %>%
  dplyr::select(-c(MacroM, NereoM, MacroP, ID, Depth_logger_m, Depth_logger_datum_m)) %>%
  ungroup() %>%
  mutate(across(c(2:5,7:17), scale))


## saving a .csv file of the scaled ecological and environmental predictors
write.csv(preds_ord_scale, "./MSc_data/Data_new/AllPredictors_scaled_2022.csv", row.names=F)


#


#### Community data: Cleaning ----


# Cleaning the RLS 'comms' frame
comms <- read_csv("./MSc_data/Data_new/RLS_2022_KDC_CMA_final.csv") %>%
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

# Cleaning up the frame columns
comms <- comms %>%
  mutate(Species = as.factor(Species), common_name = as.factor(common_name), SiteName = as.factor(site_name),abundance = as.numeric(Total), Date = as.POSIXct(Date, format="%d/%m/%Y")) %>% # Changing column formats
  dplyr::select(-c(site_name, Total)) # Removing unnecessary columns



## Filtering out the site redo surveys at Ross Islet 2, Less Dangerous Bay, & Second Beach S in Sept 2022
# and the acoustic baseline acoustic baseline site 'Sand Town'
comms <- comms %>%
  filter(Date <= as.POSIXct("2022-09-08") & SiteName != "Sand Town") %>%
  droplevels()

## Removing second beach south as an outlier site
comms <- comms %>%
  filter(SiteName != "Second Beach South") %>%
  droplevels()

## Removing the 'no kelp sites'
comms <- comms %>%
  filter(SiteName != "Wizard Islet North" & SiteName != "Less Dangerous Bay") %>%
  droplevels()



## Grouping Henricia spp. and Henricia leviuscula to a combined Henricia spp. as we can't be confident in the separation of these groups
comms_cleaned <- comms %>%
  mutate(Species = fct_collapse(Species, 'Henricia spp.'= c("Henricia spp.", "Henricia leviuscula")),
         common_name = fct_collapse(common_name, 'Blood star' = c("Pacific blood star", "Unidentified blood star"))) %>%
  droplevels()



## Exploring some general questions:

# Which site has the most species?
richness <- comms_cleaned %>%
  group_by(SiteName) %>%
  summarize(spp_rich = n_distinct(Species))%>%
  arrange(desc(spp_rich))

# Which species were most abundant overall?
abun <- comms_cleaned %>%
  group_by(Species) %>%
  summarize(totalab = sum(abundance)) %>%
  arrange(desc(totalab))

# Which species were most abundant on average per site?
abun_ave <- comms_cleaned %>%
  group_by(Species) %>%
  summarize(siteab = mean(abundance)) %>%
  arrange(desc(siteab))

# Which species were observed at the most sites?
obvs <- comms_cleaned %>%
  group_by(Species) %>%
  summarise(observs = n_distinct(SiteName)) %>%
  arrange(desc(observs))# Arranging in descending order

# Which species were observed at the most sites by forest type?
obvs_canopy <- comms_cleaned %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  mutate(canopy = case_when(SiteName == "Swiss Boy" ~ "Nereo", 
                              SiteName == "Second Beach" ~ "Nereo", 
                              SiteName == "Cable Beach (Blow Hole)" ~ "Nereo",
                              SiteName == "Bordelais Island" ~ "Nereo",
                            .default = "Macro")) %>%
  mutate(canopy = factor(canopy, levels=c("Nereo", "Macro"))) %>%
  group_by(Species, canopy) %>%
  summarise(observs = n_distinct(SiteName))



# Plot for species observations - Supp Fig
tiff(file="./MSc_plots/SuppFigs/SpeciesObvsCutoff.tiff", height = 7, width = 6, units="in", res=400)

obvsplot <- ggplot(data=obvs, aes(x=observs, y=(fct_reorder(Species, desc(-observs))), fill=observs)) + 
  geom_bar(stat="identity", position="dodge", width=0.4) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3")) +
  theme_classic() +
  theme(legend.position="none",
    axis.text.x = element_text(color="black", size="8.6"),
    axis.text.y = element_text(color="black", size="6", face="italic"),
    axis.title.y = element_text(color="black", size="9", vjust=0.5),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"), 
    ) +
  xlab("Site observations") + ylab("") +
  geom_hline(yintercept = 27.5, linetype=2) # Adding cutoff line for spp present at >= 3 sites
obvsplot

dev.off()



# Plot for species observations by canopy type - Supp Fig
tiff(file="./MSc_plots/SuppFigs/SpeciesObvsCutoff_canopy.tiff", height = 7, width = 6, units="in", res=400)

obvsplot_canopy <- ggplot(data=obvs_canopy, aes(x=observs, y=(reorder(Species, desc(-observs), sum)), fill=canopy)) + 
  geom_bar(stat="identity", position="stack", width=0.4) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL,
    labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  theme_classic() +
  theme(legend.position="top",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(color="black", size="7"),
        axis.text.x = element_text(color="black", size="8"),
        axis.text.y = element_text(color="black", size="6", face="italic"),
        axis.title.x = element_text(color="black", size="8.5", vjust=0.5),
        plot.margin = unit(c(0.5,1,0.5,0), "cm"), 
  ) +
  xlab("Site observations") + ylab("") +
  geom_hline(yintercept = 27.5, linetype=2) + # Adding cutoff line for spp present at >= 4 sites
  geom_hline(yintercept = 44.5, linetype=2, color="grey75") # Showing cutoff line for spp present at >= 5 sites
obvsplot_canopy

dev.off()



# Plot for species abundances by observation - Supp Fig
tiff(file="./MSc_plots/SuppFigs/SpeciesAbunCutoff.tiff", height = 7, width = 6, units="in", res=400)

# Joining the abundance and observation frames
abunobvs <- obvs %>%
  left_join(abun_ave, by="Species")

abun_obvsplot <- ggplot(data=abunobvs, aes(x=siteab, y=(fct_reorder(Species, desc(-observs))), fill=siteab)) + 
  geom_bar(stat="identity", position="dodge", width=0.4) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3")) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(color="black", size="8.6"),
        axis.text.y = element_text(color="black", size="6", face="italic"),
        axis.title.y = element_text(color="black", size="9", vjust=0.5),
        plot.margin = unit(c(0.5,1,0.5,0), "cm"), 
  ) +
  xlab("Average site abundance") + ylab("") +
  geom_hline(yintercept = 27.5, linetype=2) # Adding cutoff line for spp present at < 5 sites
abun_obvsplot

dev.off()



## We remove the rarely sighted species (present at >= 3 sites) from the analysis frame
# 88 species total, 27 species excluded
sppkeep <- obvs %>%
  filter(observs >= 3)

# Filtering the main frame for the retained species (n = 61)
commsclean <- comms %>%
  filter(Species %in% sppkeep$Species)


## What are the total abundance counts for all remaining species of
# interest across all sites?
comms_all <- commsclean %>%
  as.data.frame() %>%
  ungroup() %>%
  group_by(SiteName, Species, common_name, Method) %>%
  summarise(abundance = sum(abundance))



# saving a .csv file of the clean total taxa abundances by site
write.csv(comms_all, "./MSc_data/Data_new/RLS_2022.csv", row.names=F)

#


#### Community data: Site/Species groupings ----


# Loading the cleaned RLS spp abundances by site sheet
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  dplyr::select(SiteName, Species, common_name, abundance, Method) %>%
  rename(TaxaAb = abundance) %>%
  rename(CommonName = common_name) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName), CommonName = as.factor(CommonName), Species = as.factor(Species))


### Isolating fishes that utilize the water column (RLS Method 1 obvs)
taxap <- taxa %>%
  filter(Method == 1) %>%
  filter(CommonName != "Blackeye goby") %>% # Removing mislabelled obvs of benthic fish
  filter(CommonName != "Longfin sculpin") %>% # Removing mislabelled obvs of benthic fish
  filter(CommonName != "Red Irish lord") %>% # Removing mislabelled obvs of benthic fish
  droplevels() %>%
  select(-Method)
# 14 species 


### Isolating invertebrates (RLS Method 2 obvs - any demersal fishes)
M2benthicfish <- c("Snubnose sculpin", "Blackeye goby", "Crescent gunnel", "Penpoint gunnel", "Longfin sculpin", "Scalyhead sculpin", "Smoothhead sculpin", "Red Irish lord", "Buffalo sculpin")

taxab_inv <- taxa %>%
  filter(Method == 2) %>%
  filter(!CommonName %in% M2benthicfish) %>%
  droplevels()
# Removing mislabelled M1 fishes
taxab_inv <- taxab_inv %>%
  filter(CommonName != "Whitespotted greenling" & 
           CommonName != "Kelp greenling" & 
              CommonName != "Painted greenling") %>%
  droplevels() %>%
  select(-Method)
# 38 species


### Isolating fishes that remain on benthos (RLS Method 2)
taxab_fish <- taxa %>%
  filter(Method == 2) %>%
  filter(CommonName %in% M2benthicfish) %>%
  droplevels() %>%
  select(-Method)
# 9 species


#### Community data: Relative abundance tables & transformations ----


## WATER COLUMN FISH SPP
pfish_wide <- taxap %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>% # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(pfish_wide) <- pfish_wide$SiteName # Setting col as rownames
pfish_wide <- pfish_wide[,-1] # Removing the col used above
pfish_wide[is.na(pfish_wide)] <- 0 # NAs to zeros for calculations
pfish_mat <- as.matrix(pfish_wide) # As matrix

# pfish_rel <- make_relative(pfish_mat) # Relative abundance matrix of community data
pfish_hel <- decostand(pfish_mat, method = "hellinger") # Hellinger transformation of data


## BENTHIC FISH SPP
bfish_wide <- taxab_fish %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>%  # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(bfish_wide) <- bfish_wide$SiteName # Setting col as rownames
bfish_wide <- bfish_wide[,-1] # Removing the col used above
bfish_wide[is.na(bfish_wide)] <- 0 # NAs to zeros for calculations
bfish_mat <- as.matrix(bfish_wide) # As matrix

# bfish_rel <- make_relative(bfish_mat) # Relative abundance matrix of community data
bfish_hel <- decostand(bfish_mat, method = "hellinger") # Hellinger transformation of


## BENTHIC INVERT SPP
binv_wide <- taxab_inv %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>% # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(binv_wide) <- binv_wide$SiteName # Setting col as rownames
binv_wide <- binv_wide[,-1] # Removing the col used above
binv_wide[is.na(binv_wide)] <- 0 # NAs to zeros for calculations
binv_mat <- as.matrix(binv_wide) # As matrix

## binv_rel <- make_relative(binv_mat) # Relative abundance matrix of community data
binv_hel <- decostand(binv_mat, method = "hellinger") # Hellinger transformation of data

#


#### Community data: Abundance & rel abundance plots ----

## Cleaning raw table for pelagic fish
taxap_fish_frame <- taxap %>%
  ungroup() %>%
  as.data.frame()
## Grouping raw table for pelagic fish
taxap_fish_grp <- taxap %>%
  group_by(Species, CommonName) %>%
  summarise(TaxaAb_av = mean(TaxaAb), SD = sd(TaxaAb), Sites = n())


colourCount = length(unique(taxap_fish_frame$Species))
getPalette = colorRampPalette(met.brewer(name="Thomas", n=14))


# Making the regional abundance plot for pelagic fish
abPfish <- ggplot() +
  geom_jitter(data=taxap_fish_frame, position=position_jitter(width=0.3, height=0.2),
              aes(x=TaxaAb, y=as.factor(reorder(Species, TaxaAb, FUN=mean)), 
                  fill=as.factor(Species), alpha=0.2), color="grey55",
                  show.legend=FALSE, shape=21, size=1) +
  geom_pointrange(data=taxap_fish_grp, size=0.3, shape=22,
                  aes(y=as.factor(reorder(Species, TaxaAb_av)), x=TaxaAb_av, 
                  xmin=TaxaAb_av-SD, xmax=TaxaAb_av+SD,
                  color=as.factor(Species), fill=as.factor(Species))) +
  scale_fill_manual(values=getPalette(colourCount)) +
  scale_colour_manual(values=getPalette(colourCount)) +
  theme_classic() +
  theme(legend.position = "NA",
        strip.text.x = element_text(size=9, color="black", face="bold"),
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=9, color="black", face="italic"),
        axis.title.y = element_blank(),
        plot.margin = margin (5.5,20,5.5,5.5, unit="pt")
  ) +
  annotate("text", label="a", x=-50, y=14.2, fontface="bold", size=4.5)


## Grouping for number of site occurrences for each species
taxap_fish_occ <- taxap_fish_frame %>%
  dplyr::select(SiteName, Species, TaxaAb) %>%
  group_by(Species) %>%
  summarise(Mean=mean(TaxaAb), SiteNum=n()) %>%
  ungroup() %>%
  mutate(Prop = round((SiteNum/20)*100, digits=1))

# Making the prop. of occurrence plot for pelagic fish
ocPfish <- ggplot(taxap_fish_occ, aes(x=Species, y=Prop)) +
  geom_segment( aes(x=reorder(Species, Prop), xend=Species, y=0, yend=Prop), color="#453781FF", linewidth=1) +
  geom_point(color="#453781FF", size=3) +
  theme_classic() +
  scale_size_continuous("Mean Ab", breaks = c(100,200,400,600)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black", angle=65, hjust=1, size=9, face="italic"),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_blank()
  ) +
  annotate("text", label="a", x=1, y=100, fontface="bold", size=6)




## Cleaning raw table for benthic fish
taxab_fish_frame <- taxab_fish %>%
  ungroup() %>%
  as.data.frame()
## Grouping raw table for benthic fish
taxab_fish_grp <- taxab_fish %>%
  group_by(Species, CommonName) %>%
  summarise(TaxaAb_av = mean(TaxaAb), SD = sd(TaxaAb), Sites = n())

colourCount = length(unique(taxab_fish_frame$Species))
getPalette = colorRampPalette(met.brewer(name="Nattier", n=9))

# Making the regional abundance plot for demersal fish
abBfish <- ggplot() +
  geom_jitter(data=taxab_fish_frame, position=position_jitter(width=0.3, height=0.2),
              aes(x=TaxaAb, y=as.factor(reorder(Species, TaxaAb, FUN=mean)), 
                  fill=as.factor(Species), alpha=0.2), color="grey55",
              show.legend=FALSE, shape=21, size=1) +
  geom_pointrange(data=taxab_fish_grp, size=0.3, shape=22,
                  aes(y=as.factor(reorder(Species, TaxaAb_av)), x=TaxaAb_av, 
                      xmin=TaxaAb_av-SD, xmax=TaxaAb_av+SD,
                      color=as.factor(Species), fill=as.factor(Species))) +
  scale_fill_manual(values=getPalette(colourCount)) +
  scale_colour_manual(values=getPalette(colourCount)) +
  theme_classic() +
  theme(legend.position = "NA",
        strip.text.x = element_text(size=9, color="black", face="bold"),
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=9, color="black", face="italic"),
        axis.title.y = element_blank(),
        plot.margin = margin (5.5,20,5.5,5.5, unit="pt")
  ) +
  annotate("text", label="b", x=-9, y=9.15, fontface="bold", size=4.5)


## Grouping for number of site occurrences for each species
taxab_fish_occ <- taxab_fish_frame %>%
  dplyr::select(SiteName, Species, TaxaAb) %>%
  group_by(Species) %>%
  summarise(Mean=mean(TaxaAb), SiteNum=n()) %>%
  ungroup() %>%
  mutate(Prop = round((SiteNum/20)*100, digits=1))


# Making the prop. of occurrence plot for demersal fish
ocBfish <- ggplot(taxab_fish_occ, aes(x=Species, y=Prop)) +
  geom_segment( aes(x=reorder(Species, Prop), xend=Species, y=0, yend=Prop), color="#1F9A8AFF", linewidth=1) +
  geom_point(color="#1F9A8AFF", size=3) +
  theme_classic() +
  scale_size_continuous("Mean Ab", breaks = c(100,200,400,600), ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black", angle=65, hjust=1, size=9, face="italic"),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_blank()
  ) +
  annotate("text", label="b", x=1, y=100, fontface="bold", size=6)



## Grouping raw table for benthic inverts
taxab_inv_frame <- taxab_inv %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(split = factor(ifelse(Species == "Strongylocentrotus purpuratus" | 
                               Species == "Nucella lamellosa" | 
                               Species == "Patiria miniata" | 
                               Species == "Pomaulax gibberosus" | 
                               Species == "Mesocentrotus franciscanus", "large", "small")))
## Grouping raw table for benthic inverts
taxab_inv_grp <- taxab_inv %>%
  mutate(split = factor(ifelse(Species == "Strongylocentrotus purpuratus" | 
                                 Species == "Nucella lamellosa" | 
                                 Species == "Patiria miniata" | 
                                 Species == "Pomaulax gibberosus" | 
                                 Species == "Mesocentrotus franciscanus", "large", "small"))) %>%
  group_by(Species, CommonName, split) %>%
  summarise(TaxaAb_av = mean(TaxaAb), SD = sd(TaxaAb), Sites = n())
  


colourCount = length(unique(taxab_inv_frame$Species))
getPalette = colorRampPalette(met.brewer(name="Redon", n=38))


# Making the regional abundance plot for benthic inverts
# Setting the highly abundant spp on separate axis to visualize the data more clearly
abBinv <- ggplot()+
geom_jitter(data=taxab_inv_frame, position=position_jitter(width=0.3, height=0.2),
            aes(x=TaxaAb, y=as.factor(reorder(Species, TaxaAb, FUN=mean)), 
            fill=as.factor(Species), alpha=0.2), color="grey55",
            show.legend=FALSE, shape=21, size=1) +
geom_pointrange(data=taxab_inv_grp, size=0.3, shape=22,
                aes(y=as.factor(reorder(Species, TaxaAb_av)), x=TaxaAb_av, 
                xmin=TaxaAb_av-SD, xmax=TaxaAb_av+SD,
                color=as.factor(Species), fill=as.factor(Species))) +
scale_fill_manual(values=getPalette(colourCount)) +
scale_colour_manual(values=getPalette(colourCount)) +
  theme_classic() +
  theme(legend.position="NA",
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=9, color="black", face="italic"),
        axis.title.y = element_blank(),
        # plot.margin = margin(r = 10),
        axis.title.x = element_blank(),
        strip.background = element_blank(), # removing the facet boxes
        strip.text = element_blank(), # removing the facet text
        plot.margin = margin (5.5,10,5.5,5.5, unit="pt")
  ) + 
  facet_col(~split, scales="free", space="free") +
  geom_text(data = data.frame(x = -130, y = 4.2, split = "large", label = "c"), 
            aes(x = x, y = y, label = label), fontface="bold", size = 4.5)


## Grouping for number of site occurrences for each species
taxab_inv_occ <- taxab_inv_frame %>%
  dplyr::select(SiteName, Species, TaxaAb) %>%
  group_by(Species) %>%
  summarise(Mean=mean(TaxaAb), SiteNum=n()) %>%
  ungroup() %>%
  mutate(Prop = round((SiteNum/20)*100, digits=1))

# Making the prop. of occurrence plot for demersal fish
ocBinv <- ggplot(taxab_inv_occ, aes(x=Species, y=Prop)) +
  geom_segment( aes(x=reorder(Species, Prop), xend=Species, y=0, yend=Prop), color="#31688EFF", linewidth=1) +
  geom_point(color="#31688EFF", size=3) +
  theme_classic() +
  scale_size_continuous("Mean Ab", breaks = c(100,200,400,600)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(l = 14),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black", angle=65, hjust=1, size=9, face="italic"),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_blank(),
    # plot.margin = margin (5.5,5.5,5.5,, unit="pt")
  ) +
  annotate("text", label="c", x=1.5, y=100, fontface="bold", size=6)



# Arranging the smaller plots together (ABUNDANCE)
abfish <- ggarrange(abPfish + rremove("axis.title"), abBfish + rremove("axis.title"), ncol=1, align="hv")
abfish

# Setting up the layout
lay <- rbind(c(1,1,1,2,2,2,2),
             c(1,1,1,2,2,2,2))
# Rectangles to help visualize
gs <- lapply(1:2, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)


# Saving as an image object
tiff(file="./MSc_plots/SuppFigs/Ab_plots.tiff", height=7.5, width=10.5, units="in", res=300)

# End ggplot arrangement
grid <- grid.arrange(abfish, abBinv,
                     layout_matrix = lay, bottom=textGrob("Abundance"))
dev.off()



# Arranging the smaller plots together (OCCURRENCE)
ocfish <- ggarrange(ocPfish, ocBfish, ncol=2, align="hv") +
                      theme(plot.margin = margin(l = 10))
ocfish


# Setting up the layout
lay <- rbind(c(1,1),
             c(1,1),
             c(1,1),
             c(2,2),
             c(2,2),
             c(2,2),
             c(2,2))
# Rectangles to help visualize
gs <- lapply(1:2, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)

# Saving as an image object
tiff(file="./MSc_plots/SuppFigs/Occurr_plots.tiff", height=9, width=10, units="in", res=400)

# End ggplot arrangement
grid <- grid.arrange(ocfish, ocBinv, 
                     layout_matrix = lay, left=textGrob("Occurrence (% of 20 sites)", 
                                                        rot=90, vjust=1, hjust=0.02))

dev.off()


#### Pelagic fish: Generating models, stats ----

## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_scaled_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))


### CAP (kelp spp)
pfish_1 <- capscale(pfish_hel ~ Kelpdom, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray")
anova.cca(pfish_1) # Model is sig
RsquareAdj(pfish_1)
# Adjusted R squared of this model is: 7.9%



## SIMPER for dissimilarity
# Running SIMPER to test between groups
pfish_simp <- with(preds_ord, simper(pfish_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
pfish_simp_df <- as.data.frame(pfish_simp$Giant_Bull)
sum(pfish_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.562663
# Partitioned among the 14 pelagic fish spp = 0.0442699 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
pfish_simp_filt <- pfish_simp_df %>%
  filter(average > (2*0.04019021))
pfish_simp_filt

pfish_simpSPP <- "Embiotoca lateralis"

write.csv(pfish_simp_df, "C:/Users/clair/Desktop/pelagicfishsimper.csv")


### RDA (kelp forest structure)

# Generating the kelp forest structure model
pfish_2 <- capscale(pfish_hel ~ Kelpdom+DensityM+Area_m2+BiomassM+HeightM, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray") 
RsquareAdj(pfish_2)
## The structural and environmental variables explain 36.6% of variation (constrained axes)
# Adjusted R squared of this model is: 13.9%
vif.cca(pfish_2) # low vif scores

# Exploring outcomes
anova.cca(pfish_2, permutations=999) # significance of model
anova.cca(pfish_2, by = "margin", permutations=999) # significance of predictor terms
anova.cca(pfish_2, by = "axis", permutations=999) # significance by model axes



### RDA (environmental vars)

# Generating the environmental model
pfish_3 <- capscale(pfish_hel ~ Kelpdom+Tempave+Depth_datum_m+exp_36+Phardbottom+Psoftbottom+Punderstory+Pturf, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray")
RsquareAdj(pfish_3)
## Only the environmental variables explain explain 46.9% of variation (constrained axes)
# Adjusted R squared of this model is: 8.3%
vif.cca(pfish_2)

# Exploring significance
anova.cca(pfish_3, permutations = 999) # significance of model
anova.cca(pfish_3, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(pfish_3, permutations = 999, by = "axis") # significance by model axes


#


#### Pelagic fish: Plots ----

# Pelagic Fish 1 = CAP model (by kelp spp)
# Pelagic Fish 2 = RDA model (by kelp forest structure)
# Pelagic Fish 3 = RDA model (by environmental vars)


### PLOT 1: CAP (kelp spp)

# Quick ordiplot
pfish1_ordplot <- ordiplot(pfish_1, type="text", display="all")
ordisymbol(pfish1_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(pfish_1, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of sites from the model
pfish1_sites.long <- sites.long(pfish1_ordplot, env.data=preds_ord)
head(pfish1_sites.long)
# Extracting the locations of centroids from the sites.long output
pfish1_centroids.long <- centroids.long(pfish1_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish1_centroids.long)

# Creating hulls around points
hull_cyl <- pfish1_sites.long %>%
  group_by(Kelpdom) %>%
  slice(chull(axis1, axis2))

# The CAP plot
pfish1_ordplot <- ggplot(pfish1_sites.long) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (12.8%)") +
  ylab("MDS 1 (23.5%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=pfish1_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2) +
  stat_chull(data=hull_cyl, geom="polygon", aes(x=axis1, y=axis2, fill=Kelpdom), alpha=0.5) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill=guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position=c(0.25,0.14),
        legend.text=element_text(size=8.5, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(linewidth=0.5),
        axis.ticks=element_line(size=0.5)) +
  annotate("text", x = -2.1, y = 2.2, label = "a", size=5, fontface=2)
pfish1_ordplot


## Saving plot

tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ord_CAP.tiff", height = 4, width = 5, units = "in", res=400)

pfish1_ordplot

dev.off()




### Plot 2: RDA (kelp forest structure)

# Quick ordiplot
pfish2_ordplot <- ordiplot(pfish_2, type="text", display="all")
ordisymbol(pfish2_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(pfish_2, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of species from the model
pfish2_species.long <- species.long(pfish2_ordplot)
head(pfish2_species.long)

# Extracting the locations of sites from the model
pfish2_sites.long <- sites.long(pfish2_ordplot, env.data=preds_ord)
head(pfish2_sites.long)

# Extracting the locations of centroids from the sites.long output
pfish2_centroids.long <- centroids.long(pfish2_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish2_centroids.long)

# Extracting locations for env biplot correlations
pfish2_scrs <- as.data.frame(scores(pfish2_ordplot, display = "biplot"))


# Recoding fct level names
pfish2_scrs <- pfish2_scrs %>%
  mutate(vector = rownames(pfish2_scrs)) %>%
  mutate(vector = fct_recode(vector,
                             "Density" = "DensityM",
                             "Area" = "Area_m2",
                             "Height" = "HeightM",
                             "Biomass" = "BiomassM"))



# Extracting the 95 CI ellipses from the model and ordiplot
pfish2_ordsimple <- ordiplot(pfish_2)
pfish2_ordellipses <- with(preds_ord, ordiellipse(pfish2_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
pfish2_ellipses <- ordiellipse.long(pfish2_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
pfish2_species.envfit <- envfit(pfish2_ordsimple, env=pfish_hel) # For Bray distance ordiplot of Hellinger transformed spp data
pfish2_species.envfit.table <- data.frame(r=pfish2_species.envfit$vectors$r, p=pfish2_species.envfit$vectors$pvals)
pfish2_species.long.var <- species.long(pfish2_ordsimple, spec.data=pfish2_species.envfit.table)
pfish2_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
pfish2_species.long.vargrt <- pfish2_species.long.var[pfish2_species.long.var$r >= 0.5, ]
pfish2_species.long.vargrt


# Recoding spp fct levels to be shorthand scientific names
pfish2_species.long.vargrt <- pfish2_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "Shiner perch" = "Cymatogaster aggregata",
                             "Striped surfperch" = "Embiotoca lateralis",
                             "Rockfish spp." = "Sebastes spp."))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.02, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = pfish2_scrs %>% 
  mutate(r = sqrt(CAP1^2 + CAP2^2),
         theta = atan2(CAP2,CAP1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))



# The full structural plot 
pfish2_ordplot_str <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (15.5%)") +
  ylab("CAP 2 (12.8%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=pfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish2_scrs, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=CAP1*2.1, yend=CAP2*2.1),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2.1, y=ynew*2.1, label=vector),
            colour="black", size=2) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        # legend.position=c(0.12,0.93),
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(linewidth=0.3),
        axis.ticks=element_line(linewidth=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "a", size=5, fontface=2)
pfish2_ordplot_str




## The biplot for influential spp

# Radial shift function for arrow text
rshift = function(r, theta, a=0.02, b=0.09) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
spp.arrows = pfish2_species.long.vargrt %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full species plot
pfish2_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (15.5%)") +
  ylab("CAP 2 (12.8%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=pfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish2_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*2.4, yend=axis2*2.4),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=spp.arrows, # Adding in the SIMPER species labels
                  aes(x=xnew*2.4, y=ynew*2.4, label=labels),
                  colour="black", size=2, nudge_y = 0.09, nudge_x = -0.3) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "b", size=5, fontface=2)
pfish2_ordplot_spp



## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ord_RDAstr_draft.tiff", height = 3, width = 7, units = "in", res=400)

pfish2_ordplots <- ggarrange(pfish2_ordplot_str, pfish2_ordplot_spp, ncol=2, common.legend=FALSE)
pfish2_ordplots

dev.off()




### Plot 3: RDA (environmental vars)

# Quick ordiplot
pfish3_ordplot <- ordiplot(pfish_3, type="text", display="all")
ordisymbol(pfish3_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(pfish_3, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of species from the model
pfish3_species.long <- species.long(pfish3_ordplot)
head(pfish3_species.long)

# Extracting the locations of sites from the model
pfish3_sites.long <- sites.long(pfish3_ordplot, env.data=preds_ord)
head(pfish3_sites.long)

# Extracting the locations of centroids from the sites.long output
pfish3_centroids.long <- centroids.long(pfish3_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish3_centroids.long)

# Extracting locations for env biplot correlations
pfish3_env <- as.data.frame(scores(pfish3_ordplot, display = "biplot"))


# Recoding fct level names
pfish3_env <- pfish3_env %>%
  mutate(vector = rownames(pfish3_env)) %>%
  mutate(vector = fct_recode(vector,
                             "Temp" = "Tempave",
                             "Depth" = "Depth_datum_m",
                             "Wave Exp" = "exp_36",
                             "Understory %" = "Punderstory",
                             "Hard %" = "Phardbottom",
                             "Soft %" = "Psoftbottom",
                             "Turf %" = "Pturf"))


# Extracting the ellipses from the model and ordiplot
pfish3_ordsimple <- ordiplot(pfish_3)
pfish3_ordellipses <- with(preds_ord, ordiellipse(pfish3_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
pfish3_ellipses <- ordiellipse.long(pfish3_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
pfish3_species.envfit <- envfit(pfish3_ordsimple, env=pfish_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
pfish3_species.envfit.table <- data.frame(r=pfish3_species.envfit$vectors$r, p=pfish3_species.envfit$vectors$pvals)
pfish3_species.long.var <- species.long(pfish3_ordsimple, spec.data=pfish3_species.envfit.table)
pfish3_species.long.var

# Selecting for individual spp that explain > 50% of the observed variation among groups
pfish3_species.long.vargrt <- pfish3_species.long.var[pfish3_species.long.var$r >= 0.5, ]
pfish3_species.long.vargrt


# Recoding spp fct levels to be shorthand scientific names
pfish3_species.long.vargrt <- pfish3_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "Shiner perch" = "Cymatogaster aggregata",
                             "Striped surfperch" = "Embiotoca lateralis",
                             "Rockfish spp." = "Sebastes spp."
                             ))



## The biplot for environmental vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.09, b=0.05) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = pfish3_env %>% 
  mutate(r = sqrt(CAP1^2 + CAP2^2),
         theta = atan2(CAP2,CAP1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
pfish3_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (18.8%)") +
  ylab("CAP 2 (10.0%)") +   
  scale_y_continuous(limits=c(-2.5,2.5)) + 
  scale_x_continuous(limits=c(-2.5,2.5)) + 
  geom_point(data=pfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish3_env, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=CAP1*2.9, yend=CAP2*2.9),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2.9, y=ynew*2.9, label=vector),
            colour="black", size=2) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        # legend.position=c(0.14,0.93),
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "c", size=5, fontface=2)
pfish3_ordplot_env



## The biplot for influential spp

# Calculate shift of text from arrows
spp.arrows = pfish3_species.long.vargrt %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

# The full species plot
pfish3_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (18.8%)") +
  ylab("CAP 2 (10.0%)") +   
  scale_y_continuous(limits=c(-2.5,2.5)) + 
  scale_x_continuous(limits=c(-2.5,2.5)) + 
  geom_point(data=pfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish3_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*1.2, yend=axis2*1.2),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=pfish3_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*1.8, y=axis2*1.4, label=labels),
                  colour="black", size=2, angle=0, nudge_y=-0.08)+
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "d", size=5, fontface=2)
pfish3_ordplot_spp


## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ord_RDAenv_draft.tiff", height = 3, width = 7, units = "in", res=400)

pfish3_ordplots <- ggarrange(pfish3_ordplot_env, pfish3_ordplot_spp, ncol=2, common.legend=FALSE)
pfish3_ordplots

dev.off()

#


### FINAL ARRANGEMENT

# Setting up the layout
lay <- rbind(c(3,3,3), # Leaving room for images
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2))
# Rectangles to help visualize
gs <- lapply(1:3, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)

# Saving as an image object
# pdf(file="./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined.pdf", height=5.5, width=10)
tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined.tiff", height=15, width=17, units="cm", res=300)

# End ggplot arrangement
grid <- grid.arrange(pfish2_ordplots, pfish3_ordplots,
             layout_matrix = lay)

dev.off()


# Drawing in the animal outlines

# Loading spp image outlines
rfish <- image_read('./Phylopic/Fish/rockfish.png')
sperch <- image_read('./Phylopic/Fish/surfperch.png')

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined.tiff") 

# Adding in the animal outlines (scaling & location)
mmap <- image_composite(mplot, image_scale(rfish, "x90"), offset="+770+60")
mmap <- image_composite(mmap, image_scale(sperch, "x100"), offset="+1070+40")

# Saving the final fig
image_write(mmap, path = "./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined_outlines.tiff", format = "tiff", density=400)

#

#### Benthic fish: Generating ordination models, stats ----


## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_scaled_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))



### CAP (kelp spp)
bfish_1 <- capscale(bfish_hel ~ Kelpdom, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray")
anova.cca(bfish_1) # model is sig
RsquareAdj(bfish_1)
# Adjusted R squared of this model is: 27.4%



## SIMPER for dissimilarity
# Running SIMPER to test between groups (limiting to 10 most influential)
bfish_simp <- with(preds_ord, simper(bfish_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
bfish_simp_df <- as.data.frame(bfish_simp$Giant_Bull)
sum(bfish_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.6741229
# Partitioned among the 9 demersal fish spp = 0.07490254 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
bfish_simp_filt <- bfish_simp_df %>%
  filter(average > (2*0.07490254))
bfish_simp_filt

bfish_simpSPP <- c("Rhinogobiops nicholsii")


write.csv(bfish_simp_df, "C:/Users/clair/Desktop/demersalfishsimper.csv")


### RDA (kelp forest structure)

# Generating the kelp forest structure model
bfish_2 <- capscale(bfish_hel ~ Kelpdom+DensityM+Area_m2+BiomassM+HeightM, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray") 
RsquareAdj(bfish_2)
## Only the environmental variables explain 45.9% of variation (constrained axes)
# Adjusted R squared of this model is: 26.6%
vif.cca(bfish_2) # low vif scores

# Exploring outcomes
anova.cca(bfish_2) # significance of model
anova.cca(bfish_2, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(bfish_2, permutations = 999, by = "axis") # significance by model axes



### RDA (environmental vars)

# Generating the environmental model
bfish_3 <- capscale(bfish_hel ~ Kelpdom+Tempave+Depth_datum_m+exp_36+Phardbottom+Psoftbottom+Punderstory+Pturf, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray")
RsquareAdj(bfish_3)
## The structural and environmental variables explain explain 56.6% of variation (constrained axes)
# Adjusted R squared of this model is: 25.0%
vif.cca(bfish_3) # low vif scores

# Exploring model significance
anova.cca(bfish_3) # significance of model
anova.cca(bfish_3, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(bfish_3, permutations = 999, by = "axis") # significance by model axes



# 
# ### Checking for more specific correlations
# 
# library(psych)
# 
# # merging the demersal fish comm dataset with the predictors
# bfish_corrcheck <- cbind(bfish_wide, preds_ord)
# 
# 
# # checking/visualizing correlations
# # correlation coefficients in the upper right hand (size scaled to their |r|)
# pairs.panels(bfish_corrcheck[,c(1:10,12:15)], scale=T) 
# pairs.panels(bfish_corrcheck[,c(1:10,21:23)], scale=T) 
# pairs.panels(bfish_corrcheck[,c(1:10,28:34)], scale=T) 
# 
# 
# #


#### Benthic fish: Plots ----


# Pelagic Fish 1 = CAP model (by kelp spp)
# Pelagic Fish 2 = RDA model (by kelp forest structure)
# Pelagic Fish 3 = RDA model (by environmental vars)


### PLOT 1: CAP (kelp spp)

# Quick ordiplot
bfish1_ordplot <- ordiplot(bfish_1, type="text", display="all")
ordisymbol(bfish1_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(bfish_1, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of sites from the model
bfish1_sites.long <- sites.long(bfish1_ordplot, env.data=preds_ord)
head(bfish1_sites.long)

# Extracting the locations of centroids from the sites.long output
bfish1_centroids.long <- centroids.long(bfish1_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(bfish1_centroids.long)

# Creating hulls around the points
hull_cyl <- bfish1_sites.long %>%
  group_by(Kelpdom) %>%
  slice(chull(axis1, axis2))


# The CAP plot
bfish1_ordplot <- ggplot(bfish1_sites.long) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (31.2%)") +
  ylab("MDS 1 (36.4%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=bfish1_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2) +
  stat_chull(data=hull_cyl, geom="polygon", aes(x=axis1, y=axis2, fill=Kelpdom), alpha=0.5) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill=guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position=c(0.25,0.14),
        legend.text=element_text(size=8.5, color="black"),
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(linewidth=0.5),
        axis.ticks=element_line(linewidth=0.5)) +
  annotate("text", x = -2.1, y = 2.2, label = "b", size=5, fontface=2)
bfish1_ordplot


## Saving plot

tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ord_CAP.tiff", height = 4, width = 5, units = "in", res=400)

bfish1_ordplot

dev.off()




### Plot 2: RDA (kelp forest structure)

# Quick ordiplot
bfish2_ordplot <- ordiplot(bfish_2, type="text", display="all")
ordisymbol(bfish2_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(bfish_2, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of species from the model
bfish2_species.long <- species.long(bfish2_ordplot)
head(bfish2_species.long)

# Extracting the locations of sites from the model
bfish2_sites.long <- sites.long(bfish2_ordplot, env.data=preds_ord)
head(bfish2_sites.long)

# Extracting the locations of centroids from the sites.long output
bfish2_centroids.long <- centroids.long(bfish2_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(bfish2_centroids.long)

# Extracting locations for env biplot correlations
bfish2_scrs <- as.data.frame(scores(bfish2_ordplot, display = "biplot"))


# Recoding fct level names
bfish2_scrs <- bfish2_scrs %>%
  mutate(vector = rownames(bfish2_scrs)) %>%
  mutate(vector = fct_recode(vector,
                             "Density" = "DensityM",
                             "Area" = "Area_m2",
                             "Height" = "HeightM",
                             "Biomass" = "BiomassM"))


# Extracting the ellipses from the model and ordiplot
bfish2_ordsimple <- ordiplot(bfish_2)
bfish2_ordellipses <- with(preds_ord, ordiellipse(bfish2_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
bfish2_ellipses <- ordiellipse.long(bfish2_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
bfish2_species.envfit <- envfit(bfish2_ordsimple, env=bfish_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
bfish2_species.envfit.table <- data.frame(r=bfish2_species.envfit$vectors$r, p=bfish2_species.envfit$vectors$pvals)
bfish2_species.long.var <- species.long(bfish2_ordsimple, spec.data=bfish2_species.envfit.table)
bfish2_species.long.var

# Selecting for individual spp that explain > 50% of the observed variation among groups
bfish2_species.long.vargrt <- bfish2_species.long.var[bfish2_species.long.var$r >= 0.5, ]
bfish2_species.long.vargrt


# Recoding spp fct levels to be shorthand scientific names
bfish2_species.long.vargrt <- bfish2_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "Scalyhead sculpin" = "Artedius harringtoni",
                             "Smoothhead sculpin" = "Artedius lateralis",
                             "Longfin sculpin" = "Jordania zonope",
                             "Blackeye goby" = "Rhinogobiops nicholsii"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.08, b=0.1) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = bfish2_scrs %>% 
  mutate(r = sqrt(CAP1^2 + CAP2^2),
         theta = atan2(CAP2,CAP1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))



# The full structural plot 
bfish2_ordplot_str <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (36.6%)") + # CAP 1
  ylab("CAP 2 (2.9%)") + # CAP 2
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=bfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish2_scrs, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=CAP1*2.8, yend=CAP2*2.8),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2.8, y=ynew*2.8, label=vector),
            colour="black", size=2) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        # legend.position=c(0.12,0.93),
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "a", size=5, fontface=2)

bfish2_ordplot_str



## The biplot for influential spp

# Radial shift function for arrow text
rshift = function(r, theta, a=0.08, b=0.04) {
  r + a + b*abs(cos(theta))
}


# Calculate shift of text from arrows
spp.arrows = bfish2_species.long.vargrt %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full species plot
bfish2_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (36.6%)") + # CAP 1
  ylab("CAP 2 (2.9%)") + # CAP 2
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=bfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish2_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=spp.arrows, # Adding in the SIMPER species labels
            aes(x=xnew*2.3, y=ynew*2.1, label=labels),
            colour="black", size=2, nudge_x=0.36, nudge_y=-0.18) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "b", size=5, fontface=2)
bfish2_ordplot_spp


## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ord_RDAstr_draft.tiff", height = 3.5, width = 7, units = "in", res=400)

bfish2_ordplots <- ggarrange(bfish2_ordplot_str, bfish2_ordplot_spp, ncol=2, common.legend=FALSE)
bfish2_ordplots

dev.off()




### Plot 3: RDA (environmental vars)

# Quick ordiplot
bfish3_ordplot <- ordiplot(bfish_3, type="text", display="all")
ordisymbol(bfish3_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(bfish_3, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of species from the model
bfish3_species.long <- species.long(bfish3_ordplot)
head(bfish3_species.long)

# Extracting the locations of sites from the model
bfish3_sites.long <- sites.long(bfish3_ordplot, env.data=preds_ord)
head(bfish3_sites.long)

# Extracting the locations of centroids from the sites.long output
bfish3_centroids.long <- centroids.long(bfish3_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(bfish3_centroids.long)

# Extracting locations for env biplot correlations
bfish3_env <- as.data.frame(scores(bfish3_ordplot, display = "biplot"))
bfish3_env


# Recoding fct level names
bfish3_env <- bfish3_env %>%
  mutate(vector = rownames(bfish3_env)) %>%
  mutate(vector = fct_recode(vector,
                             "Temp" = "Tempave",
                             "Depth" = "Depth_datum_m",
                             "Wave Exp" = "exp_36",
                             "Understory %" = "Punderstory",
                             "Hard %" = "Phardbottom",
                             "Soft %" = "Psoftbottom",
                             "Turf %" = "Pturf"))


# Extracting the ellipses from the model and ordiplot
bfish3_ordsimple <- ordiplot(bfish_3)
bfish3_ordellipses <- with(preds_ord, ordiellipse(bfish3_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
bfish3_ellipses <- ordiellipse.long(bfish3_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
bfish3_species.envfit <- envfit(bfish3_ordsimple, env=bfish_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
bfish3_species.envfit.table <- data.frame(r=bfish3_species.envfit$vectors$r, p=bfish3_species.envfit$vectors$pvals)
bfish3_species.long.var <- species.long(bfish3_ordsimple, spec.data=bfish3_species.envfit.table)
bfish3_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
bfish3_species.long.vargrt <- bfish3_species.long.var[bfish3_species.long.var$r >= 0.5, ]
bfish3_species.long.vargrt


# Recoding spp fct levels to be shorthand scientific names
bfish3_species.long.vargrt <- bfish3_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "Scalyhead sculpin" = "Artedius harringtoni",
                             "Smoothhead sculpin" = "Artedius lateralis",
                             "Longfin sculpin" = "Jordania zonope",
                             "Blackeye goby" = "Rhinogobiops nicholsii"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.26, b=0.28) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = bfish3_env %>% 
  mutate(r = sqrt(CAP1^2 + CAP2^2),
         theta = atan2(CAP2,CAP1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
bfish3_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (39.6%)") + # CAP 1
  ylab("CAP 2 (6.7%)") + # CAP 2
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  geom_point(data=bfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish3_env, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=CAP1*2.8, yend=CAP2*2.8),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2, y=ynew*2, label=vector),
            colour="black", size=2, hjust=0, vjust=0, box.padding=0.0415, nudge_x=-0.02, nudge_y=0.026) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        # legend.position=c(0.14,0.93),
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "c", size=5, fontface=2)
bfish3_ordplot_env



## The biplot for influential spp


# Radial shift function for arrow text
rshift = function(r, theta, a=0.13, b=0.13) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
spp.arrows = bfish3_species.long.vargrt %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

# The full species plot
bfish3_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (39.6%)") + # CAP 1
  ylab("CAP 2 (6.7%)") + # CAP 2  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  geom_point(data=bfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish3_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*1.6, yend=axis2*1.6),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=spp.arrows, # Adding in the SIMPER species labels
            aes(x=xnew*1.7, y=ynew*1.8, label=labels),
            colour="black", size=2, nudge_x=0.2) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "d", size=5, fontface=2)
bfish3_ordplot_spp


## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ord_RDAenv_draft.tiff", height = 3.5, width = 7, units = "in", res=400)

bfish3_ordplots <- ggarrange(bfish3_ordplot_env, bfish3_ordplot_spp, ncol=2, common.legend=FALSE)
bfish3_ordplots

dev.off()

#


### FINAL ARRANGEMENT

# Setting up the layout
lay <- rbind(c(3,3,3), # Leaving room for images
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2))
# Rectangles to help visualize
gs <- lapply(1:3, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)


tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ords_combined.tiff", height=15, width=17, units="cm", res=300)

# End ggplot arrangement
grid.arrange(bfish2_ordplots, bfish3_ordplots,
             layout_matrix = lay)

dev.off()

#

# Drawing in the animal outlines

# Loading spp image outlines
bgoby <- image_read('./Phylopic/Fish/blackeyegoby.png')
sculp <- image_read('./Phylopic/Fish/sculpin.png')

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/PaperFigs/Bfish/Bfish_ords_combined.tiff") 

# Adding in the animal outlines (scaling & location)
mmap <- image_composite(mplot, image_scale(bgoby, "x70"), offset="+780+60")
mmap <- image_composite(mmap, image_rotate(image_scale(sculp, "x70"), 340), offset="+1070+05")

# Saving the final fig
image_write(mmap, path = "./MSc_plots/PaperFigs/Bfish/Bfish_ords_combined_outlines.tiff", format = "tiff", density=400)

#


#### Benthic inverts: Generating ordination models, stats ----


## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_scaled_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))


### CAP (kelp spp)
binv_1 <- capscale(binv_hel ~ Kelpdom, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray") 
anova.cca(binv_1) # model is sig
RsquareAdj(binv_1)
# Adjusted R squared of this model is: 6.7%



## SIMPER for dissimilarity 
# Running SIMPER to test between groups (limiting to 10 most influential)
binv_simp <- with(preds_ord, simper(binv_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
binv_simp_df <- as.data.frame(binv_simp$Giant_Bull)
sum(binv_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.4679515
# Partitioned among the 38 invertebrate spp = 0.01231451 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
binv_simp_filt <- binv_simp_df %>%
  filter(average > (2*0.01231451))
binv_simp_filt

binv_simpSPP <- c("Pomaulax gibberosus", "Strongylocentrotus purpuratus",
                  "Mesocentrotus franciscanus", "Patiria miniata")


write.csv(binv_simp_df, "C:/Users/clair/Desktop/invertsimper.csv")



### RDA (kelp forest structure)

# Generating the kelp forest structure model
binv_2 <- capscale(binv_hel ~ Kelpdom+DensityM+Area_m2+BiomassM+HeightM, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray") 
RsquareAdj(binv_2)
## Only the structural variables explain 44.7% of variation (constrained axes)
# Adjusted R squared of this model is: 24.9%
vif.cca(binv_2) # low vif scores

# Exploring model significance
anova.cca(binv_2) # significance of model
anova.cca(binv_2, step = 1000, by = "margin") # significance of predictor terms
anova.cca(binv_2, step = 1000, by = "axis") # significance by model axes



### RDA (environmental vars)

# Generating the environmental model
binv_3 <- capscale(binv_hel ~ Kelpdom+Tempave+Depth_datum_m+exp_36+Phardbottom+Psoftbottom+Punderstory+Pturf, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray")
RsquareAdj(binv_3)
## The environmental variables explain explain 58.6% of variation (constrained axes)
# Adjusted R squared of this model is: 28.4%
vif.cca(binv_3)

# Exploring model significance
anova.cca(binv_3) # significance of model
anova.cca(binv_3, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(binv_3, permutations = 999, by = "axis") # significance by model axes


#


#### Benthic inverts: Plots ----


# Benthic Invert 1 = CAP model (by kelp spp)
# Benthic Invert 2 = RDA model (by kelp forest structure)
# Benthic Invert 3 = RDA model (by environmental vars)


### PLOT 1: CAP (kelp spp)

# Quick ordiplot
binv1_ordplot <- ordiplot(binv_1, type="text", display="all")
ordisymbol(binv1_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(binv_1, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of sites from the model
binv1_sites.long <- sites.long(binv1_ordplot, env.data=preds_ord)
head(binv1_sites.long)

# Extracting the locations of centroids from the sites.long output
binv1_centroids.long <- centroids.long(binv1_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(binv1_centroids.long)

# Creating hulls around the data points
hull_cyl <- binv1_sites.long %>%
  group_by(Kelpdom) %>%
  slice(chull(axis1, axis2))


# The CAP plot
binv1_ordplot <- ggplot(binv1_sites.long) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (11.7%)") +
  ylab("MDS 1 (27.9%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=binv1_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2) +
  stat_chull(data=hull_cyl, geom="polygon", aes(x=axis1, y=axis2, fill=Kelpdom), alpha=0.5) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill=guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position=c(0.25,0.14),
        legend.text=element_text(size=8.5, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(linewidth=0.5),
        axis.ticks=element_line(size=0.5)) +
  annotate("text", x = -2.1, y = 2.2, label = "c", size=5, fontface=2)
binv1_ordplot


## Saving plot

tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ord_CAP.tiff", height = 4, width = 5, units = "in", res=400)

binv1_ordplot

dev.off()




### Plot 2: RDA (kelp forest structure)

# Quick ordiplot
binv2_ordplot <- ordiplot(binv_2, type="text", display="all")
ordisymbol(binv2_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(binv_2, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of species from the model
binv2_species.long <- species.long(binv2_ordplot)
head(binv2_species.long)

# Extracting the locations of sites from the model
binv2_sites.long <- sites.long(binv2_ordplot, env.data=preds_ord)
head(binv2_sites.long)

# Extracting the locations of centroids from the sites.long output
binv2_centroids.long <- centroids.long(binv2_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(binv2_centroids.long)

# Extracting locations for env biplot correlations
binv2_scrs <- as.data.frame(scores(binv2_ordplot, display = "biplot"))


# Recoding fct level names
binv2_scrs <- binv2_scrs %>%
  mutate(vector = rownames(binv2_scrs)) %>%
  mutate(vector = fct_recode(vector,
                             "Density" = "DensityM",
                             "Area" = "Area_m2",
                             "Height" = "HeightM",
                             "Biomass" = "BiomassM"))


# Extracting the ellipses from the model and ordiplot
binv2_ordsimple <- ordiplot(binv_2)
binv2_ordellipses <- with(preds_ord, ordiellipse(binv2_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
binv2_ellipses <- ordiellipse.long(binv2_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
binv2_species.envfit <- envfit(binv2_ordsimple, env=binv_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
binv2_species.envfit.table <- data.frame(r=binv2_species.envfit$vectors$r, p=binv2_species.envfit$vectors$pvals)
binv2_species.long.var <- species.long(binv2_ordsimple, spec.data=binv2_species.envfit.table)
binv2_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
binv2_species.long.vargrt <- binv2_species.long.var[binv2_species.long.var$r >= 0.5, ]
binv2_species.long.vargrt


# Recoding spp fct levels to be shorthand scientific names
binv2_species.long.vargrt <- binv2_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "Red rock crab" = "Cancer productus",
                             "Leafy hornmouth" = "Ceratostoma foliatum",
                             "Leather star" = "Dermasterias imbricata",
                             "Leopard dorid" = "Diaulula odonoghuei",
                             "Mottled star" = "Evasterias troschelii",
                             "Red sea urchin" = "Mesocentrotus franciscanus",
                             "Graceful crab" = "Oregonia gracilis",
                             "Red turban snail" = "Pomaulax gibberosus",
                             "Sunflower star" = "Pycnopodia helianthoides"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = binv2_scrs %>% 
  mutate(r = sqrt(CAP1^2 + CAP2^2),
         theta = atan2(CAP2,CAP1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full structural plot 
binv2_ordplot_str <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (21.3%)") +
  ylab("CAP 2 (13.5%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=binv2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv2_scrs, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=CAP1*2, yend=CAP2*2),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2, y=ynew*2, label=vector),
            colour="black", size=2) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        # legend.position=c(0.12,0.93),
        legend.text=element_text(size=8.5, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "a", size=5, fontface=2)
binv2_ordplot_str



## The biplot for influential spp

# Calculate shift of text from arrows
spp.arrows = binv2_species.long.vargrt %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full species plot
binv2_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (21.3%)") +
  ylab("CAP 2 (13.5%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=binv2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv2_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*2.2, yend=axis2*2.2),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=binv2_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*2.2, y=axis2*2.2, label=labels),
                  colour="black",size=2, angle=0, hjust=0, vjust=0,
                  box.padding=0.44, nudge_y=0.84, nudge_x=0.42, segment.linetype=2,
                  segment.size=0.3, segment.color="grey40") +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=8.5, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "b", size=5, fontface=2)

binv2_ordplot_spp



## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ord_RDAstr_draft.tiff", height = 3.5, width = 7, units = "in", res=400)

binv2_ordplots <- ggarrange(binv2_ordplot_str, binv2_ordplot_spp, ncol=2, common.legend=FALSE)
binv2_ordplots

dev.off()




### Plot 3: RDA (environmental vars)

# Quick ordiplot
binv3_ordplot <- ordiplot(binv_3, type="text", display="all")
ordisymbol(binv3_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(binv_3, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the locations of species from the model
binv3_species.long <- species.long(binv3_ordplot)
head(binv3_species.long)

# Extracting the locations of sites from the model
binv3_sites.long <- sites.long(binv3_ordplot, env.data=preds_ord)
head(binv3_sites.long)

# Extracting the locations of centroids from the sites.long output
binv3_centroids.long <- centroids.long(binv3_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(binv3_centroids.long)

# Extracting locations for env biplot correlations
binv3_env <- as.data.frame(scores(binv3_ordplot, display = "biplot"))
binv3_env


# Recoding fct level names
binv3_env <- binv3_env %>%
  mutate(vector = rownames(binv3_env)) %>%
  mutate(vector = fct_recode(vector,
                             "Temp" = "Tempave",
                             "Depth" = "Depth_datum_m",
                             "Wave Exp" = "exp_36",
                             "Understory %" = "Punderstory",
                             "Hard %" = "Phardbottom",
                             "Soft %" = "Psoftbottom",
                             "Turf %" = "Pturf"))


# Extracting the ellipses from the model and ordiplot
binv3_ordsimple <- ordiplot(binv_3)
binv3_ordellipses <- with(preds_ord, ordiellipse(binv3_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
binv3_ellipses <- ordiellipse.long(binv3_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
binv3_species.envfit <- envfit(binv3_ordsimple, env=binv_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
binv3_species.envfit.table <- data.frame(r=binv3_species.envfit$vectors$r, p=binv3_species.envfit$vectors$pvals)
binv3_species.long.var <- species.long(binv3_ordsimple, spec.data=binv3_species.envfit.table)
binv3_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
binv3_species.long.vargrt <- binv3_species.long.var[binv3_species.long.var$r >= 0.5, ]
binv3_species.long.vargrt

# Recoding spp fct levels to be shorthand scientific names
binv3_species.long.vargrt <- binv3_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "Red rock crab" = "Cancer productus",
                             "Leafy hornmouth" = "Ceratostoma foliatum",
                             "Leather star" = "Dermasterias imbricata",
                             "Leopard dorid" = "Diaulula odonoghuei",
                             "Mottled star" = "Evasterias troschelii",
                             "Red sea urchin" = "Mesocentrotus franciscanus",
                             "Graceful crab" = "Oregonia gracilis",
                             "Red turban snail" = "Pomaulax gibberosus",
                             "Sunflower star" = "Pycnopodia helianthoides",
                             ))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.08, b=0.1) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = binv3_env %>% 
  mutate(r = sqrt(CAP1^2 + CAP2^2),
         theta = atan2(CAP2,CAP1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
binv3_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (21.2%)") +
  ylab("CAP 2 (19.9%)") +  
  scale_y_continuous(limits=c(-2.2,2.2)) + 
  scale_x_continuous(limits=c(-2.2,2.2)) + 
  geom_point(data=binv3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv3_env, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=CAP1*2, yend=CAP2*2),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2, y=ynew*2, label=vector),
            colour="black", size=2, nudge_y=-0.05) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=23))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        # legend.position=c(0.14,0.93),
        legend.text=element_text(size=8.5, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -2.2, y = 2.2, label = "c", size=5, fontface=2)
binv3_ordplot_env



## The biplot for influential spp

# Calculate shift of text from arrows
spp.arrows = binv3_species.long.vargrt %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

# The full species plot
binv3_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (21.2%)") +
  ylab("CAP 2 (19.9%)") +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=binv3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv3_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*3.5, yend=axis2*3.5),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=binv3_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*3.5, y=axis2*3.5, label=labels),
                  colour="black", fontface="italic", size=2,
                  box.padding=0.31, nudge_y=-1.2, nudge_x=1.55, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=10, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=9, color="black"),
        axis.title=element_text(size=9.5, color="black"),
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3.5, y = 3.5, label = "d", size=5, fontface=2)

binv3_ordplot_spp



## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ord_RDAenv_draft.tiff", height = 3.5, width = 7, units = "in", res=400)

binv3_ordplots <- ggarrange(binv3_ordplot_env, binv3_ordplot_spp, ncol=2, common.legend=FALSE)
binv3_ordplots

dev.off()

#


### FINAL ARRANGEMENT

# Setting up the layout
lay <- rbind(c(3,3,3), # Leaving room for images
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2),
             c(2,2,2))
# Rectangles to help visualize
gs <- lapply(1:3, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)


tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ords_combined.tiff", height=15, width=17, units="cm", res=300)

# End ggplot arrangement
grid.arrange(binv2_ordplots, binv3_ordplots,
             layout_matrix = lay)

dev.off()

#

# Drawing in the animal outlines

# Loading spp image outlines
purchin <- image_read('./Phylopic/EladTrace/purpleurchin.png')
turban <- image_read('./Phylopic/EladTrace/turbansnail.png')
bstar <- image_read('./Phylopic/EladTrace/bloodstar.png')

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/PaperFigs/Invert/Invert_ords_combined.tiff") 

# Adding in the animal outlines (scaling & location)
mmap <- image_composite(mplot, image_scale(purchin, "x200"), offset="+690+08")
mmap <- image_composite(mmap, image_scale(bstar, "x150"), offset="+900+05")
mmap <- image_composite(mmap, image_rotate(image_scale(turban, "x130"), 330), offset="+1090+05")

# Saving the final fig
image_write(mmap, path = "./MSc_plots/PaperFigs/Invert/Invert_ords_combined_outlines.tiff", format = "tiff", density=400)

#


#### Trends in overall abundance & richness ----

library(DHARMa)
library(glmmTMB)
library(car)
library(performance)
library(mgcv)
library(tidymv)
library(dotwhisker)

# Loading the RLS comms files
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  dplyr::select(SiteName, Species, common_name, abundance, Method) %>%
  rename(TaxaAb = abundance) %>%
  rename(CommonName = common_name) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName), CommonName = as.factor(CommonName), Species = as.factor(Species)) %>%
  filter(SiteName != "Second Beach South")

# Loading the env predictors file
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom), Cluster = as.factor(Cluster)) %>%
  mutate(Area_scaled = scale(Area_m2)) %>%
  dplyr::select(SiteName, Kelpdom, Depth_datum_m, exp_36, Tempave, Phardbottom, Psoftbottom, Punderstory, DensityM, BiomassM, Area_m2, HeightM) %>%
  mutate(across(Depth_datum_m:HeightM, scale)) # standardizing all predictors prior to GLMMs
  

####

### Grouping for: Pelagic fish: Total richness, total abundance
taxapfish_site <- taxap %>%
  group_by(SiteName) %>%
  summarise(TotalAb = sum(TaxaAb), rich = n_distinct(Species)) %>%
  merge(preds_ord, by="SiteName") 
  
## Generating the environmental model for total pelagic fish abundance
pfish_envmod <- glmmTMB(TotalAb ~ Kelpdom + Depth_datum_m + Tempave + exp_36 + Phardbottom + (1|SiteName),
        family="poisson", data=taxapfish_site)
check_collinearity(pfish_envmod) # no high collinearity of terms
hist(residuals(pfish_envmod)) # looking good
testDispersion(pfish_envmod) 
simulationOutput <- simulateResiduals(fittedModel = pfish_envmod, plot = F)
plot(simulationOutput) # looks relatively normal and homogeneous 
summary(pfish_envmod)
Anova(pfish_envmod)

sig <-  # sig vars
# Making the coefficient plot
pfish_envdw <- dwplot(pfish_envmod, ci=0.95, effects="fixed") +
  geom_point(aes(x=estimate, y=term), color="black", size=2) +
  geom_segment(aes(x=conf.low, y=term, xend=conf.high, yend=term), color="black", size=1) +
  scale_x_continuous(limits=c(-3,3)) +
  theme_classic() +
  geom_vline(xintercept=0, linetype=2) +
  theme(legend.position="none",
        axis.text=element_text(size=11, color="black"),
        axis.title.y=element_text(size=14, color="black", face="bold")) +
  scale_y_discrete(labels=rev(c("Kelp canopy spp", "Depth", "Temp", "Wave exp", "Hard benthos"))) +
  annotate("text", x = -3, y = 5.5, label = "a", size=5, fontface=2) +
  ylab("Environmental")
pfish_envdw

pfish_envdw


## Generating the kelp structural model for total pelagic fish abundance
pfish_strmod <- glmmTMB(TotalAb ~ Kelpdom + DensityM + Area_m2 + BiomassM + HeightM + (1|SiteName),
                        family="poisson", data=taxapfish_site)
check_collinearity(pfish_strmod) # no high collinearity of terms
hist(residuals(pfish_strmod)) # looking good
testDispersion(pfish_strmod) 
simulationOutput <- simulateResiduals(fittedModel = pfish_strmod, plot = F)
plot(simulationOutput) # looks relatively normal and homogeneous 
summary(pfish_strmod)
Anova(pfish_strmod)

# Making the coefficient plot
pfish_strdw <- dwplot(pfish_strmod, ci=0.95, effects="fixed") +
  geom_point(aes(x=estimate, y=term), col="black", size=2) +
  geom_segment(aes(x=conf.low, y=term, xend=conf.high, yend=term), col="black", size=1) +
  scale_x_continuous(limits=c(-3,3)) +
  theme_classic() +
  geom_vline(xintercept=0, linetype=2) +
  theme(legend.position="none",
        axis.text=element_text(size=11, color="black"),
        axis.title.y=element_text(size=14, color="black", face="bold")) +
  scale_y_discrete(labels=rev(c("Kelp canopy spp", "Density", "Area", "Biomass", "Height"))) +
  annotate("text", x = -3, y = 5.5, label = "d", size=5, fontface=2) +
  ylab("Structural")
pfish_strdw




### Grouping for: Demersal fish: Total richness, total abundance
taxabfish_site <- taxab_fish %>%
  group_by(SiteName) %>%
  summarise(TotalAb = sum(TaxaAb), rich = n_distinct(Species)) %>%
  merge(preds_ord, by="SiteName")

## Generating the environmental model for total pelagic fish abundance
bfish_envmod <- glmmTMB(TotalAb ~ Kelpdom + Depth_datum_m + Tempave + exp_36 + Phardbottom + (1|SiteName),
                        family="poisson", data=taxabfish_site)
check_collinearity(bfish_envmod) # no high collinearity of terms
hist(residuals(bfish_envmod)) # looking good
testDispersion(bfish_envmod) 
simulationOutput <- simulateResiduals(fittedModel = bfish_envmod, plot = F)
plot(simulationOutput) # looks relatively normal and homogeneous 
summary(bfish_envmod)
Anova(bfish_envmod)

# Making the coefficient plot
bfish_envdw <- dwplot(bfish_envmod, ci=0.95, effects="fixed") +
  geom_point(aes(x=estimate, y=term), col="black", size=2) +
  geom_segment(aes(x=conf.low, y=term, xend=conf.high, yend=term), col="black", size=1) +
  scale_x_continuous(limits=c(-3,3)) +
  theme_classic() +
  geom_vline(xintercept=0, linetype=2) +
  theme(legend.position="none",
        axis.text.y=element_text(size=11, color="white"),
        axis.text.x=element_text(size=11, color="black"),
        axis.title.y=element_blank()) +
  scale_y_discrete(labels=rev(c("Kelp canopy spp", "Depth", "Temp", "Wave exp", "Hard benthos"))) +
  annotate("text", x = -3, y = 5.5, label = "b", size=5, fontface=2)
bfish_envdw


## Generating the kelp structural model for total pelagic fish abundance
bfish_strmod <- glmmTMB(TotalAb ~ Kelpdom + DensityM + Area_m2 + BiomassM + HeightM + (1|SiteName),
                        family="poisson", data=taxabfish_site)
check_collinearity(bfish_strmod) # no high collinearity of terms
hist(residuals(bfish_strmod)) # ehh
testDispersion(bfish_strmod) 
simulationOutput <- simulateResiduals(fittedModel = bfish_strmod, plot = F)
plot(simulationOutput) # looks relatively normal and homogeneous 
summary(bfish_strmod)
Anova(bfish_strmod)

# Making the coefficient plot
bfish_strdw <- dwplot(bfish_strmod, ci=0.95, effects="fixed") +
  geom_point(aes(x=estimate, y=term), col="black", size=2) +
  geom_segment(aes(x=conf.low, y=term, xend=conf.high, yend=term), col="black", size=1) +
  scale_x_continuous(limits=c(-3,3)) +
  theme_classic() +
  geom_vline(xintercept=0, linetype=2) +
  theme(legend.position="none",
        axis.text.y=element_text(size=11, color="white"),
        axis.text.x=element_text(size=11, color="black"),
        axis.title.y=element_blank()) +
  scale_y_discrete(labels=rev(c("Kelp canopy spp", "Density", "Area", "Biomass", "Height"))) +
  annotate("text", x = -3, y = 5.5, label = "e", size=5, fontface=2)
bfish_strdw


test <- lm(HeightM ~ Tempave + exp_36 + Depth_datum_m, data=taxabinv_site)
summary(test)
plot(taxabinv_site$HeightM ~ taxabinv_site$exp_36)

### Grouping for: Benthic inverts: Total richness, total abundance
taxabinv_site <- taxab_inv %>%
  group_by(SiteName) %>%
  summarise(TotalAb = sum(TaxaAb), rich = n_distinct(Species)) %>%
  merge(preds_ord, by="SiteName")

## Generating the environmental model for total invert abundance
binv_envmod <- glmmTMB(TotalAb ~ factor(Kelpdom) + Depth_datum_m + Tempave + exp_36 + Phardbottom,
                        family="nbinom1", data=taxabinv_site) # negative binomial (quasi-poisson) because of overdispersion in the residuals, random effect removed because of extremely low variation (singularity)
check_collinearity(binv_envmod) # no high collinearity of terms
hist(residuals(binv_envmod)) # looking good
testDispersion(binv_envmod) 
simulationOutput <- simulateResiduals(fittedModel = binv_envmod, plot = F)
plot(simulationOutput) # looks relatively normal and homogeneous now
summary(binv_envmod)
Anova(binv_envmod)

# making the coefficient plot
binv_envdw <- dwplot(binv_envmod, ci=0.95, effects="fixed") +
  geom_point(aes(x=estimate, y=term), col="black", size=2) +
  geom_segment(aes(x=conf.low, y=term, xend=conf.high, yend=term), col="black", size=1) +
  scale_x_continuous(limits=c(-3,3)) +
  theme_classic() +
  geom_vline(xintercept=0, linetype=2) +
  theme(legend.position="none",
        axis.text.y=element_text(size=11, color="white"),
        axis.text.x=element_text(size=11, color="black"),
        axis.title.y=element_blank()) +
  scale_y_discrete(labels=rev(c("Kelp canopy spp", "Depth", "Temp", "Wave exp", "Hard benthos"))) +
  annotate("text", x = -3, y = 5.5, label = "c", size=5, fontface=2)
binv_envdw


## Generating the kelp structural model for total invert abundance
binv_strmod <- glmmTMB(TotalAb ~ Kelpdom + DensityM + Area_m2 + BiomassM + HeightM + (1|SiteName),
                        family="poisson", data=taxabinv_site) # negative binomial (quasi-poisson) because of overdispersion in the residual
check_collinearity(binv_strmod) # no high collinearity of terms
hist(residuals(binv_strmod)) # looking good
testDispersion(binv_strmod) 
simulationOutput <- simulateResiduals(fittedModel = binv_strmod, plot = F)
plot(simulationOutput) # looks relatively normal and homogeneous 
summary(binv_strmod)
Anova(binv_strmod)

# Making the coefficient plot
binv_strdw <- dwplot(binv_strmod, ci=0.95, effects="fixed") +
  geom_point(aes(x=estimate, y=term), col="black", size=2) +
  geom_segment(aes(x=conf.low, y=term, xend=conf.high, yend=term), col="black", size=1) +
  scale_x_continuous(limits=c(-3,3)) +
  theme_classic() +
  geom_vline(xintercept=0, linetype=2) +
  theme(legend.position="none",
        axis.text.y=element_text(size=11, color="white"),
        axis.text.x=element_text(size=11, color="black"),
        axis.title.y=element_blank()) +
  scale_y_discrete(labels=rev(c("Kelp canopy spp", "Density", "Area", "Biomass", "Height"))) +
  annotate("text", x = -3, y = 5.5, label = "f", size=5, fontface=2)
binv_strdw


### Putting them all together

tiff(file="./MSc_plots/PaperFigs/GLMM_coefs.tiff", height = 7, width = 13, units = "in", res=400)

arr <- ggarrange(pfish_envdw + rremove("x.text"), bfish_envdw + rremove("x.text"), binv_envdw + rremove("x.text"),
          pfish_strdw, bfish_strdw, binv_strdw,
          ncol=3, nrow=2, align="hv")
annotate_figure(arr, bottom = text_grob("Estimates", size=14, face="bold"))

dev.off()



#### SIMPER plots ----

### SIMPER PLOT (to go with CAP plots)


# Joining all spp identified as influential by SIMPER
simp_all <- rbind(pfish_simp_filt, bfish_simp_filt, binv_simp_filt)



# Arranging by descending order of average dissimilarity 
# Making the values neg. for spp that are biased to Macro in their dissimilarity (plotting purposes)
simp_arr <- simp_all %>%
  arrange(desc(average)) %>%
  mutate(average_plot = ifelse(ava > avb, average*(-1), average)) %>% 
  ungroup() %>%
  as.data.frame() %>%
  mutate(species = factor(species, levels=c("Embiotoca lateralis", "Rhinogobiops nicholsii", "Mesocentrotus franciscanus", "Patiria miniata", "Pomaulax gibberosus", "Strongylocentrotus purpuratus"), order=TRUE))



# Recoding the spp to shorthand genus spp
# Adding new column for faunal groups
# Adding new column for kelp groups (i.e., which group a spp is more relatively abundant in)
simp_arr <- simp_arr %>%
  mutate(species_common = fct_recode(species,
                             "Striped surfperch" = "Embiotoca lateralis",
                             "Blackeye goby" = "Rhinogobiops nicholsii",
                             "Red sea urchin" = "Mesocentrotus franciscanus",
                             "Bat star" = "Patiria miniata",
                             "Red turban snail" = "Pomaulax gibberosus",
                             "Purple sea urchin" = "Strongylocentrotus purpuratus")) %>%
  mutate(faunalgroup = case_when(species_common == "Striped surfperch" ~ "Pelagic fish",
                                 species_common == "Blackeye goby" ~ "Demersal fish",
                                 species_common == "Red sea urchin" | species_common == "Bat star" | 
                                 species_common == "Red turban snail" | species_common == "Purple sea urchin" ~ "Benthic invert")) %>%
  mutate(faunalgroup = factor(faunalgroup, levels=c("Pelagic fish", "Demersal fish", "Benthic invert"), ordered=TRUE)) %>%
  mutate(species_common = factor(species_common, levels=c("Striped surfperch", "Purple sea urchin", "Red sea urchin", "Bat star", "Red turban snail", "Blackeye goby"), ordered=TRUE)) %>%
  mutate(kelpgroup = ifelse(ava > avb, "Giant kelp", "Bull kelp"))




# Loading the kelp icon .png files
macroicon <- readPNG('./Phylopic/EladTrace/macrokelp.png')
bullicon <- readPNG('./Phylopic/EladTrace/bullkelp.png')


tiff(file="./MSc_plots/PaperFigs/simper.tiff", height=4, width=6, units="in", res=500)

# Making the plot of average SIMPER dissimilarities
# Only the spp that contribute > 2x their expected dissimilarity
simpplot <- ggplot(data=simp_arr, aes(y=species_common, x=average_plot, color=kelpgroup)) +
  # geom_rect(aes(xmin = 0.008,
  #               xmax = 0.25,
  #               ymin = 2.9,
  #               ymax = 6.4), fill = '#3b7c70', alpha = 0.03, color=NA) +
  # geom_rect(aes(xmin = -0.24,
  #               xmax = -0.008,
  #               ymin = 0.5,
  #               ymax = 2.8), fill = '#ce9642', alpha = 0.03, color=NA) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin=average_plot-sd, xmax=average_plot+sd), height=0.1, linewidth=1) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(limits=c(-0.35, 0.35)) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color="black", size=9, angle=25),
        axis.text.x = element_text(color="black", size=9),
        axis.title.x = element_text(color="black", size=9.5),
        legend.position="none",
        legend.title = element_blank()) +
  scale_color_manual(values=met.brewer("Kandinsky", n=2)) +
  geom_vline(xintercept=0, linetype="dashed") +
  xlab("Average dissimilarity") +
  annotate("text", x = -0.35, y = 6.5, label = "d", size=6, fontface=2)
simpplot

dev.off()


# Adding in the kelp icons
simpplot_draw <- simpplot + 
  annotation_raster(macroicon, ymin = 1.5, ymax= 2.3, xmin = -0.26, xmax = -0.18) +
  annotation_raster(bullicon, ymin = 4.5, ymax= 5.3, xmin = 0.21, xmax = 0.29)
simpplot_draw



## Arranging the CAP plots (loaded from the faunal group code chunks above)

# Arranging the smaller plots together (ABUNDANCE)
capplots <- ggarrange(pfish1_ordplot + rremove("legend"), bfish1_ordplot + rremove("legend"), 
                      binv1_ordplot + rremove("legend"), ncol=1, align="hv")
capplots

# Setting up the layout
lay <- rbind(c(1,1,1,1,2,2,2,2,2,2,2),
             c(1,1,1,1,2,2,2,2,2,2,2))
# Rectangles to help visualize
gs <- lapply(1:2, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)


# Saving as an image object
tiff(file="./MSc_plots/PaperFigs/cap_ords.tiff", height=6.5, width=8.5, units="in", res=400)

# End ggplot arrangement
grid <- grid.arrange(capplots, simpplot_draw,
                     layout_matrix = lay)
dev.off()


#


#



#




#### Testing correlations ----

testdata <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))

plot(x=testdata$Psoftbottom, y=testdata$Area_m2)

testmod <- glmmTMB(Area_m2 ~ Psoftbottom, data=testdata)
testmod <- glmmTMB(MacroM ~ Psoftbottom, data=testdata)
testmod <- glmmTMB(BiomassM ~ Psoftbottom, data=testdata)
testmod <- glmmTMB(HeightM ~ Psoftbottom, data=testdata)
summary(testmod)




