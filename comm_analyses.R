## Barkley Sound community analyses
## Author: Claire Attridge
## Origin date: Feb 2023

# Loading base packages
library(tidyverse)
library(ggpubr)
library(grid)
library(gridExtra)
library(MetBrewer)
library(wesanderson)
library(vegan)
library(cowplot)
library(funrar)
library(ecodist)
library(BiodiversityR)
library(ggforce)
library(ggrepel)
library(magick)

#

#### Loading all data sheets (kelp forest structure & environmental) ----

# Wave exp
rei <- read_csv("./MSc_data/Data_new/REI_2022.csv") %>%
  dplyr::select(-c(x,y)) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName), exp_scaled = scale(exp_36)) %>% # Scaling REI
  mutate(exp_scaled = as.numeric(exp_scaled))

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


# Removing second beach south as an outlier site 
preds_sites <- preds %>%
  filter(SiteName != "Second Beach South")
#### Removing the 'kelp free sites'
preds_sites <- preds_sites %>%
  filter(SiteName != "Wizard Islet North" & SiteName != "Less Dangerous Bay")

# Adding Cluster column to dataframe for later use (see section below for determination of Cluster groups C1-C4) & dominant kelp column
preds_ord <- preds %>%
  # mutate(Cluster = case_when(SiteName == "Between Scotts and Bradys" | SiteName == "Danvers Danger Rock" | SiteName == "Dixon Island Back (Bay)" | SiteName == "Dodger Channel 1" | SiteName == "Flemming 112" | SiteName == "Flemming 114" | SiteName == "North Helby Rock" | SiteName == "Tzartus 116" ~ "C2",
  #                            SiteName == "Ed King East Inside" | SiteName == "Ross Islet 2" | SiteName == "Ross Islet Slug Island" | SiteName == "Taylor Rock" | SiteName == "Turf Island 2" | SiteName == "Wizard Islet South" ~ "C1",
  #                            SiteName == "Bordelais Island" | SiteName == "Second Beach" | SiteName == "Dodger Channel 2" | SiteName == "Nanat Bay" ~ "C3",
  #                            SiteName == "Cable Beach (Blow Hole)" | SiteName == "Less Dangerous Bay" | SiteName == "Swiss Boy" | SiteName == "Wizard Islet North" ~ "C4",
  #                            TRUE ~ "Aux")) %>%
  # mutate(Cluster = as.factor(Cluster), Composition = as.factor(Composition)) %>%
  mutate(Kelpdom = case_when(Composition == "Macro" | Composition == "Mixed" ~ "Giant", # Macro is dominant at mixed sites (too few mixed sites to statistically compare)
                             Composition == "Nereo" ~ "Bull",
                             Composition == "None" ~ "Giant")) %>% # Re. Starko (2022) - used to have M. pyrifera present around the area of Less Dangerous bay
  mutate(Kelpdom = as.factor(Kelpdom)) %>%
  droplevels()


# Removing second beach south as an outlier site 
preds_ord <- preds_ord %>%
  filter(SiteName != "Second Beach South") %>%
  droplevels()
#### Removing the 'kelp free sites'
preds_ord <- preds_ord %>%
  filter(SiteName != "Wizard Islet North" & SiteName != "Less Dangerous Bay") %>%
  droplevels()


# saving a .csv file of the ecological and environmental predictors of interest for our community abundance data
write.csv(preds_ord, "./MSc_data/Data_new/AllPredictors_2022.csv", row.names=F)


#### Cluster analysis (kelp forest structure) *IGNORE* ----

# Plotting out individual structural variable relationships
# Biomass - Density
mod <- lm(DensityM ~ BiomassM, data=preds)
modpred <- predict(mod)  
plot(DensityM ~ BiomassM, data=preds)
lines(preds$BiomassM, modpred)
# Biomass - Height
mod <- lm(HeightM ~ BiomassM, data=preds)
modpred <- predict(mod)  
plot(HeightM ~ BiomassM, data=preds)
lines(preds$BiomassM, modpred)
# Density - Height
mod <- lm(HeightM ~ DensityM, data=preds)
modpred <- predict(mod)  
plot(HeightM ~ DensityM, data=preds)
lines(preds$DensityM, modpred)

## Strong correlation between Biomass & Density
# Vars likely to be collinear and may overemphasize the same underlying contribution
# Biomass removed for the following cluster analysis


# ## PCA (for density & biomass)
# pca <- prcomp(na.omit(preds[,c(3:4)]), center = TRUE, scale = TRUE)
# fviz_eig(pca) # Scree plot: PC1 explains vast majority of variation (89.9%)
# fviz_pca_biplot(pca) # Biplot: Pos. PC1 values = higher density & biomass
# 
# pc1 <- as.vector(pca$x[,1]) # Extracting the pc1 values for further use
# 
# preds$PC1 <- pc1 # Adding pc1 values to the working dataframe
  

## K MEANS CLUSTERING
preds <- preds %>% # Scaling all kelp forest numeric predictors (Euc dists are sens to scale)
  mutate_at(c(2,3,4,5,6,7,9,10), funs(c(scale(.))))

# optimizing for 'k' cluster values
rng <- 2:15 # K from 2 to 15
tries <- 100 #Run the K Means algorithm 100 times
avg.totw.ss <-integer(length(rng)) #Set up an empty vector to hold all of points

for(v in rng){ # For each value of the range variable
  v.totw.ss <- integer(tries) # Set up an empty vector to hold the 100 tries
  for(i in 1:tries){
    k.temp <- kmeans(na.omit(preds[,c(3,5,9,10)]), centers=v) # Run kmeans (3 = height, 5 = density, 9 = area, 10 = proportion M. pyrifera)
    v.totw.ss[i] <-k.temp$tot.withinss # Store the total withinss
  }
  avg.totw.ss[v-1] <-mean(v.totw.ss) #Average the 100 total withinss
}

# plot to show ave within SS as K increases (when is within SS reduced as we ^ clusters)
plot(rng,avg.totw.ss,type="b", main="Total Within SS by Various K",
     ylab="Average Total Within Sum of Squares",
     xlab="Value of K")
# 4 or more clusters will substantially minimize (>50% decrease) error within clusters


# only using the kelp forest attributes for clustering sites
k <- kmeans(na.omit(preds[,c(3,5,9,10)]), centers=4, nstart=20) # specifying 4 clusters
k$centers # looking at the centroid (group) mean values
table(k$cluster) # looking at number of sites per group


cols <- c("#997700", "#225555", "#4477aa", "#e4632d")

tiff(file="./MSc_plots/Clusters_kmeans4.tiff", height = 4.5, width = 8, units = "in", res=400)

# plotting the clustered groups 
clusts <- fviz_cluster(k, data=na.omit(preds[,c(3,5,9,10)]), pointsize=2, repel=TRUE) +
  theme_classic() +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme(
    legend.position="top",
    axis.text = element_text(color="black", size="12"),
    axis.title = element_text(color="black", size="13", vjust=1),
    plot.margin = unit(c(0,1,0.5,0.85), "cm"), # Selecting margins
    # panel.border = element_rect(colour="black", fill=NA, linewidth=1), # Border option
    axis.line = element_line(color="black", linewidth=0.5)) + 
  ggtitle("") +
  annotate("text", label = "C1", size=6, fontface=2, x = -2.2, y = 0.9) +
  annotate("text", label = "C2", size=6, fontface=2, x = -0.2, y = 1.4) +
  annotate("text", label = "C4", size=6, fontface=2, x = 1.5, y = 1.6) +
  annotate("text", label = "C3", size=6, fontface=2, x = -1.2, y = -2.1) 
  # annotate("text", label = "C5", size=6, fontface=2, x = 1.8, y = -2.0) 
# Dim1 - Dens & prop Macro (reversed)
# Dim2 - Height & area (reversed)
clusts

dev.off()

## Clusters for- density, area, height, macro (proportion)

# # C1:
# Ed King East Inside, Ross Islet 2, Ross Islet Slug Island, Taylor Rock, Turf Island 2, Wizard Islet South
# # C2: 
# Between Scotts and Bradys, Danvers Danger Rock, Dixon Island Back (Bay), Dodger Channel 1, Flemming 112, Flemming 114, North Helby Rock, Tzartus 116
# # C3: 
# Cable Beach (Blow Hole), Less Dangerous Bay, Swiss Boy, Wizard Islet North
# # C4:
# Bordelais Island, Second Beach, Dodger Channel 2, Nanat Bay


## IMAGE PLOTS FOR CLUSTERS

library(imager)
library(magick)
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html

#C1
img1 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251965.JPG')
img1bd <- img1 %>%
  magick::image_border(color="#225555", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C2
img2 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251975.JPG')
img2bd <- img2 %>%
  magick::image_border(color="#4477aa", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C3
img3 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251989.JPG')
img3bd <- img3 %>%
  magick::image_border(color="#e4632d", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C5
img4 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251994.JPG')
img4bd <- img4 %>%
  magick::image_border(color="#997700", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

imgs <- c(img1bd,img2bd,img3bd,img4bd) # Combining all imgs to magick vector

imgstk <- magick::image_append(image_scale(imgs, "480"), stack=TRUE) # Stacking the images vertically


# # optional formats to shift to

# imgras <- as.raster(imgstk) 

# imggrb <- rasterGrob(imgras)


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

# Filtering out the site redo surveys at Ross Islet 2, Less Dangerous Bay, & Second Beach S in Sept 2022, and the acoustic baseline acoustic baseline site 'Sand Town'
comms <- comms %>%
  filter(Date <= as.POSIXct("2022-09-08") & SiteName != "Sand Town") %>%
  droplevels()

# Removing second beach south as an outlier site 
comms <- comms %>%
  filter(SiteName != "Second Beach South") %>%
  droplevels()
#### Removing the 'kelp free sites'
comms <- comms %>%
  filter(SiteName != "Wizard Islet North" & SiteName != "Less Dangerous Bay") %>%
  droplevels()

# Grouping Henricia spp. and Henricia leviuscula to a combined Henricia spp. as we can't be confident in the separation of these groups
comms <- comms %>%
  mutate(Species = fct_collapse(Species, 'Henricia spp.'= c("Henricia spp.", "Henricia leviuscula")),
         common_name = fct_collapse(common_name, 'Blood star' = c("Pacific blood star", "Unidentified blood star"))) %>%
  droplevels()


## Looking into some general questions first:

# Which site has the most species?
richness <- comms %>%
  group_by(SiteName) %>%
  summarize(spp_rich = n_distinct(Species))%>%
  arrange(desc(spp_rich))

# Which species were most abundant overall?
abun <- comms %>%
  group_by(Species) %>%
  summarize(totalab = sum(abundance)) %>%
  arrange(desc(totalab))

# Which species had the most site observations?
obvs <- comms %>%
  group_by(Species) %>%
  summarize(observs = length(Species)) %>% # The number of rows (observations) per spp common name
  arrange(desc(observs)) %>% # Arranging in descending order
  ungroup() %>%
  as.data.frame()


# Plot for species observations - Supp Fig
tiff(file="./MSc_plots/SuppFigs/SpeciesAbCutoff.tiff", height = 9, width = 7, units="in", res=400)

obvsplot <- ggplot(data=obvs, aes(x=observs, y=(fct_reorder(Species, desc(-observs))), fill=observs)) + 
  geom_bar(stat="identity", position="dodge", width=0.4) +
  scale_fill_gradientn(colors=met.brewer("VanGogh3")) +
  theme_classic() +
  theme(legend.position="none",
    axis.text.x = element_text(color="black", size="10"),
    axis.text.y = element_text(color="black", size="7"),
    axis.title.y = element_text(color="black", size="11", vjust=0.5),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"), 
    ) +
  xlab("Observations") + ylab("") +
  geom_hline(yintercept = 55.5, linetype=2) # Adding in the cutoff line for =< 5 obvs
obvsplot

dev.off()


## With this ^^ in mind: Removing rarely sighted species (< 5 site obvs) from the analysis frame (n = 38 removed)
sppkeep <- obvs %>%
  filter(observs >= 5)


# Remaining species (n = 50 retained)
commsclean <- comms %>%
  filter(Species %in% sppkeep$Species)


# Finally: What are the total abundance counts for all species at all sites?
comms_all <- commsclean %>%
  as.data.frame() %>%
  ungroup() %>%
  group_by(SiteName, Species, common_name, Method) %>%
  summarise(abundance = sum(abundance))


# saving a .csv file of the clean total taxa abundances by site
write.csv(comms_all, "./MSc_data/Data_new/RLS_2022.csv", row.names=F)


#### Community data: Site/Species groupings ----

# Loading the RLS comms file
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  dplyr::select(SiteName, Species, common_name, abundance, Method) %>%
  rename(TaxaAb = abundance) %>%
  rename(CommonName = common_name) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName), CommonName = as.factor(CommonName), Species = as.factor(Species)) %>%
  filter(SiteName != "Second Beach South")


### For all species in community together
taxa_all <- taxa %>%
  dplyr::select(-Method) %>% # Removing RLS method col when using all species together
  group_by(SiteName, Species, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb)) %>%
  ungroup() %>%
  as.data.frame()
# Cleaning up the all spp frame
taxa_all <- taxa_all %>%
  mutate(CommonName = fct_recode(CommonName, "Green sea urchin" = "Northern sea urchin"),
         CommonName = fct_recode(CommonName, "Chiton" = "Unidentified chiton"),
         CommonName = fct_recode(CommonName, "Hermit crab" = "Unidentified hermit crab")) %>%
  # filter(CommonName != "Unidentified rockfish") %>% # Removing unidentified rockfish as a spp.
  droplevels()


### For just Pelagic fish spp (RLS method 1)
taxap <- taxa %>%
  filter(Method == 1) %>%
  # filter(CommonName != "Unidentified rockfish") %>% # Removing unidentified rockfish as a spp.
  filter(CommonName != "Blackeye goby") %>% # Removing a benthic associated spp
  filter(CommonName != "Longfin sculpin") %>% # Removing a benthic associated spp
  filter(CommonName != "Red Irish lord") %>% # Removing a benthic associated spp
  filter(CommonName != "Painted greenling") %>% # Removing a benthic associated spp
  filter(CommonName != "Kelp greenling") %>% # Removing a benthic associated spp
  filter(CommonName != "Whitespotted greenling") %>% # Removing a benthic associated spp
  filter(CommonName != "Lingcod") %>% # Removing a benthic associated spp
  droplevels() %>%
  dplyr::select(-Method) 


### For just demersal fish & benthic invert spp (RLS Method 2)
taxab <- taxa %>%
  filter(Method == 2) %>%
  mutate(CommonName = fct_recode(CommonName, "Green sea urchin" = "Northern sea urchin"),
         CommonName = fct_recode(CommonName, "Chiton" = "Unidentified chiton"),
         CommonName = fct_recode(CommonName, "Hermit crab" = "Unidentified hermit crab")) %>%
  droplevels() %>%
  dplyr::select(-Method)


### For just benthic invert spp (RLS Method 2 - demersal fish)
M2benthicfish <- c("Blackeye goby", "Crescent gunnel", "Kelp greenling", "Longfin sculpin", "Penpoint gunnel", "Red Irish lord", "Scalyhead sculpin", "Smoothhead sculpin", "Whitespotted greenling", "Lingcod", "Painted greenling")

taxab_inv <- taxab %>%
  filter(!CommonName %in% M2benthicfish) %>%
  droplevels()
taxab_inv <- taxab_inv %>%
  group_by(SiteName, Species, CommonName) %>%
    summarise(TaxaAb = sum(TaxaAb))


### For just demersal fish spp
taxap_benthic <- taxa %>% # Isolating the demersal fish spp obvs from M1 surveys
  filter(Method == 1) %>%
  filter(CommonName %in% M2benthicfish) %>%
  dplyr::select(-Method) %>%
  ungroup()

taxab_fish <- taxab %>%
  rbind(taxap_benthic) %>% # Adding on the M1 obvs of demersal fish spp
  filter(CommonName %in% M2benthicfish) %>%
  droplevels() %>%
  group_by(SiteName, Species, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb)) # Summing M1 & M2 pelagic fish data to sites


## For all fish spp (pelagic + demersal fish)
taxa_fish <- taxap %>%
  rbind(taxab_fish) %>%
  group_by(SiteName, Species, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb))


#### Community data: Relative abundance tables & transformations ----

## ALL SPP
allspp_wide <- taxa_all %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>% # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(allspp_wide) <- allspp_wide$SiteName # Setting col as rownames
allspp_wide <- allspp_wide[,-1] # Removing the col used above
allspp_wide[is.na(allspp_wide)] <- 0 # (no NAs)
allspp_mat <- as.matrix(allspp_wide) # As matrix

allspp_rel <- make_relative(allspp_mat) # Relative abundance matrix of community data
allspp_hel <- decostand(allspp_mat, method = "hellinger") # Hellinger transformation of data

## PELAGIC FISH SPP
pfish_wide <- taxap %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>% # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(pfish_wide) <- pfish_wide$SiteName # Setting col as rownames
pfish_wide <- pfish_wide[,-1] # Removing the col used above
pfish_wide[is.na(pfish_wide)] <- 0 # (no NAs)
pfish_mat <- as.matrix(pfish_wide) # As matrix

pfish_rel <- make_relative(pfish_mat) # Relative abundance matrix of community data
pfish_hel <- decostand(pfish_mat, method = "hellinger") # Hellinger transformation of data

## BENTHIC FISH SPP
bfish_wide <- taxab_fish %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>%  # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(bfish_wide) <- bfish_wide$SiteName # Setting col as rownames
bfish_wide <- bfish_wide[,-1] # Removing the col used above
bfish_wide[is.na(bfish_wide)] <- 0 # (no NAs)
bfish_mat <- as.matrix(bfish_wide) # As matrix

bfish_rel <- make_relative(bfish_mat) # Relative abundance matrix of community data
bfish_hel <- decostand(bfish_mat, method = "hellinger") # Hellinger transformation of square root

## ALL FISH SPP
fish_wide <- taxa_fish %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>% # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(fish_wide) <- fish_wide$SiteName # Setting col as rownames
fish_wide <- fish_wide[,-1] # Removing the col used above
fish_wide[is.na(fish_wide)] <- 0 # (no NAs)
fish_mat <- as.matrix(fish_wide) # As matrix

fish_rel <- make_relative(fish_mat) # Relative abundance matrix of community data
fish_hel <- decostand(fish_mat, method = "hellinger") # Hellinger transformation of square root

## BENTHIC INVERT SPP
binv_wide <- taxab_inv %>% # Spreading the data to wide format
  dplyr::select(-CommonName) %>% # Removing common names for now
  spread(key = Species, value = TaxaAb) %>%
  as.data.frame()
rownames(binv_wide) <- binv_wide$SiteName # Setting col as rownames
binv_wide <- binv_wide[,-1] # Removing the col used above
binv_wide[is.na(binv_wide)] <- 0 # (no NAs)
binv_mat <- as.matrix(binv_wide) # As matrix

binv_rel <- make_relative(binv_mat) # Relative abundance matrix of community data
binv_hel <- decostand(binv_mat, method = "hellinger") # Hellinger transformation of data


#### Community data: Abundance & rel abundance plots ----

## Grouping the table for pelagic fish
taxap_fish_frame <- taxap %>%
  ungroup() %>%
  as.data.frame()
# group_by(Species, CommonName) %>%
# summarise(TaxaAb_av = mean(TaxaAb), SD = sd(TaxaAb), Sites = n())

colourCount = length(unique(taxap_fish_frame$Species))
getPalette = colorRampPalette(met.brewer(name="Thomas", n=8))

# Making the regional abundance plot for pelagic fish
abPfish <- ggplot(taxap_fish_frame, 
       aes(y=as.factor(reorder(Species, TaxaAb, FUN=median)), x=TaxaAb, fill=as.factor(Species))) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.2),
              aes(fill=as.factor(Species), alpha=0.9), color="grey55", show.legend=FALSE, shape=21) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.shape=NA, fatten=1) + # hiding outlier points since already showing with the jittered raw data
  scale_fill_manual(values=getPalette(colourCount)) +
  theme_classic() +
  theme(strip.text.x = element_text(size=9, color="black", face="bold"),
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=9, color="black"),
        axis.title.y = element_blank()
  )


## Grouping for number of site occurrences for each species
taxap_fish_occ <- taxap_fish_frame %>%
  dplyr::select(SiteName, Species, TaxaAb) %>%
  group_by(Species) %>%
  summarise(Mean=mean(TaxaAb), SiteNum=n()) %>%
  ungroup() %>%
  mutate(Prop = round((SiteNum/22)*100, digits=1))

# Making the prop occurrence plot for pelagic fish
ocPfish <- ggplot(taxap_fish_occ, aes(x=Species, y=Prop)) +
  geom_segment( aes(x=reorder(Species, Prop), xend=Species, y=0, yend=Prop), color="#453781FF", linewidth=1) +
  geom_point(color="#453781FF", size=3) +
  theme_classic() +
  scale_size_continuous("Mean Ab", breaks = c(100,200,400,600), ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black", angle=60, hjust=1, size=9),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_blank()
  ) +
  annotate("text", label="A", x=1, y=100, fontface="bold", size=6)



## Grouping the table for demersal fish
taxab_fish_frame <- taxab_fish %>%
  ungroup() %>%
  as.data.frame()
  # group_by(Species, CommonName) %>%
  # summarise(TaxaAb_av = mean(TaxaAb), SD = sd(TaxaAb), Sites = n())

colourCount = length(unique(taxab_fish_frame$Species))
getPalette = colorRampPalette(met.brewer(name="Nattier", n=8))

# Making the regional abundance plot for demersal fish
abBfish <- ggplot(taxab_fish_frame, 
       aes(y=as.factor(reorder(Species, TaxaAb, FUN=median)), x=TaxaAb, fill=as.factor(Species))) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.2),
              aes(fill=as.factor(Species), alpha=0.9), color="grey55", show.legend=FALSE, shape=21) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.shape=NA, fatten=1) + # hiding outlier points since already showing with the jittered raw data
  scale_fill_manual(values=getPalette(colourCount)) +
  theme_classic() +
  theme(strip.text.x = element_text(size=9, color="black", face="bold"),
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=9, color="black"),
        axis.title.y = element_blank()
  )


## Grouping for number of site occurrences for each species
taxab_fish_occ <- taxab_fish_frame %>%
  dplyr::select(SiteName, Species, TaxaAb) %>%
  group_by(Species) %>%
  summarise(Mean=mean(TaxaAb), SiteNum=n()) %>%
  ungroup() %>%
  mutate(Prop = round((SiteNum/22)*100, digits=1))

# Making the prop occurrence plot for demersal fish
ocBfish <- ggplot(taxab_fish_occ, aes(x=Species, y=Prop)) +
  geom_segment( aes(x=reorder(Species, Prop), xend=Species, y=0, yend=Prop), color="#1F9A8AFF", linewidth=1) +
  geom_point(color="#1F9A8AFF", size=3) +
  theme_classic() +
  scale_size_continuous("Mean Ab", breaks = c(100,200,400,600), ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black", angle=60, hjust=1, size=9),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_blank()
  ) +
  annotate("text", label="B", x=1, y=100, fontface="bold", size=6)



## Grouping the table for benthic inverts
taxab_inv_frame <- taxab_inv %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(split = factor(ifelse(Species == "Strongylocentrotus purpuratus" | Species == "Nucella lamellosa" | Species == "Patiria miniata" | Species == "Pomaulax gibberosus" | Species == "Mesocentrotus franciscanus", "large", "small")))

colourCount = length(unique(taxab_inv_frame$Species))
getPalette = colorRampPalette(met.brewer(name="Redon", n=12))


# Making the regional abundance plot for benthic inverts
abBinv <- ggplot(taxab_inv_frame, 
       aes(y=as.factor(reorder(Species, TaxaAb, FUN=median)), x=TaxaAb, fill=as.factor(Species))) + 
  geom_jitter(position=position_jitter(width=0.3, height=0.2),
              aes(fill=as.factor(Species), alpha=0.9), color="grey55", show.legend=FALSE, shape=21) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.shape=NA, fatten=1) + # hiding outlier points since already showing with the jittered raw data
  scale_x_continuous()+
  scale_fill_manual(values=getPalette(colourCount)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=9, color="black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(), # removing the facet boxes
        strip.text = element_blank(), # removing the facet text
  ) + 
  facet_col(~split, scales="free", space="free")


## Grouping for number of site occurrences for each species
taxab_inv_occ <- taxab_inv_frame %>%
  dplyr::select(SiteName, Species, TaxaAb) %>%
  group_by(Species) %>%
  summarise(Mean=mean(TaxaAb), SiteNum=n()) %>%
  ungroup() %>%
  mutate(Prop = round((SiteNum/22)*100, digits=1))

# Making the prop occurrence plot for demersal fish
ocBinv <- ggplot(taxab_inv_occ, aes(x=Species, y=Prop)) +
  geom_segment( aes(x=reorder(Species, Prop), xend=Species, y=0, yend=Prop), color="#31688EFF", linewidth=1) +
  geom_point(color="#31688EFF", size=3) +
  theme_classic() +
  scale_size_continuous("Mean Ab", breaks = c(100,200,400,600)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color="black", angle=60, hjust=1, size=9),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_blank()
  ) +
  annotate("text", label="C", x=1.5, y=100, fontface="bold", size=6)



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
tiff(file="./MSc_plots/SuppFigs/Ab_plots.tiff", height=8, width=8.5, units="in", res=400)

# End ggplot arrangement
grid <- grid.arrange(abfish, abBinv,
                     layout_matrix = lay, bottom=textGrob("Abundance"))
dev.off()



# Arranging the smaller plots together (OCCURRENCE)
ocfish <- ggarrange(ocPfish, ocBfish, ncol=2, align="hv")
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
                     layout_matrix = lay, left=textGrob("Occurrence (% of 22 sites)", rot=90, vjust=1, hjust=0.02))

dev.off()


#### Community data: Variable selection *IGNORE* ----

## ALL SPECIES

# Specifying the full and null models
fullord <- dbrda(binv_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray")
nullord <- dbrda(binv_hel ~ 1, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray")

# Testing whether the full model is significant
anova.cca(fullord) # Significant (p=0.007)
adjR2_ord <- RsquareAdj(fullord)$adj.r.squared
adjR2_ord # The adj R2 explained by all variables is 23.7%

# Running forward selection for all ecological/environmental variables
sel.ord <- ordiR2step(nullord, scope = formula(fullord), R2scope =
                        adjR2_ord, direction = 'forward', permutations = 999)
sel.ord$anova # Looking at the summary table


## PELAGIC FISH SPP

## BENTHIC FISH SPP





#### All species: Generating ordination models, stats ----

# Proportion of constrained variation: The amount of variation in spp abundances explained across sites by the included predictor variables 
# Adjusted R Squared: Amount of variation in the relative abundances of spp explained by the predictor variables when corrected for the number of predictor variables 
# For anova.cca() by term: https://stats.stackexchange.com/questions/405024/what-null-hypothesis-does-vegananova-ccaby-terms-test-against

## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))

# Running PERMANOVA to test for initial grouping factors
allspp_perm <- adonis2(allspp_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom significant (0.002)


### CAP (kelp spp)

### CAP (kelp spp)
allspp_1 <- capscale(allspp_hel ~ Kelpdom, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray") 
RsquareAdj(allspp_1)
## Only the environmental variables explain 11.2% of variation (constrained axes)
# Adjusted R squared of this model is: 6.1%

## SIMPER for dissimilarity
# Running SIMPER to test between groups (limiting to 10 most influential)
allspp_simp <- with(preds_ord, simper(allspp_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
allspp_simp_df <- as.data.frame(allspp_simp$Giant_Bull)
sum(allspp_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.5023465
# Partitioned among the 55 spp = 0.009133573 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
allspp_simp_filt <- allspp_simp_df %>%
  filter(average > (2*0.009133573))

allspp_simpSPP <- c("Clupea pallasii", "Cymatogaster aggregata", "Mesocentrotus franciscanus", "Patiria miniata", "Pomaulax gibberosus", "Rhinogobiops nicholsii", "Strongylocentrotus purpuratus")


### RDA (kelp forest structure)

# Generating the kelp forest structure model
allspp_2 <- capscale(allspp_hel ~ Kelpdom+DensityM+Area_m2+HeightM+BiomassM, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray") 
RsquareAdj(allspp_2)
## The structural and environmental variables explain 37.3% of variation (constrained axes)
# Adjusted R squared of this model is: 18.7%
vif.cca(allspp_2) # low vif scores

# Exploring model significance
anova.cca(allspp_2, permutations=999) # significance of model
anova.cca(allspp_2, step = 1000, by = "term", permutations=999) # significance of predictor terms
anova.cca(allspp_2, step = 1000, by = "axis", permutations=999) # significance by model axes


### RDA (environmental vars)

# Generating the environmental model
allspp_3 <- capscale(allspp_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray")
RsquareAdj(allspp_3)
## Only the environmental variables explain explain 46.4% of variation (constrained axes)
# Adjusted R squared of this model is: 21.0%
vif.cca(allspp_3) # low vif scores

# Exploring model significance
anova.cca(allspp_3) # significance of model
anova.cca(allspp_3, step = 1000, by = "term") # significance of predictor terms
anova.cca(allspp_3, step = 1000, by = "axis") # significance by model axes


### ### ###

# checking/visualizing continuous predictor correlations (if all in a single model together, i.e., kelp forest structure & environmental vars)
pairs.panels(preds_ord[,c(3:5,9)], scale=T) # kelp forest structure vars
pairs.panels(preds_ord[,c(2,11,14,16,18,19)], scale=T) # environmental vars
# Correlation coefficients in the upper right hand (size scaled to their |r|)
# Biomass & Density most strongly correlated (0.80)


# Saving plot of variable correlations for kelp forest structure
pdf(file="./MSc_plots/SuppFigs/ord_var_corrs_str.pdf", height=8, width=10)
pairs.panels(preds_ord[,c(3:5,9)], scale=T)
dev.off()
# Saving plot of variable correlations for environmental
pdf(file="./MSc_plots/SuppFigs/ord_var_corrs_env.pdf", height=8, width=10)
pairs.panels(preds_ord[,c(2,11,14,16,18,19)], scale=T)
dev.off()


#### All species: Plots
#### All species: Plots
#### Pelagic fish: Generating models, stats ----

## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))

# # Running PERMANOVA to test for initial grouping factors
# pfish_perm <- adonis2(pfish_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom not significant (marginally)


### CAP (kelp spp)
pfish_1 <- capscale(pfish_hel ~ Kelpdom, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray")
anova.cca(pfish_1) # Model is sig
RsquareAdj(pfish_1)
# Adjusted R squared of this model is: 11.6%

## SIMPER for dissimilarity
# Running SIMPER to test between groups
pfish_simp <- with(preds_ord, simper(pfish_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
pfish_simp_df <- as.data.frame(pfish_simp$Giant_Bull)
sum(pfish_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.6027097
# Partitioned among the 10 pelagic fish spp = 0.06027097 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
pfish_simp_filt <- pfish_simp_df %>%
  filter(average > (2*0.06027097))

pfish_simpSPP <- "Embiotoca lateralis"


### RDA (kelp forest structure)

# Generating the kelp forest structure model
pfish_2 <- capscale(pfish_hel ~ Kelpdom+DensityM+Area_m2+BiomassM+HeightM, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray") 
RsquareAdj(pfish_2)
## The structural and environmental variables explain 36.2% of variation (constrained axes)
# Adjusted R squared of this model is: 19.2%
vif.cca(pfish_2) # low vif scores

# Exploring significance
anova.cca(pfish_2, permutations=999) # significance of model
anova.cca(pfish_2, by = "margin", permutations=999) # significance of predictor terms
anova.cca(pfish_2, by = "axis", permutations=999) # significance by model axes


### RDA (environmental vars)

# Generating the environmental model
pfish_3 <- capscale(pfish_hel ~ Kelpdom+Tempave+Depth_datum_m+exp_36+Phardbottom+Psoftbottom+Punderstory, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray")
RsquareAdj(pfish_3)
## Only the environmental variables explain explain 40.7% of variation (constrained axes)
# Adjusted R squared of this model is: 13.7%
vif.cca(pfish_2)

# Exploring significance
anova.cca(pfish_3, permutations=999) # significance of model
anova.cca(pfish_3, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(pfish_3, permutations = 999, by = "axis") # significance by model axes

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

# Extracting the axes data from the model
pfish1_axis.long <- axis.long(pfish_1, choices=c(1,2))
pfish1_axis.long
# Extracting the locations of sites from the model
pfish1_sites.long <- sites.long(pfish1_ordplot, env.data=preds_ord)
head(pfish1_sites.long)
# Extracting the locations of centroids from the sites.long output
pfish1_centroids.long <- centroids.long(pfish1_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish1_centroids.long)

hull_cyl <- pfish1_sites.long %>%
  group_by(Kelpdom) %>%
  slice(chull(axis1, axis2))

# The CAP plot
pfish1_ordplot <- ggplot(pfish1_sites.long) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (14.6%)") +
  ylab("MDS 1 (25.1%)") +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
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
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "a", size=5, fontface=2)
pfish1_ordplot


## Saving plot

tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ord_CAP.tiff", height = 4, width = 5, units = "in", res=400)

pfish1_ordplot

dev.off()



### SIMPER PLOT (to go with CAP plot)

pfish_simp_arr <- pfish_simp_filt %>%
  arrange(desc(average)) %>%
  # mutate(ava = ava*(-1)) %>% # Making the macro relative ab column neg for plotting purposes
  dplyr::select(species, average, sd, ava, avb) %>%
  pivot_longer(cols=c(ava, avb)) %>% # Pivoting the ava (macro) & avb (nereo) columns into one
  ungroup() %>%
  as.data.frame() %>%
  mutate(name = as.factor(name), species = as.factor(species))


# pfish_simplot <- 
ggplot(data=pfish_simp_arr, aes(x=species, y=value, group=name)) +
  geom_errorbar(aes(ymin=value+sd, ymax=(value+sd), group=name), position=position_dodge(width=0.2)) +
  geom_point(aes(color=name), position=position_dodge(0.5)) +
  theme_classic()


ggplot(data=pfish_simp_arr, aes(x=species, y=value, group=name)) +
  geom_pointrange(aes(ymin=(value-sd), ymax=(value+sd), color=name), position=position_dodge(0.3))



### Plot 2: RDA (kelp forest structure)

# Quick ordiplot
pfish2_ordplot <- ordiplot(pfish_2, type="text", display="all")
ordisymbol(pfish2_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(pfish_2, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Prepping for the ggplot biplots

# Extracting the axes data from the model
pfish2_axis.long <- axis.long(pfish_2, choices=c(1,2))
pfish2_axis.long
# Extracting the locations of species from the model
pfish2_species.long <- species.long(pfish2_ordplot)
head(pfish2_species.long)
# Extracting the locations of sites from the model
pfish2_sites.long <- sites.long(pfish2_ordplot, env.data=preds_ord)
head(pfish2_sites.long)
# Extracting the locations of centroids from the sites.long output
pfish2_centroids.long <- centroids.long(pfish2_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish2_centroids.long)
# Extracting env biplot correlations from the model
pfish2_envfit <- envfit(ord=pfish2_ordplot, env=preds_ord)
# Defining env vectors of interest
pfish2_envvectors <- c("DensityM", "Area_m2", "HeightM", "BiomassM")
pfish2_vectorfit.long <- vectorfit.long(pfish2_envfit) %>%
  filter(vector %in% pfish2_envvectors)
head(pfish2_vectorfit.long)

# Recoding spp fct levels to be shorthand scientific names
pfish2_vectorfit.long <- pfish2_vectorfit.long %>%
  mutate(vector = fct_recode(vector,
                             "Density" = "DensityM",
                             "Area" = "Area_m2",
                             "Height" = "HeightM",
                             "Biomass" = "BiomassM"
                             # "Ave Temp" = "Tempave",
                             # "Depth" = "Depth_datum_m",
                             # "Wave Exp" = "exp_36",
                             # "Understory %" = "Punderstory",
                             # "Hard %" = "Phardbottom",
                             # "Soft %" = "Psoftbottom"
  ))


# Extracting the ellipses from the model and ordiplot
pfish2_ordsimple <- ordiplot(pfish_2)
pfish2_ordellipses <- with(preds_ord, ordiellipse(pfish2_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
pfish2_ellipses <- ordiellipse.long(pfish2_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
pfish2_species.envfit <- envfit(pfish2_ordsimple, env=pfish_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
pfish2_species.envfit.table <- data.frame(r=pfish2_species.envfit$vectors$r, p=pfish2_species.envfit$vectors$pvals)
pfish2_species.long.var <- species.long(pfish2_ordsimple, spec.data=pfish2_species.envfit.table)
pfish2_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
pfish2_species.long.vargrt <- pfish2_species.long.var[pfish2_species.long.var$r >= 0.5, ]
pfish2_species.long.vargrt

# # Selecting for SIMPER spp that explain > 70% of the total observed variation among groups
# binv_species_SIMPER <- binv_species.long.var %>%
#   filter(labels %in% binv_simpSPP)
# binv_species_SIMPER

# Recoding spp fct levels to be shorthand scientific names
pfish2_species.long.vargrt <- pfish2_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "C. aggregata" = "Cymatogaster aggregata",
                             "E. lateralis" = "Embiotoca lateralis",
                             "Sebastes" = "Sebastes spp."))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = pfish2_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full structural plot 
pfish2_ordplot_str <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab(pfish2_axis.long[1, "label"]) +
  ylab(pfish2_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=pfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish2_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2.5, yend=axis2*2.5),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=pfish2_vectorfit.long, # Adding in the biplot correlation labels
            aes(x=axis1*3.15, y=axis2*2.9, label=vector),
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
  annotate("text", x = -3.5, y = 3.5, label = "b", size=5, fontface=2)
pfish2_ordplot_str


## The biplot for influential spp

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
  xlab(pfish2_axis.long[1, "label"]) +
  ylab(pfish2_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=pfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish2_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*3.5, yend=axis2*3.5),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=pfish2_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*3.5, y=axis2*3.5, label=labels),
                  colour="black", fontface="italic", size=2, angle=0, nudge_x=0.35, nudge_y=0.1, box.padding=0.5, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "c", size=5, fontface=2)
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

# Extracting the axes data from the model
pfish3_axis.long <- axis.long(pfish_3, choices=c(1,2))
pfish3_axis.long
# Extracting the locations of species from the model
pfish3_species.long <- species.long(pfish3_ordplot)
head(pfish3_species.long)
# Extracting the locations of sites from the model
pfish3_sites.long <- sites.long(pfish3_ordplot, env.data=preds_ord)
head(pfish3_sites.long)
# Extracting the locations of centroids from the sites.long output
pfish3_centroids.long <- centroids.long(pfish3_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish3_centroids.long)
# Extracting env biplot correlations from the model
pfish3_envfit <- envfit(ord=pfish3_ordplot, env=preds_ord)
# Defining env vectors of interest
pfish3_envvectors <- c("Tempave", "Depth_datum_m", "exp_36", "Punderstory", "Phardbottom", "Psoftbottom")
pfish3_vectorfit.long <- vectorfit.long(pfish3_envfit) %>%
  filter(vector %in% pfish3_envvectors)
head(pfish3_vectorfit.long)

# Recoding spp fct levels to be shorthand scientific names
pfish3_vectorfit.long <- pfish3_vectorfit.long %>%
  mutate(vector = fct_recode(vector,
                             "Ave Temp" = "Tempave",
                             "Depth" = "Depth_datum_m",
                             "Wave Exp" = "exp_36",
                             "Understory %" = "Punderstory",
                             "Hard %" = "Phardbottom",
                             "Soft %" = "Psoftbottom"
  ))


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

# # Selecting for SIMPER spp that explain > 70% of the total observed variation among groups
# pfish_species_SIMPER <- pfish_species.long.var %>%
#   filter(labels %in% pfish_simpSPP)
# pfish_species_SIMPER

# Recoding spp fct levels to be shorthand scientific names
pfish3_species.long.vargrt <- pfish3_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "B. frenatus" = "Brachyistius frenatus",
                             "C. pallasii" = "Clupea pallasii",
                             "E. lateralis" = "Embiotoca lateralis",
                             "S. caurinus" = "Sebastes caurinus"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = pfish3_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
pfish3_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab(pfish3_axis.long[1, "label"]) + 
  ylab(pfish3_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=pfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish3_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2.5, yend=axis2*2.5),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=pfish3_vectorfit.long, # Adding in the biplot correlation labels
            aes(x=axis1*3.25, y=axis2*2.9, label=vector),
            colour="black", size=2, hjust=0.4, vjust=0.6) +
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
  annotate("text", x = -3.5, y = 3.5, label = "d", size=5, fontface=2)
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
  xlab(pfish3_axis.long[1, "label"]) +
  ylab(pfish3_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=pfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=pfish3_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*3.5, yend=axis2*3.5),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=pfish3_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*3.5, y=axis2*3.5, label=labels),
                  colour="black", fontface="italic", size=2, angle=0, box.padding=0.4, nudge_y=-0.9, nudge_x=-0.05, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "e", size=5, fontface=2)
pfish3_ordplot_spp


## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ord_RDAenv_draft.tiff", height = 3, width = 7, units = "in", res=400)

pfish3_ordplots <- ggarrange(pfish3_ordplot_env, pfish3_ordplot_spp, ncol=2, common.legend=FALSE)
pfish3_ordplots

dev.off()

#


### FINAL ARRANGEMENT

# Setting up the layout
lay <- rbind(c(1,1,2,2,2,2),
             c(1,1,2,2,2,2),
             c(1,1,3,3,3,3),
             c(4,5,3,3,3,3))
# Rectangles to help visualize
gs <- lapply(1:5, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)

# Saving as an image object
pdf(file="./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined.pdf", height=5.5, width=10)
tiff(file="./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined.tiff", height=12, width=21, units="cm", res=300)

# End ggplot arrangement
grid <- grid.arrange(pfish1_ordplot, pfish2_ordplots, pfish3_ordplots,
             layout_matrix = lay)

dev.off()


# Drawing in the animal outlines

# Loading spp image outlines
rfish <- image_read('./Phylopic/Fish/rockfish.png')
sperch <- image_read('./Phylopic/Fish/surfperch.png')

# Calling map plot as magick image
mplot <- image_read("./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined.tiff") 

# Adding in the animal outlines (scaling & location)
mmap <- image_composite(mplot, image_scale(rfish, "x100"), offset="+200+1100")
mmap <- image_composite(mmap, image_scale(sperch, "x120"), offset="+500+1200")

# Saving the final fig
image_write(mmap, path = "./MSc_plots/PaperFigs/Pfish/Pfish_ords_combined_outlines.tiff", format = "tiff", density=400)

#

#### Demersal fish: Generating ordination models, stats ----

## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))

# # Running PERMANOVA to test for initial grouping factors
# bfish_perm <- adonis2(bfish_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray")  # Kelpdom significant (0.006)


### CAP (kelp spp)
bfish_1 <- capscale(bfish_hel ~ Kelpdom, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray")
anova.cca(bfish_1) # model is sig
RsquareAdj(bfish_1)
## Only the environmental variables explain 22.1% of variation (constrained axes)
# Adjusted R squared of this model is: 20.7%

## SIMPER for dissimilarity
# Running SIMPER to test between groups (limiting to 10 most influential)
bfish_simp <- with(preds_ord, simper(bfish_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
bfish_simp_df <- as.data.frame(summary(bfish_simp)$Giant_Bull)
sum(bfish_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.4860742
# Partitioned among the 10 demersal fish spp = 0.04860742 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
bfish_simp_filt <- bfish_simp_df %>%
  filter(average > (2*0.04860742))

bfish_simpSPP <- c("Artedius harringtoni", "Rhinogobiops nicholsii")


### RDA (kelp forest structure)

# Generating the kelp forest structure model
bfish_2 <- capscale(bfish_hel ~ Kelpdom+DensityM+Area_m2+BiomassM+HeightM, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray") 
RsquareAdj(bfish_2)
## Only the environmental variables explain 43.8% of variation (constrained axes)
# Adjusted R squared of this model is: 31.0%
vif.cca(bfish_2) # low vif scores

# Exploring model significance
anova.cca(bfish_2) # significance of model
anova.cca(bfish_2, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(bfish_2, permutations = 999, by = "axis") # significance by model axes


### RDA (environmental vars)

# Generating the environmental model
bfish_3 <- capscale(bfish_hel ~ Kelpdom+Tempave+Depth_datum_m+exp_36+Phardbottom+Psoftbottom+Punderstory, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray")
RsquareAdj(bfish_3)
## The structural and environmental variables explain explain 43.8% of variation (constrained axes)
# Adjusted R squared of this model is: 45.3%
vif.cca(bfish_3) # low vif scores

# Exploring model significance
anova.cca(bfish_3) # significance of model
anova.cca(bfish_3, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(bfish_3, permutations = 999, by = "axis") # significance by model axes


#### Demersal fish: Plots ----

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

# Extracting the axes data from the model
bfish1_axis.long <- axis.long(bfish_1, choices=c(1,2))
bfish1_axis.long
# Extracting the locations of sites from the model
bfish1_sites.long <- sites.long(bfish1_ordplot, env.data=preds_ord)
head(bfish1_sites.long)
# Extracting the locations of centroids from the sites.long output
bfish1_centroids.long <- centroids.long(bfish1_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(bfish1_centroids.long)

hull_cyl <- bfish1_sites.long %>%
  group_by(Kelpdom) %>%
  slice(chull(axis1, axis2))

# The CAP plot
bfish1_ordplot <- ggplot(bfish1_sites.long) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (19.9%)") +
  ylab("MDS 1 (4.4%)") +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  geom_point(data=bfish1_sites.long, # Adding in the main site data points
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
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "a", size=5, fontface=2)
bfish1_ordplot


## Saving plot

tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ord_CAP.tiff", height = 4, width = 5, units = "in", res=400)

bfish1_ordplot

dev.off()



### SIMPER PLOT (to go with CAP plot)

bfish_simp_arr <- bfish_simp_filt %>%
  arrange(desc(average)) %>%
  # mutate(ava = ava*(-1)) %>% # Making the macro relative ab column neg for plotting purposes
  dplyr::select(species, average, sd, ava, avb) %>%
  pivot_longer(cols=c(ava, avb)) %>% # Pivoting the ava (macro) & avb (nereo) columns into one
  ungroup() %>%
  as.data.frame() %>%
  mutate(name = as.factor(name), species = as.factor(species))


# pfish_simplot <- 
ggplot(data=bfish_simp_arr, aes(x=species, y=value, group=name)) +
  geom_errorbar(aes(ymin=value+sd, ymax=(value+sd), group=name), position=position_dodge(width=0.2)) +
  geom_point(aes(color=name), position=position_dodge(0.5)) +
  theme_classic()


ggplot(data=bfish_simp_arr, aes(x=species, y=value, group=name)) +
  geom_pointrange(aes(ymin=(value-sd), ymax=(value+sd), color=name), position=position_dodge(0.3))



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
# Extracting env biplot correlations from the model
bfish2_envfit <- envfit(ord=bfish2_ordplot, env=preds_ord)
# Defining env vectors of interest
bfish2_envvectors <- c("DensityM", "Area_m2", "HeightM", "BiomassM")
bfish2_vectorfit.long <- vectorfit.long(bfish2_envfit) %>%
  filter(vector %in% bfish2_envvectors)
head(bfish2_vectorfit.long)

# Recoding spp fct levels to be shorthand scientific names
bfish2_vectorfit.long <- bfish2_vectorfit.long %>%
  mutate(vector = fct_recode(vector,
                             "Density" = "DensityM",
                             "Area" = "Area_m2",
                             "Height" = "HeightM",
                             "Biomass" = "BiomassM"
                             # "Ave Temp" = "Tempave",
                             # "Depth" = "Depth_datum_m",
                             # "Wave Exp" = "exp_36",
                             # "Understory %" = "Punderstory",
                             # "Hard %" = "Phardbottom",
                             # "Soft %" = "Psoftbottom"
  ))


# # Extracting the ellipses from the model and ordiplot
# bfish2_ordsimple <- ordiplot(bfish_2)
# bfish2_ordellipses <- with(preds_ord, ordiellipse(bfish2_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
# bfish2_ellipses <- ordiellipse.long(bfish2_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
bfish2_species.envfit <- envfit(bfish2_ordsimple, env=bfish_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
bfish2_species.envfit.table <- data.frame(r=bfish2_species.envfit$vectors$r, p=bfish2_species.envfit$vectors$pvals)
bfish2_species.long.var <- species.long(bfish2_ordsimple, spec.data=bfish2_species.envfit.table)
bfish2_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
bfish2_species.long.vargrt <- bfish2_species.long.var[bfish2_species.long.var$r >= 0.5, ]
bfish2_species.long.vargrt

# # Selecting for SIMPER spp that explain > 70% of the total observed variation among groups
# binv_species_SIMPER <- binv_species.long.var %>%
#   filter(labels %in% binv_simpSPP)
# binv_species_SIMPER

# Recoding spp fct levels to be shorthand scientific names
bfish2_species.long.vargrt <- bfish2_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "A. harringtoni" = "Artedius harringtoni",
                             "H. decagrammos" = "Hexagrammos decagrammus",
                             "J. zonope" = "Jordania zonope",
                             "R. nicholsii" = "Rhinogobiops nicholsii"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = bfish2_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full structural plot 
bfish2_ordplot_str <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (35.6%)") + # CAP 1
  ylab("CAP 2 (2.9%)") + # CAP 2
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=bfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish2_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
            aes(x=xnew*2.4, y=ynew*2.3, label=vector),
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
  annotate("text", x = -3.5, y = 3.5, label = "b", size=5, fontface=2)
bfish2_ordplot_str


## The biplot for influential spp

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
  xlab("CAP 1 (35.6%)") + # CAP 1
  ylab("CAP 2 (2.9%)") + # CAP 2
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=bfish2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish2_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=bfish2_species.long.vargrt, # Adding in the SIMPER species labels
            aes(x=axis1*2, y=axis2*2, label=labels),
            colour="black", fontface="italic", size=2, angle=0, box.padding=0.2, nudge_x=1.1, nudge_y=-0.8, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "c", size=5, fontface=2)
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
# Extracting env biplot correlations from the model
bfish3_envfit <- envfit(ord=bfish3_ordplot, env=preds_ord)
# Defining env vectors of interest
bfish3_envvectors <- c("Tempave", "Depth_datum_m", "exp_36", "Punderstory", "Phardbottom", "Psoftbottom")
bfish3_vectorfit.long <- vectorfit.long(bfish3_envfit) %>%
  filter(vector %in% bfish3_envvectors)
head(bfish3_vectorfit.long)

# Recoding spp fct levels to be shorthand scientific names
bfish3_vectorfit.long <- bfish3_vectorfit.long %>%
  mutate(vector = fct_recode(vector,
                             "Ave Temp" = "Tempave",
                             "Depth" = "Depth_datum_m",
                             "Wave Exp" = "exp_36",
                             "Understory %" = "Punderstory",
                             "Hard %" = "Phardbottom",
                             "Soft %" = "Psoftbottom"))


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
                             "A. harringtoni" = "Artedius harringtoni",
                             "H. decagrammos" = "Hexagrammos decagrammus",
                             "J. zonope" = "Jordania zonope",
                             "R. nicholsii" = "Rhinogobiops nicholsii"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = bfish3_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
bfish3_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (40.1%)") + # CAP 1
  ylab("CAP 2 (7.9%)") + # CAP 2
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=bfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish3_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2.5, yend=axis2*2.5),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=bfish3_vectorfit.long, # Adding in the biplot correlation labels
            aes(x=axis1*2.5, y=axis2*2.5, label=vector),
            colour="black", size=2, hjust=0, vjust=0, nudge_y=-0.65, nudge_x=-0.22, box.padding=0.93, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "d", size=5, fontface=2)
bfish3_ordplot_env


## The biplot for influential spp

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
  xlab("CAP 1 (40.1%)") + # CAP 1
  ylab("CAP 2 (7.9%)") + # CAP 2  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=bfish3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=bfish3_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=bfish3_species.long.vargrt, # Adding in the SIMPER species labels
            aes(x=axis1*2, y=axis2*2, label=labels),
            colour="black", fontface="italic", size=2, angle=0, nudge_x=0.8, nudge_y=-0.2, box.padding=0.6, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "e", size=5, fontface=2)
bfish3_ordplot_spp


## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ord_RDAenv_draft.tiff", height = 3.5, width = 7, units = "in", res=400)

bfish3_ordplots <- ggarrange(bfish3_ordplot_env, bfish3_ordplot_spp, ncol=2, common.legend=FALSE)
bfish3_ordplots

dev.off()

#


### FINAL ARRANGEMENT

# Setting up the layout
lay <- rbind(c(1,1,2,2,2,2),
             c(1,1,2,2,2,2),
             c(1,1,3,3,3,3),
             c(4,5,3,3,3,3))
# Rectangles to help visualize
gs <- lapply(1:5, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)


tiff(file="./MSc_plots/PaperFigs/Bfish/Bfish_ords_combined.tiff", height=12, width=21, units="cm", res=300)

# End ggplot arrangement
grid.arrange(bfish1_ordplot, bfish2_ordplots, bfish3_ordplots,
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
mmap <- image_composite(mplot, image_scale(bgoby, "x90"), offset="+200+1130")
mmap <- image_composite(mmap, image_rotate(image_scale(sculp, "x90"), 330), offset="+500+1160")

# Saving the final fig
image_write(mmap, path = "./MSc_plots/PaperFigs/Bfish/Bfish_ords_combined_outlines.tiff", format = "tiff", density=400)

#


#### All fish: Generating ordination models, stats ----

###
#### Models / stats
###

## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom), Cluster = as.factor(Cluster))

# Running PERMANOVA to test for grouping factors
fish_perm <- adonis2(fish_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom significant (0.015)


### CAP (kelp spp)

# Generating the CAP kelp spp model
fish_1 <- capscale(fish_hel ~ Kelpdom, data=preds_ord, comm=fish_wide, add=FALSE, distance="bray") 
RsquareAdj(fish_1)
## The structural and environmental variables explain 13.2% of variation (constrained axes)
# Adjusted R squared of this model is: 10.3%

## SIMPER for dissimilarity 
# Running SIMPER to test between groups (limiting to 10 most influential)
fish_simp <- with(preds_ord, simper(fish_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
fish_simp_df <- as.data.frame(fish_simp$Giant_Bull)
sum(fish_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.5705693
# Partitioned among the 19 fish spp = 0.03002996 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
fish_simp_filt <- fish_simp_df %>%
  filter(average > (2*0.03002996))

fish_simpSPP <- c("Embiotoca lateralis", "Rhinogobiops nicholsii")


### RDA (kelp forest structure)

# Generating the kelp forest structure model
fish_2 <- capscale(fish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+BiomassM, data=preds_ord, comm=fish_wide, add=FALSE, distance="bray") 
RsquareAdj(fish_2)
## The structural and environmental variables explain 29.4% of variation (constrained axes)
# Adjusted R squared of this model is: 11.2%
vif.cca(fish_2) # low vif scores

# Exploring model significance
anova.cca(fish_2) # significance of model
anova.cca(fish_2, step = 1000, by = "term") # significance of predictor terms
anova.cca(fish_2, step = 1000, by = "axis") # significance by model axes


### RDA (environmental vars)

# Generating the environmental model
fish_3 <- capscale(fish_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=fish_wide, add=FALSE, distance="bray")
RsquareAdj(fish_3)
## Only the environmental variables explain explain 39.3% of variation (constrained axes)
# Adjusted R squared of this model is: 14.8%
vif.cca(fish_3) # low vif scores

# Exploring model significance
anova.cca(fish_3) # significance of model
anova.cca(fish_3, step = 1000, by = "term") # significance of predictor terms
anova.cca(fish_3, step = 1000, by = "axis") # significance by model axes

#### All fish: Plots
#### Benthic inverts: Generating ordination models, stats ----

## Calling the predictor variables sheet
preds_ord <- read_csv("./MSc_data/Data_new/AllPredictors_2022.csv") %>%
  mutate(Kelpdom = as.factor(Kelpdom))

# # Running PERMANOVA to test for grouping factors
# binv_perm <- adonis2(binv_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom & Cluster significant

### CAP (kelp spp)

# Generating the CAP kelp spp model
binv_1 <- capscale(binv_hel ~ Kelpdom, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray") 
anova.cca(binv_1) # model is sig
RsquareAdj(binv_1)
## Only the kelp spp explains 11.7% of variation (constrained axes)
# Adjusted R squared of this model is: 7.3%

## SIMPER for dissimilarity 
# Running SIMPER to test between groups (limiting to 10 most influential)
binv_simp <- with(preds_ord, simper(binv_hel, group=Kelpdom, permutations=999))

# Making table of the SIMPER output
binv_simp_df <- as.data.frame(summary(binv_simp)$Giant_Bull)
sum(binv_simp_df$average) # Total dissimilarity b/w kelp spp groups = 0.4547604
# Partitioned among the 30 invertebrate spp = 0.01515868 (assuming all contribute equally)

# Filter for spp that contribute >= 2x their expected dissimilarity among groups
binv_simp_filt <- binv_simp_df %>%
  filter(average > (2*0.01515868))
binv_simp_filt

binv_simps <- c("Pomaulax gibberosus", "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus")
binv_simpSPP <- c("P. gibberosus", "S. purpuratus", "M. franciscanus")


### RDA (kelp forest structure)

# Generating the kelp forest structure model
binv_2 <- capscale(binv_hel ~ Kelpdom+DensityM+Area_m2+BiomassM+HeightM, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray") 
RsquareAdj(binv_2)
## Only the structural variables explain 45.2% of variation (constrained axes)
# Adjusted R squared of this model is: 28.2%
vif.cca(binv_2) # low vif scores

# Exploring model significance
anova.cca(binv_2) # significance of model
anova.cca(binv_2, step = 1000, by = "margin") # significance of predictor terms
anova.cca(binv_2, step = 1000, by = "axis") # significance by model axes


### RDA (environmental vars)

# Generating the environmental model
binv_3 <- capscale(binv_hel ~ Kelpdom+Tempave+Depth_datum_m+exp_36+Phardbottom+Psoftbottom+Punderstory, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray")
RsquareAdj(binv_3)
## The environmental variables explain explain 54.8% of variation (constrained axes)
# Adjusted R squared of this model is: 32.1%
vif.cca(binv_3)

# Exploring model significance
anova.cca(binv_3) # significance of model
anova.cca(binv_3, permutations = 999, by = "margin") # significance of predictor terms
anova.cca(binv_3, permutations = 999, by = "axis") # significance by model axes


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

hull_cyl <- binv1_sites.long %>%
  group_by(Kelpdom) %>%
  slice(chull(axis1, axis2))

# The CAP plot
binv1_ordplot <- ggplot(binv1_sites.long) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (11.2%)") +
  ylab("MDS 1 (2.6%)") +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
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
        axis.line=element_line(size=0.3),
        axis.ticks=element_line(size=0.3)) +
  annotate("text", x = -3, y = 3, label = "a", size=5, fontface=2)
binv1_ordplot


## Saving plot

tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ord_CAP.tiff", height = 4, width = 5, units = "in", res=400)

binv1_ordplot

dev.off()



### SIMPER PLOT (to go with CAP plot)

binv_simp_arr <- binv_simp_filt %>%
  arrange(desc(average)) %>%
  # mutate(ava = ava*(-1)) %>% # Making the macro relative ab column neg for plotting purposes
  dplyr::select(species, average, sd, ava, avb) %>%
  pivot_longer(cols=c(ava, avb)) %>% # Pivoting the ava (macro) & avb (nereo) columns into one
  ungroup() %>%
  as.data.frame() %>%
  mutate(name = as.factor(name), species = as.factor(species))


# binv1_simplot <- 
  ggplot(data=binv_simp_arr, aes(x=species, y=value, group=name)) +
    geom_errorbar(aes(ymin=value+sd, ymax=(value+sd), group=name), position=position_dodge(width=0.2)) +
    geom_point(aes(color=name), position=position_dodge(0.5)) +
    theme_classic()
  

  ggplot(data=binv_simp_arr, aes(x=species, y=value, group=name)) +
    geom_pointrange(aes(ymin=(value-sd), ymax=(value+sd), color=name), position=position_dodge(0.3))



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
# Extracting env biplot correlations from the model
binv2_envfit <- envfit(ord=binv2_ordplot, env=preds_ord)
# Defining env vectors of interest
binv2_envvectors <- c("DensityM", "Area_m2", "HeightM", "BiomassM")
binv2_vectorfit.long <- vectorfit.long(binv2_envfit) %>%
  filter(vector %in% binv2_envvectors)
head(binv2_vectorfit.long)

# Recoding spp fct levels to be shorthand scientific names
binv2_vectorfit.long <- binv2_vectorfit.long %>%
  mutate(vector = fct_recode(vector,
                             "Density" = "DensityM",
                             "Area" = "Area_m2",
                             "Height" = "HeightM",
                             "Biomass" = "BiomassM"
                             # "Ave Temp" = "Tempave",
                             # "Depth" = "Depth_datum_m",
                             # "Wave Exp" = "exp_36",
                             # "Understory %" = "Punderstory",
                             # "Hard %" = "Phardbottom",
                             # "Soft %" = "Psoftbottom"
         ))


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

# # Selecting for SIMPER spp that explain > 70% of the total observed variation among groups
# binv_species_SIMPER <- binv_species.long.var %>%
#   filter(labels %in% binv_simpSPP)
# binv_species_SIMPER

# Recoding spp fct levels to be shorthand scientific names
binv2_species.long.vargrt <- binv2_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "C. foliatum" = "Ceratostoma foliatum",
                             "D. odonoghuei" = "Diaulula odonoghuei",
                             "H. leviuscula" = "Henricia leviuscula",
                             "H. stylus" = "Heptacarpus stylus",
                             "M. franciscanus" = "Mesocentrotus franciscanus",
                             "O. gracilis" = "Oregonia gracilis",
                             "P. gibberosus" = "Pomaulax gibberosus",
                             "P. helianthoides" = "Pycnopodia helianthoides"))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = binv2_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full structural plot 
binv2_ordplot_str <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (19.6%)") +
  ylab("CAP 2 (12.3%)") +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=binv2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv2_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=binv2_vectorfit.long, # Adding in the biplot correlation labels
            aes(x=axis1*2.7, y=axis2*2.4, label=vector),
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
  annotate("text", x = -3.5, y = 3.5, label = "b", size=5, fontface=2)
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
  xlab("CAP 1 (19.6%)") +
  ylab("CAP 2 (12.3%)") +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=binv2_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv2_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*3, yend=axis2*3),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=binv2_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*3, y=axis2*3, label=labels),
                  colour="black", fontface="italic", size=1.8, angle=0, hjust=0, vjust=0, box.padding=0.3, nudge_y=2.2, nudge_x=0.62, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "c", size=5, fontface=2)
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
# Extracting env biplot correlations from the model
binv3_envfit <- envfit(ord=binv3_ordplot, env=preds_ord)
# Defining env vectors of interest
binv3_envvectors <- c("Tempave", "Depth_datum_m", "exp_36", "Punderstory", "Phardbottom", "Psoftbottom")
binv3_vectorfit.long <- vectorfit.long(binv3_envfit) %>%
  filter(vector %in% binv3_envvectors)
head(binv3_vectorfit.long)

# Recoding spp fct levels to be shorthand scientific names
binv3_vectorfit.long <- binv3_vectorfit.long %>%
  mutate(vector = fct_recode(vector,
                             "Ave Temp" = "Tempave",
                             "Depth" = "Depth_datum_m",
                             "Wave Exp" = "exp_36",
                             "Understory %" = "Punderstory",
                             "Hard %" = "Phardbottom",
                             "Soft %" = "Psoftbottom"
  ))


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

# # Selecting for SIMPER spp that explain > 70% of the total observed variation among groups
# binv_species_SIMPER <- binv_species.long.var %>%
#   filter(labels %in% binv_simpSPP)
# binv_species_SIMPER

# Recoding spp fct levels to be shorthand scientific names
binv3_species.long.vargrt <- binv3_species.long.vargrt %>%
  mutate(labels = fct_recode(labels,
                             "C. foliatum" = "Ceratostoma foliatum",
                             "D. odonoghuei" = "Diaulula odonoghuei",
                             "H. leviuscula" = "Henricia leviuscula",
                             "H. stylus" = "Heptacarpus stylus",
                             "M. franciscanus" = "Mesocentrotus franciscanus",
                             "O. gracilis" = "Oregonia gracilis",
                             "P. gibberosus" = "Pomaulax gibberosus",
                             "P. helianthoides" = "Pycnopodia helianthoides",
                             "E. troschelii" = "Evasterias troschelii",
                             "C. productus" = "Cancer productus"
                             ))


## The biplot for structure vars

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = binv3_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
binv3_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 5, linewidth=0.3) +  
  xlab("CAP 1 (18.9%)") +
  ylab("CAP 2 (16.0%)") +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=binv3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv3_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.3, arrow=arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text(data=binv3_vectorfit.long, # Adding in the biplot correlation labels
            aes(x=axis1*2.7, y=axis2*2.3, label=vector),
            colour="black", size=2) +
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
  annotate("text", x = -3.5, y = 3.5, label = "d", size=5, fontface=2)
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
  xlab("CAP 1 (18.9%)") +
  ylab("CAP 2 (16.0%)") +  
  scale_y_continuous(limits=c(-3.5,3.5)) + 
  scale_x_continuous(limits=c(-3.5,3.5)) + 
  geom_point(data=binv3_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=2.5, shape=23, stroke=0.2, alpha=0.8) +
  geom_segment(data=binv3_species.long.vargrt, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*3, yend=axis2*3),
               colour="grey20", linewidth=0.3, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
  geom_text_repel(data=binv3_species.long.vargrt, # Adding in the SIMPER species labels
                  aes(x=axis1*3, y=axis2*3, label=labels),
                  colour="black", fontface="italic", size=1.8, angle=0, hjust=0, vjust=0, box.padding=0.33, nudge_y=-2.4, nudge_x=1.6, segment.linetype=2, segment.size=0.3, segment.color="grey40") +
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
  annotate("text", x = -3.5, y = 3.5, label = "e", size=5, fontface=2)
binv3_ordplot_spp



## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ord_RDAenv_draft.tiff", height = 3.5, width = 7, units = "in", res=400)

binv3_ordplots <- ggarrange(binv3_ordplot_env, binv3_ordplot_spp, ncol=2, common.legend=FALSE)
binv3_ordplots

dev.off()

#


### FINAL ARRANGEMENT

# Setting up the layout
lay <- rbind(c(1,1,2,2,2,2),
             c(1,1,2,2,2,2),
             c(1,1,3,3,3,3),
             c(4,5,3,3,3,3))
# Rectangles to help visualize
gs <- lapply(1:5, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
grid.arrange(grobs=gs,
             layout_matrix = lay)


tiff(file="./MSc_plots/PaperFigs/Invert/Invert_ords_combined.tiff", height=12, width=21, units="cm", res=300)

# End ggplot arrangement
grid.arrange(binv1_ordplot, binv2_ordplots, binv3_ordplots,
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
mmap <- image_composite(mplot, image_scale(purchin, "x230"), offset="+80+1020")
mmap <- image_composite(mmap, image_rotate(image_scale(turban, "x160"), 330), offset="+500+1070")
mmap <- image_composite(mmap, image_scale(bstar, "x200"), offset="+280+1160")

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


