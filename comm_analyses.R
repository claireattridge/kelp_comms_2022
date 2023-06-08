## Barkley Sound community analyses
## Author: Claire Attridge
## Origin date: Feb 2023

# Loading base packages
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(MetBrewer)
library(stats)
library(car)
library(factoextra)
library(wesanderson)
library(vegan)
library(cowplot)

library(funrar)
library(ecodist)
library(BiodiversityR)
library(ggforce)
library(ggrepel)
library(pairwiseAdonis)

#### Loading all data sheets ----

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
kelps <- read_csv("./MSc_data/Data_new/kelpmetrics_2022.csv") %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  dplyr::select(-c(HeightSD, BiomassSD, DensitySD, MacroSD, NereoSD))
kelps <- kelps %>%
  mutate(MacroP = MacroM/DensityM) # Creating a column for proportion Macro at a given site


# Listing the variables of interest together
dfs1 <- list(kelps, temps, rei, subs)

# Joining the variables into one dataframe
preds <- reduce(dfs1, dplyr::left_join, by="SiteName")
preds[is.na(preds)] <- 0 # Translating the NAs to zeros for clustering analysis

preds <- preds[-17,] # Removing second beach south as an outlier site 

# Polishing dataframe with environmental and ecological data to use for later use in ordination
preds_ord <- preds %>%
  mutate(Cluster = case_when(SiteName == "Between Scotts and Bradys" | SiteName == "Danvers Danger Rock" | SiteName == "Dixon Island Back (Bay)" | SiteName == "Dodger Channel 1" | SiteName == "Flemming 112" | SiteName == "Flemming 114" | SiteName == "North Helby Rock" | SiteName == "Tzartus 116" ~ "C2",
                             SiteName == "Ed King East Inside" | SiteName == "Ross Islet 2" | SiteName == "Ross Islet Slug Island" | SiteName == "Taylor Rock" | SiteName == "Turf Island 2" | SiteName == "Wizard Islet South" ~ "C1",
                             SiteName == "Bordelais Island" | SiteName == "Second Beach" | SiteName == "Dodger Channel 2" | SiteName == "Nanat Bay" ~ "C3",
                             SiteName == "Cable Beach (Blow Hole)" | SiteName == "Less Dangerous Bay" | SiteName == "Swiss Boy" | SiteName == "Wizard Islet North" ~ "C4",
                             TRUE ~ "Aux")) %>%
  mutate(Cluster = as.factor(Cluster), Composition = as.factor(Composition)) %>%
  mutate(Kelpdom = case_when(Composition == "Macro" | Composition == "Mixed" ~ "Giant", # Macro is dominant at mixed sites (too few mixed sites to statistically compare)
                             Composition == "Nereo" ~ "Bull",
                             Composition == "None" ~ "Giant")) %>% # Re. Starko (2022) - used to have M. pyrifera present around the area of Less Dangerous bay
  mutate(Kelpdom = as.factor(Kelpdom)) %>%
  droplevels()

# saving a .csv file of the ecological and environmental predictors of interest for our community abundance data
write.csv(preds_ord, "./MSc_data/Data_new/AllPredictors_2022.csv", row.names=F)


#### Cluster analysis (kelp forest structure) ----

# Checking for collinearity of variables
mod <- lm(SpeciesAb ~ BiomassM + DensityM + HeightM + Area_m2 + MacroP, data=comms)
vif(mod)
# Biomass is the highly correlated with other kelp predictors in my dataset
# Once removed the VIF scores of remaining predictors are within a reasonable range

plot(BiomassM ~ DensityM, data=preds) # Appears most correlated with density
plot(BiomassM ~ HeightM, data=preds) 
plot(BiomassM ~ Area_m2, data=preds)
plot(BiomassM ~ MacroP, data=preds)

# Plotting out the relationships to look for correlations
mod <- lm(HeightM ~ Depth_datum_m, data=preds)
modpred <- predict(mod)  
plot(HeightM ~ Depth_datum_m, data=preds)
lines(preds$Depth_datum_m, modpred)

## Density & Biomass -> very correlated
## Height & Depth -> very correlated


# ## PCA (for density & biomass)
# pca <- prcomp(na.omit(preds[,c(3:4)]), center = TRUE, scale = TRUE)
# fviz_eig(pca) # Scree plot: PC1 explains vast majority of variation (89.9%)
# fviz_pca_biplot(pca) # Biplot: Pos. PC1 values = higher density & biomass
# 
# pc1 <- as.vector(pca$x[,1]) # Extracting the pc1 values for further use
# 
# preds$PC1 <- pc1 # Adding pc1 values to the working dataframe
  

## K MEANS

preds <- preds %>% # Scaling all predictors prior to clustering (Euc dists are sens to scale)
  mutate_at(c(2,3,4,5,6,7,9,10), funs(c(scale(.))))


# optimizing for k cluster values
rng <- 2:15 # K from 2 to 15
tries <- 100 #Run the K Means algorithm 100 times
avg.totw.ss <-integer(length(rng)) #Set up an empty vector to hold all of points

for(v in rng){ # For each value of the range variable
  v.totw.ss <- integer(tries) # Set up an empty vector to hold the 100 tries
  for(i in 1:tries){
    k.temp <- kmeans(na.omit(preds[,c(3,5,9,10)]), centers=v) # Run kmeans
    v.totw.ss[i] <-k.temp$tot.withinss # Store the total withinss
  }
  avg.totw.ss[v-1] <-mean(v.totw.ss) #Average the 100 total withinss
}

# plot to show ave within SS as K increases (when is within SS reduced as we ^ clusters)
plot(rng,avg.totw.ss,type="b", main="Total Within SS by Various K",
     ylab="Average Total Within Sum of Squares",
     xlab="Value of K")
# looks like 4 clusters will substantially minimize error within clusters


# only using the kelp forest attributes for clustering sites
k <- kmeans(na.omit(preds[,c(3,5,9,10)]), centers=4, nstart=20)
k$centers # looking at the centroid mean values
table(k$cluster) # looking at # sites per


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
# 
# imgras <- as.raster(imgstk) 
# 
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

# Filtering out the site redo surveys at Ross Islet 2, Less Dangerous Bay, & Second Beach S in Sept 2022
# Also removing the acoustic baseline site of 'Sand Town'
comms <- comms %>%
  filter(Date <= as.POSIXct("2022-09-08") & SiteName != "Sand Town") %>%
  droplevels()
  

## Base done, ready to go!

## Let's look into some questions:

qs <- comms %>%
  filter(Method==1) %>%
  droplevels()

# Which site has the most species?
richness <- comms %>%
  group_by(SiteName) %>%
  summarize(spp_rich = n_distinct(Species))%>%
  arrange(desc(spp_rich))

# Which species were most abundant overall?
abun <- comms %>%
  group_by(common_name) %>%
  summarize(totalab = sum(abundance)) %>%
  arrange(desc(totalab))


# Which species had the most site observations?
obvs <- comms %>%
  group_by(common_name) %>%
  summarize(observs = length(common_name)) %>% # The number of rows (observations) per spp common name
  arrange(desc(observs)) %>% # Arranging in descending order
  ungroup() %>%
  as.data.frame()


# Plot for species observations - Supp Fig
pdf(file="./MSc_plots/SuppFigs/SpeciesAbCutoff.pdf", height=10, width=7.5)

obvsplot <- ggplot(data=obvs, aes(x=observs, y=(fct_reorder(common_name, desc(observs))), fill=observs)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_gradientn(colors=met.brewer("VanGogh3")) +
  theme_classic() +
  theme(legend.position="none",
    axis.text.x = element_text(color="black", size="10"),
    axis.text.y = element_text(color="black", size="7"),
    axis.title.y = element_text(color="black", size="11", vjust=0.5),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"), 
    ) +
  xlab("Unique sightings") + ylab("") +
  geom_hline(yintercept = 56.5, linetype=2) # Adding in the cutoff line for =< 5 obvs
obvsplot

dev.off()


## With this ^^ in mind: Removing uncommon species (=< 5 site obvs) from the analysis frame (n = 37 removed)
sppkeep <- obvs %>%
  filter(observs >= 5)

# Remaining species (n = 55 retained)
commsclean <- comms %>%
  filter(common_name %in% sppkeep$common_name)


# Now: What are the total counts for all species at all sites?
comms_all <- commsclean %>%
  as.data.frame() %>%
  ungroup() %>%
  group_by(SiteName, Species, common_name, Method) %>%
  summarise(abundance = sum(abundance))


# saving a .csv file of the clean total taxa abundances by site
write.csv(comms_all, "./MSc_data/Data_new/RLS_2022.csv", row.names=F)


#### Community data: Species groupings ----

# Loading the RLS comms file
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  dplyr::select(SiteName, common_name, abundance, Method) %>%
  rename(TaxaAb = abundance) %>%
  rename(CommonName = common_name) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName), CommonName = as.factor(CommonName)) %>%
  filter(SiteName != "Second Beach South")

## For all species in community together
taxa_all <- taxa %>%
  dplyr::select(-Method) %>% # Removing RLS method col when using all species together
  group_by(SiteName, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb)) %>%
  ungroup() %>%
  as.data.frame()
# Cleaning up the all spp. frame
taxa_all <- taxa_all %>%
  mutate(CommonName = fct_recode(CommonName, "Green sea urchin" = "Northern sea urchin"),
         CommonName = fct_recode(CommonName, "Blood star" = "Unidentified blood star"),
         CommonName = fct_recode(CommonName, "Chiton" = "Unidentified chiton"),
         CommonName = fct_recode(CommonName, "Hermit crab" = "Unidentified hermit crab")) %>%
  filter(CommonName != "Unidentified rockfish") %>% # Removing unidentified rockfish as a spp.
  droplevels()

## For just RLS Method 1 species (pelagic fish)
taxap <- taxa %>%
  filter(Method == 1) %>%
  filter(CommonName != "Unidentified rockfish") %>% # Removing unidentified rockfish as a spp.
  filter(CommonName != "Blackeye goby") %>% # Removing a benthic associated spp
  filter(CommonName != "Longfin sculpin") %>% # Removing a benthic associated spp
  filter(CommonName != "Red Irish lord") %>% # Removing a benthic associated spp
  filter(CommonName != "Painted greenling") %>% # Removing a benthic associated spp
  filter(CommonName != "Kelp greenling") %>% # Removing a benthic associated spp
  filter(CommonName != "Whitespotted greenling") %>% # Removing a benthic associated spp
  filter(CommonName != "Lingcod") %>% # Removing a benthic associated spp
  droplevels() %>%
  dplyr::select(-Method) 

## For just spp. complexes of pelagic fish
taxap_grp <- taxa %>%
  filter(Method == 1) %>%
  filter(CommonName != "Unidentified rockfish") %>% # Removing unidentified rockfish as a spp.
  mutate(CommonName = fct_collapse(CommonName,
                                   Rockfish = c("Black rockfish", "Copper rockfish", "Yellowtail rockfish"),
                                   Surfperch = c("Kelp perch", "Shiner perch", "Striped seaperch", "Pile perch"),

                                   Schooling = c("Tube-snout", "Pacific Herring"))) %>%
  group_by(SiteName, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb)) %>%
  as.data.frame() 

## For just RLS Method 2 species (benthic fish (demersal) & inverts)
taxab <- taxa %>%
  filter(Method == 2) %>%
  mutate(CommonName = fct_recode(CommonName, "Green sea urchin" = "Northern sea urchin"),
         CommonName = fct_recode(CommonName, "Blood star" = "Unidentified blood star"),
         CommonName = fct_recode(CommonName, "Chiton" = "Unidentified chiton"),
         CommonName = fct_recode(CommonName, "Hermit crab" = "Unidentified hermit crab")) %>%
  droplevels() %>%
  dplyr::select(-Method)

## For just spp. of all inverts (benthic - demersal fish)
M2benthicfish <- c("Blackeye goby", "Crescent gunnel", "Kelp greenling", "Longfin sculpin", "Penpoint gunnel", "Red Irish lord", "Scalyhead sculpin", "Smoothhead sculpin", "Whitespotted greenling", "Lingcod")

taxab_inv <- taxab %>%
  filter(!CommonName %in% M2benthicfish) %>%
  droplevels() 

## For just spp. of all benthic fish
taxap_benthic <- taxa %>%
  filter(Method == 1) %>%
  filter(CommonName %in% M2benthicfish) %>%
  dplyr::select(-Method) %>%
  ungroup()

taxab_fish <- taxab %>%
  rbind(taxap_benthic) %>% # Adding on the M1 observations of benthic fish
  filter(CommonName %in% M2benthicfish) %>%
  droplevels()

## For just spp. of all fish (pelagic + demersal fish)
taxa_fish <- taxap %>%
  rbind(taxab_fish) %>%
  group_by(SiteName, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb))


#### Community data: Relative abundance tables & transformations ----

## ALL SPP
allspp_wide <- taxa_all %>% # Spreading the data to wide format
  spread(key = CommonName, value = TaxaAb) %>%
  as.data.frame()
rownames(allspp_wide) <- allspp_wide$SiteName # Setting col as rownames
allspp_wide <- allspp_wide[,-1] # Removing the col used above
allspp_wide[is.na(allspp_wide)] <- 0 # (no NAs)
allspp_mat <- as.matrix(allspp_wide) # As matrix

allspp_rel <- make_relative(allspp_mat) # Relative abundance matrix of community data
allspp_hel <- decostand(allspp_mat, method = "hellinger") # Hellinger transformation of square root

## PELAGIC FISH SPP
pfish_wide <- taxap %>% # Spreading the data to wide format
  spread(key = CommonName, value = TaxaAb) %>%
  as.data.frame()
rownames(pfish_wide) <- pfish_wide$SiteName # Setting col as rownames
pfish_wide <- pfish_wide[,-1] # Removing the col used above
pfish_wide[is.na(pfish_wide)] <- 0 # (no NAs)
pfish_mat <- as.matrix(pfish_wide) # As matrix

pfish_rel <- make_relative(pfish_mat) # Relative abundance matrix of community data
pfish_hel <- decostand(pfish_mat, method = "hellinger") # Hellinger transformation of square root

## BENTHIC FISH SPP
bfish_wide <- taxab_fish %>% # Spreading the data to wide format
  spread(key = CommonName, value = TaxaAb) %>%
  as.data.frame()
rownames(bfish_wide) <- bfish_wide$SiteName # Setting col as rownames
bfish_wide <- bfish_wide[,-1] # Removing the col used above
bfish_wide[is.na(bfish_wide)] <- 0 # (no NAs)
bfish_mat <- as.matrix(bfish_wide) # As matrix

bfish_rel <- make_relative(bfish_mat) # Relative abundance matrix of community data
bfish_hel <- decostand(bfish_mat, method = "hellinger") # Hellinger transformation of square root

## ALL FISH SPP
fish_wide <- taxa_fish %>% # Spreading the data to wide format
  spread(key = CommonName, value = TaxaAb) %>%
  as.data.frame()
rownames(fish_wide) <- fish_wide$SiteName # Setting col as rownames
fish_wide <- fish_wide[,-1] # Removing the col used above
fish_wide[is.na(fish_wide)] <- 0 # (no NAs)
fish_mat <- as.matrix(fish_wide) # As matrix

fish_rel <- make_relative(fish_mat) # Relative abundance matrix of community data
fish_hel <- decostand(fish_mat, method = "hellinger") # Hellinger transformation of square root

## BENTHIC INVERT SPP
binv_wide <- taxab_inv %>% # Spreading the data to wide format
  spread(key = CommonName, value = TaxaAb) %>%
  as.data.frame()
rownames(binv_wide) <- binv_wide$SiteName # Setting col as rownames
binv_wide <- binv_wide[,-1] # Removing the col used above
binv_wide[is.na(binv_wide)] <- 0 # (no NAs)
binv_mat <- as.matrix(binv_wide) # As matrix

binv_rel <- make_relative(binv_mat) # Relative abundance matrix of community data
binv_hel <- decostand(binv_mat, method = "hellinger") # Hellinger transformation of square root


#### Community data: Variable selection? ----

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





#### All species: Generating models, stats ----

# Proportion of constrained variation: The amount of variation in spp abundances explained across sites by the included predictor variables 
# Adjusted R Squared: Amount of variation in the relative abundances of spp explained by the predictor variables when corrected for the number of predictor variables 
# For anova.cca() by term: https://stats.stackexchange.com/questions/405024/what-null-hypothesis-does-vegananova-ccaby-terms-test-against

# Running PERMANOVA to test for initial grouping factors
allspp_perm <- adonis2(allspp_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom significant (0.003)

# Generating the restricted model
allspp_1 <- capscale(allspp_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray") 
RsquareAdj(allspp_1)
## The structural and environmental variables explain 46.4% of variation (constrained axes)
# Adjusted R squared of this model is: 21%

# Generating the full model
allspp_2 <- capscale(allspp_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=binv_wide, add=FALSE, distance="bray")
RsquareAdj(allspp_2)
## Only the environmental variables explain explain 57.9% of variation (constrained axes)
# Adjusted R squared of this model is: 21.9%

anova.cca(allspp_2) # significance of model
anova.cca(allspp_2, step = 1000, by = "term") # significance of predictor terms
anova.cca(allspp_1, step = 1000, by = "axis") # significance by model axes

# Running PERMANOVA for effects of full model predictors
allspp_perm <- adonis2(allspp_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, permutations = 999, method="bray")

# Running SIMPER to test between groups
allspp_simp <- with(preds_ord, simper(allspp_hel, group=Kelpdom, permutations=999))

#### Pelagic fish: Generating models, stats ----

# Running PERMANOVA to test for initial grouping factors
pfish_perm <- adonis2(pfish_hel ~ Kelpdom+Tempave+Tempmax+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, permutations = 999, method="bray")

# Generating the restricted model
pfish_1 <- capscale(pfish_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray") 
RsquareAdj(pfish_1)
## The structural and environmental variables explain 34.6% of variation (constrained axes)
# Adjusted R squared of this model is: 8.1%

# Generating the full model
pfish_2 <- capscale(pfish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=pfish_wide, add=FALSE, distance="bray")
RsquareAdj(pfish_2)
## Only the environmental variables explain explain 50.1% of variation (constrained axes)
# Adjusted R squared of this model is: 19.3%

anova.cca(pfish_2) # significance of model
anova.cca(pfish_2, step = 1000, by = "term") # significance of predictor terms
anova.cca(pfish_2, step = 1000, by = "axis") # significance by model axes

# Running PERMANOVA for effects of full model predictors
pfish_perm <- adonis2(pfish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, permutations = 999, method="bray")

# Running SIMPER to test between groups
pfish_simp <- with(preds_ord, simper(pfish_hel, group=Kelpdom, permutations=999))
pfish_simpSPP <- c("Striped seaperch", "Shiner perch", "Black rockfish", "Pacific herring", "Tube-snout", "Copper rockfish")

#### Benthic fish: Generating models, stats ----

# Running PERMANOVA to test for initial grouping factors
bfish_perm <- adonis2(bfish_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray")

# Generating the restricted model
bfish_1 <- capscale(bfish_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray") 
RsquareAdj(bfish_1)
## The structural and environmental variables explain 48.7% of variation (constrained axes)
# Adjusted R squared of this model is: 32.4%

# Generating the full model
bfish_2 <- capscale(bfish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=bfish_wide, add=FALSE, distance="bray")
RsquareAdj(bfish_2)
## Only the environmental variables explain explain 59.5% of variation (constrained axes)
# Adjusted R squared of this model is: 37.3%

anova.cca(bfish_2) # significance of model
anova.cca(bfish_2, step = 1000, by = "term") # significance of predictor terms
anova.cca(bfish_2, step = 1000, by = "axis") # significance by model axes

# Running PERMANOVA for effects of full model predictors
bfish_perm <- adonis2(bfish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, permutations = 999, method="bray")

# Running SIMPER to test between groups
bfish_simp <- with(preds_ord, simper(bfish_hel, group=Kelpdom, permutations=999))
bfish_simpSPP <- c("Blackeye goby", "Scalyhead sculpin", "Longfin sculpin", "Smoothhead sculpin")


#### All fish: Generating models, stats ----

# Running PERMANOVA for initial grouping factors
fish_perm <- adonis2(fish_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom significant (0.015)

# Generating the restricted model
fish_1 <- capscale(fish_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=fish_wide, add=FALSE, distance="bray") 
RsquareAdj(fish_1)
## The structural and environmental variables explain 36.9% of variation (constrained axes)
# Adjusted R squared of this model is: 10.5

# Generating the full model
fish_2 <- capscale(fish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=fish_wide, add=FALSE, distance="bray")
RsquareAdj(fish_2)
## Only the environmental variables explain explain 50.1% of variation (constrained axes)
# Adjusted R squared of this model is: 13.6%

anova.cca(fish_2) # significance of model
anova.cca(fish_2, step = 1000, by = "term") # significance of predictor terms
anova.cca(fish_2, step = 1000, by = "axis") # significance by model axes

# Running PERMANOVA for effects of full model predictors
fish_perm <- adonis2(fish_hel ~ Kelpdom+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, permutations = 999, method="bray")

# Running SIMPER to test between groups
fish_simp <- with(preds_ord, simper(fish_hel, group=Kelpdom, permutations=999))
fish_simpSPP <- c("Blackeye goby", "Striped seaperch", "Shiner perch", "Scalyhead sculpin", "Pacific Herring", "Black rockfish", "Tube-snout", "Copper rockfish")

  
#### All inverts: Generating models, stats, and plots ----

# Running PERMANOVA to test for grouping factors
binv_perm <- adonis2(binv_hel ~ Kelpdom+Cluster+Kelpdom:Cluster, data=preds_ord, permutations = 999, method="bray") # Kelpdom & Cluster significant (0.012)
binv_perm_pair <- pairwise.adonis2(binv_hel ~ Cluster, data=preds_ord) # C3 vs C1 (0.004)

# Generating the restricted model
binv_1 <- capscale(binv_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray") 
RsquareAdj(binv_1)
## The structural and environmental variables explain 53.7% of variation (constrained axes)
# Adjusted R squared of this model is: 30.6%

# Generating the full model
binv_2 <- capscale(binv_hel ~ Kelpdom+Cluster+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray")
RsquareAdj(binv_2)
## Only the environmental variables explain explain 65% of variation (constrained axes)
# Adjusted R squared of this model is: 33.2%

anova.cca(binv_2) # significance of model
anova.cca(binv_2, step = 1000, by = "term") # significance of predictor terms
anova.cca(binv_2, step = 1000, by = "axis") # significance by model axes

# Running PERMANOVA for effects of full model predictors
binv_perm <- adonis2(binv_hel ~ Kelpdom+Cluster+DensityM+Area_m2+HeightM+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, permutations = 999, method="bray")

# Running SIMPER to test between groups (limiting to 10 most influential)
binv_simp <- with(preds_ord, simper(binv_hel, group=Kelpdom, permutations=999))
binv_simpSPP <- c("Red turban snail", "Purple sea urchin", "Red sea urchin", "Bat star", "Purple sea star", "Whitecap limpet", "Pinto abalone", "California sea cucumber", "Rainbow star", "Leather star")



### Plotting the visuals!

# Quick ordiplot
binv_ordplot <- ordiplot(binv_2, type="text", display="all")
ordisymbol(binv_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(binv_2, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)

## Prepping for the ggplot biplots

# Extracting the axes data from the model
binv_axis.long <- axis.long(binv_2, choices=c(1,2))
binv_axis.long
# Extracting the locations of species from the model
binv_species.long <- species.long(binv_ordplot)
head(binv_species.long)
# Extracting the locations of sites from the model
binv_sites.long <- sites.long(binv_ordplot, env.data=preds_ord)
head(allspp_sites.long)
# Extracting the locations of centroids from the sites.long output
binv_centroids.long <- centroids.long(binv_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(binv_centroids.long)
# Extracting env biplot correlations from the model
binv_envfit <- envfit(ord=binv_ordplot, env=preds_ord)
# Defining env vectors of interest
binv_envvectors <- c("DensityM", "Area_m2", "HeightM", "Tempave", "Depth_datum_m", "exp_36", "Punderstory", "Phardbottom", "Psoftbottom")
binv_vectorfit.long <- vectorfit.long(binv_envfit) %>%
  filter(vector %in% binv_envvectors)
head(binv_vectorfit.long)

# Extracting the most influential spp from the model and ordiplot
binv_species.envfit <- envfit(binv_ordsimple, env=binv_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
binv_species.envfit.table <- data.frame(r=binv_species.envfit$vectors$r, p=binv_species.envfit$vectors$pvals)
binv_species.long.var <- species.long(binv_ordsimple, spec.data=binv_species.envfit.table)
binv_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
binv_species.long.vargrt <- binv_species.long.var[binv_species.long.var$r >= 0.5, ]
binv_species.long.vargrt

# Selecting for SIMPER spp that explain > 70% of the observed variation among groups
binv_species_SIMPER <- binv_species.long.var %>%
  filter(labels %in% binv_simpSPP)
binv_species_SIMPER


## The biplot for environmental

# Radial shift function for arrow text
rshift = function(r, theta, a=0.04, b=0.08) {
  r + a + b*abs(cos(theta))
}

# Calculate shift of text from arrows
env.arrows = binv_vectorfit.long %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))


# The full environmental plot 
binv_ordplot_env <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(binv_axis.long[1, "label"]) +
  ylab(binv_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3.2,3)) + 
  scale_x_continuous(limits=c(-3.2,3)) + 
  geom_point(data=binv_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=4.5) +
  geom_segment(data=env.arrows, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="grey20", linewidth=0.5, arrow=arrow(length = unit(0.2, "cm"), type="closed")) +
  geom_text(data=env.arrows, # Adding in the biplot correlation labels
                  aes(x=xnew*2.25, y=ynew*2, label=vector),
                  colour="black", size=3.5) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values=c(21,22,23,24)) +
  theme_classic() +
  theme(legend.position=c(0.1,0.82),
        legend.text=element_text(size=11, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=11, color="black"),
        axis.title=element_text(size=12, color="black"))
binv_ordplot_env


## The biplot for influential spp

# Calculate shift of text from arrows
spp.arrows = binv_species_SIMPER %>% 
  mutate(r = sqrt(axis1^2 + axis2^2),
         theta = atan2(axis2,axis1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

# The full species plot
binv_ordplot_spp <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(binv_axis.long[1, "label"]) +
  ylab(binv_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3.2,3)) + 
  scale_x_continuous(limits=c(-3.2,3)) + 
  geom_point(data=binv_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, shape=Cluster, fill=Kelpdom),
             size=4.5) +
  geom_segment(data=spp.arrows, # Adding in the SIMPER arrows
               aes(x=0, y=0, xend=axis1*3, yend=axis2*3),
               colour="grey20", linewidth=0.5, arrow = arrow(length = unit(0.2, "cm"), type="closed")) +
  geom_text(data=spp.arrows, # Adding in the SIMPER species labels
                  aes(x=xnew*2.9, y=ynew*3, label=labels),
                  colour="black", fontface="italic", size=3) +
  scale_fill_manual(values=met.brewer("Kandinsky"), name=NULL, labels=c(expression(italic("N. luetkeana")), expression(italic("M. pyrifera")))) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values=c(21,22,23,24)) +
  theme_classic() +
  theme(legend.position="none", # Hiding legend for this plot
        legend.text=element_text(size=11, color="black"),
        legend.text.align = 0,
        legend.spacing.y = unit(0, "lines"),
        axis.text=element_text(size=11, color="black"),
        axis.title=element_text(size=12, color="black"))
binv_ordplot_spp

binv_ordplot_spp_draw <- binv_ordplot_spp + draw_image(
  "./Phylopic/Inverts/urchin.png",
  x = -3, y = 1.2, width = 0.6, height = 0.6) + draw_image(
    "./Phylopic/Inverts/urchin.png",
    x = -3.2, y = -1.8, width = 0.6, height = 0.6) + draw_image(
    "./Phylopic/Inverts/seastar1.png",
    x = -2.1, y = -0.3, width = 0.6, height = 0.6) + draw_image(
      "./Phylopic/Inverts/seastar2.png",
      x = 1.1, y = 0.7, width = 0.6, height = 0.6) + draw_image(
        "./Phylopic/Inverts/seastar3.png",
        x = 1, y = -2, width = 0.6, height = 0.6) + draw_image(
          "./Phylopic/Inverts/seacucumber.png",
          x = 0.9, y = -1, width = 0.5, height = 0.5) + draw_image(
            "./Phylopic/Inverts/marinesnail.png",
            x = 0.2, y = -3, width = 0.5, height = 0.5) + draw_image(
              "./Phylopic/Inverts/abalone.png",
              x = -0.3, y = -1.4, width = 0.4, height = 0.4) + draw_image(
                "./Phylopic/Inverts/limpet.png",
                x = -1.1, y = -0.8, width = 0.4, height = 0.4)
  


## Putting the two together

tiff(file="./MSc_plots/PaperFigs/Invert_ord_draft.tiff", height = 11.5, width = 7.5, units = "in", res=400)

binv_ordplots <- ggarrange(binv_ordplot_env, binv_ordplot_spp_draw, ncol=1, common.legend=FALSE)
binv_ordplots

dev.off()


#### Community data ordination: Plotting in ordination space ----

# For the ggplots (https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html)

colsmap <- c("#225555", "#4477aa", "#e4632d", "#997700")
pchs <- c(21,22,23,24)



### ALL SPP

## OPTION 1: HELLINGER
allspp_capmodel <- dbrda(allspp_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray") 
## The structural and environmental variables explain 58% of variation (constrained axes)
# Adjusted R squared of this model is: 19.8%
allspp_capmodel <- dbrda(allspp_hel ~ Cluster+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=allspp_wide, add=FALSE, distance="bray") 
## Only the environmental variables explain explain 55.8% of variation (constrained axes)
# Adjusted R squared of this model is: 22.7%

anova.cca(allspp_capmodel) # significance of model
anova.cca(allspp_capmodel, step = 1000, by = "term") # significance of predictor terms
anova.cca(allspp_capmodel, step = 1000, by = "axis") # significance by model axes



## Quick ordiplot to visualize
allspp_ordplot <- ordiplot(allspp_capmodel, type="text", display="all")
ordisymbol(allspp_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(allspp_capmodel, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)

## GGplot ordiplot to finalize

# Extracting the axes data from the model
allspp_axis.long <- axis.long(allspp_capmodel, choices=c(1,2))
allspp_axis.long
# Extracting the locations of species from the model
allspp_species.long <- species.long(allspp_ordplot)
head(allspp_species.long)
# Extracting the locations of sites from the model
allspp_sites.long <- sites.long(allspp_ordplot, env.data=preds_ord)
head(allspp_sites.long)
# Extracting the locations of centroids from the sites.long output
allspp_centroids.long <- centroids.long(allspp_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(allspp_centroids.long)
# Extracting env biplot correlations from the model
allspp_envfit <- envfit(ord=allspp_ordplot, env=preds_ord)
allspp_envvectors <- c("DensityM", "Area_m2", "HeightM", "Tempave", "Depth_datum_m", "exp_36", "Punderstory", "Phardbottom", "Psoftbottom") # Defining env vectors of interest
allspp_vectorfit.long <- vectorfit.long(allspp_envfit) %>%
  filter(vector %in% allspp_envvectors)
head(allspp_vectorfit.long)

# Extracting the ellipses from the model and ordiplot
allspp_ordsimple <- ordiplot(allspp_capmodel)
allspp_ordellipses <- with(preds_ord, ordiellipse(allspp_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
allspp_ellipses <- ordiellipse.long(allspp_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
allspp_species.envfit <- envfit(allspp_ordsimple, env=allspp_hel) # For my Bray distance ordiplot of Hellinger transformed spp data
allspp_species.envfit.table <- data.frame(r=allspp_species.envfit$vectors$r, p=allspp_species.envfit$vectors$pvals)
allspp_species.long.var <- species.long(allspp_ordsimple, spec.data=allspp_species.envfit.table)
allspp_species.long.var
# Selecting for individual spp that explain > 50% of the observed variation among groups
allspp_species.long.vargrt <- allspp_species.long.var[allspp_species.long.var$r >= 0.5, ]
allspp_species.long.vargrt

# Selecting for SIMPER spp that explain > 70% of the observed variation among groups
allspp_species_SIMPER <- allspp_species.long.var %>%
  filter(labels %in% allspp_simpSPP)
allspp_species_SIMPER


# The ggplot
allspp_ordplot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(allspp_axis.long[1, "label"]) +
  ylab(allspp_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  # geom_mark_ellipse(data=allspp_sites.long, # Ellipses hulls for all group points
  #                   aes(x=axis1, y=axis2, colour=Kelpdom,
  #                       fill=after_scale(alpha(colour, 0.2))),
  #                   expand=0, linewidth=0.2, show.legend=FALSE) +
  # geom_polygon(data=allspp_ellipses,
  #              aes(x=axis1, y=axis2, colour=Kelpdom,
  #                  fill=after_scale(alpha(colour, 0.2))),
  #              size=0.2, show.legend=FALSE) +
  geom_point(data=allspp_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom),
             size=4, shape=21) +
  # geom_point(data=allspp_centroids.long, # Adding in the centroids for groups
  #            aes(x=axis1c, y=axis2c, color=Centroid), size=4, shape=4, show.legend=F) +
  geom_segment(data=allspp_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1, yend=axis2),
               colour="black", size=0.7, arrow=arrow()) +
  geom_text_repel(data=allspp_vectorfit.long, # Adding in the biplot correlation labels
                  aes(x=axis1, y=axis2, label=vector),
                  colour="black") +
  # geom_segment(data=allspp_species.long.vargrt, # Adding in the influential species data
  #              aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
  #              colour="black", linewidth=0.7, arrow=arrow()) +
  geom_text_repel(data=allspp_species.long.vargrt, # Adding in the species labels
                  aes(x=axis1*4, y=axis2*4, label=labels),
                  colour="black", fontface="italic") +
  scale_fill_manual(values=colsmap) +
  theme_minimal()
allspp_ordplot



### PELAGIC FISH SPP

## OPTION 1: HELLINGER
pfish_capmodel <- capscale(pfish_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=pfish_wide, add=FALSE, dist="bray") 

## Quick ordiplot to visualize
pfish_ordplot <- ordiplot(pfish_capmodel, type="text", display="all")
ordisymbol(pfish_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(pfish_capmodel, Kelpdom, col=colsmap, kind = "se", conf=0.95, label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)

## GGplot ordiplot to finalize

# Extracting the axes data from the model
pfish_axis.long <- axis.long(pfish_capmodel, choices=c(1,2))
pfish_axis.long
# Extracting the locations of species from the model
pfish_species.long <- species.long(pfish_ordplot)
head(pfish_species.long)
# Extracting the locations of sites from the model
pfish_sites.long <- sites.long(pfish_ordplot, env.data=preds_ord)
head(pfish_sites.long)
# Extracting the locations of centroids from the sites.long output
pfish_centroids.long <- centroids.long(pfish_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(pfish_centroids.long)
# Extracting env biplot correlations from the model
pfish_envfit <- envfit(ord=pfish_ordplot, env=preds_ord)
pfish_envvectors <- c("Tempave", "Depth_datum_m", "exp_36") # Defining env vectors of interest
pfish_vectorfit.long <- vectorfit.long(pfish_envfit) %>%
  filter(vector %in% pfish_envvectors)
head(pfish_vectorfit.long)

# Extracting the ellipses from the model and ordiplot
pfish_ordsimple <- ordiplot(pfish_capmodel)
pfish_ordellipses <- with(preds_ord, ordiellipse(pfish_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
pfish_ellipses <- ordiellipse.long(pfish_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
pfish_species.envfit <- envfit(pfish_ordsimple, env=pfish_hel) # The simple ordiplot above & Hellinger transformed spp data
pfish_species.envfit.table <- data.frame(r=pfish_species.envfit$vectors$r, p=pfish_species.envfit$vectors$pvals)
pfish_species.long.var <- species.long(pfish_ordsimple, spec.data=pfish_species.envfit.table)
pfish_species.long.var
# Selecting for spp that explain > 60% of the observed variation among groups
pfish_species.long.vargrt <- pfish_species.long.var[pfish_species.long.var$r >= 0.6, ]
pfish_species.long.vargrt

# Selecting for SIMPER spp that explain > 70% of the observed variation among groups
pfish_species_SIMPER <- pfish_species.long.var %>%
  filter(labels %in% pfish_simpSPP)
pfish_species_SIMPER


# The ggplot
pfish_ordplot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(pfish_axis.long[1, "label"]) +
  ylab(pfish_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  geom_mark_ellipse(data=pfish_sites.long, # Ellipses hulls for all group points
                    aes(x=axis1, y=axis2, colour=Kelpdom,
                        fill=after_scale(alpha(colour, 0.2))),
                    expand=0, linewidth=0.2, show.legend=FALSE) +
  # geom_polygon(data=pfish_ellipses, # Ellipses hulls for CIs
  #              aes(x=axis1, y=axis2, colour=Kelpdom,
  #                  fill=after_scale(alpha(colour, 0.2))),
  #              size=0.2, show.legend=FALSE) +
  geom_point(data=pfish_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, colour=Kelpdom),
             size=4) +
  geom_point(data=pfish_centroids.long, # Adding in the centroids for groups
             aes(x=axis1c, y=axis2c, color=Centroid), size=4, shape=4, show.legend=F) +
  geom_segment(data=pfish_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
               colour="red", linewidth=0.7, arrow=arrow()) +
  geom_text_repel(data=pfish_vectorfit.long, # Adding in the biplot correlation labels
                  aes(x=axis1, y=axis2, label=vector),
                  colour="red", size=5) +
  # geom_segment(data=pfish_species_SIMPER, # Adding in the SIMPER arrows
  #              aes(x=0, y=0, xend=axis1*2, yend=axis2*2), 
  #              colour="black", linewidth=0.7, arrow=arrow()) +
  geom_text_repel(data=pfish_species_SIMPER, # Adding in the SIMPER species labels
                  aes(x=axis1*2, y=axis2*2, label=labels),
                  colour="black", fontface="italic", size=4) +
  scale_color_manual(values=colsmap) +
  scale_shape_manual(values=pchs) +
  theme_minimal()
pfish_ordplot



### BENTHIC FISH SPP

## OPTION 1: HELLINGER
bfish_capmodel <- rda(bfish_hel ~ Kelpdom+Tempave+Tempmin+exp_36+Depth_datum_m+Punderstory+Psoftbottom+Phardbottom, data=preds_ord, comm=bfish_wide, add=FALSE, dist="bray") 

## Quick ordiplot to visualize
bfish_ordplot <- ordiplot(bfish_capmodel, type="text", display="all")
ordisymbol(bfish_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(bfish_capmodel, Kelpdom, col=colsmap, kind = "se", conf=0.95, label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)

## GGplot ordiplot to finalize

# Extracting the axes data from the model
bfish_axis.long <- axis.long(bfish_capmodel, choices=c(1,2))
bfish_axis.long
# Extracting the locations of species from the model
bfish_species.long <- species.long(bfish_ordplot)
head(bfish_species.long)
# Extracting the locations of sites from the model
bfish_sites.long <- sites.long(bfish_ordplot, env.data=preds_ord)
head(bfish_sites.long)
# Extracting the locations of centroids from the sites.long output
bfish_centroids.long <- centroids.long(bfish_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(bfish_centroids.long)
# Extracting env biplot correlations from the model
bfish_envfit <- envfit(ord=bfish_ordplot, env=preds_ord)
bfish_envvectors <- c("Tempave", "Tempmin", "Depth_datum_m", "exp_36", "Punderstory", "Phardbottom", "Psoftbottom") # Defining env vectors of interest
bfish_vectorfit.long <- vectorfit.long(bfish_envfit) %>%
  filter(vector %in% bfish_envvectors)
head(bfish_vectorfit.long)

# Extracting the ellipses from the model and ordiplot
bfish_ordsimple <- ordiplot(bfish_capmodel)
bfish_ordellipses <- with(preds_ord, ordiellipse(bfish_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
bfish_ellipses <- ordiellipse.long(bfish_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
bfish_species.envfit <- envfit(bfish_ordsimple, env=bfish_hel) # The simple ordiplot above & Hellinger transformed spp data
bfish_species.envfit.table <- data.frame(r=bfish_species.envfit$vectors$r, p=bfish_species.envfit$vectors$pvals)
bfish_species.long.var <- species.long(bfish_ordsimple, spec.data=bfish_species.envfit.table)
bfish_species.long.var
# Selecting for spp that explain > 60% of the observed variation among groups
bfish_species.long.vargrt <- bfish_species.long.var[bfish_species.long.var$r >= 0.5, ]
bfish_species.long.vargrt

# Selecting for SIMPER spp that explain > 70% of the observed variation among groups
bfish_species_SIMPER <- bfish_species.long.var %>%
  filter(labels %in% bfish_simpSPP)
bfish_species_SIMPER


# The ggplot
bfish_ordplot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(bfish_axis.long[1, "label"]) +
  ylab(bfish_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  # geom_mark_ellipse(data=bfish_sites.long, # Ellipses hulls for all group points
  #                   aes(x=axis1, y=axis2, colour=Kelpdom,
  #                       fill=after_scale(alpha(colour, 0.2))),
                    # expand=0, linewidth=0.2, show.legend=FALSE) +
  # geom_polygon(data=bfish_ellipses, # Ellipses hulls for CIs
  #              aes(x=axis1, y=axis2, colour=Kelpdom,
  #                  fill=after_scale(alpha(colour, 0.2))),
  #              size=0.2, show.legend=FALSE) +
  geom_point(data=bfish_sites.long, size=4.5, shape=21, # Adding in the main site data points
             aes(x=axis1, y=axis2, fill=Kelpdom)) +
  # geom_point(data=bfish_centroids.long, # Adding in the centroids for groups
  #            aes(x=axis1c, y=axis2c, color=Centroid), size=4, shape=4, show.legend=F) +
  geom_segment(data=bfish_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1, yend=axis2),
               colour="black", size=0.7, arrow=arrow()) +
  geom_text_repel(data=bfish_vectorfit.long, # Adding in the biplot correlation labels
                  aes(x=axis1, y=axis2, label=vector),
                  colour="black", size=5) +
  # geom_segment(data=bfish_species_SIMPER, # Adding in the SIMPER arrows
  #              aes(x=0, y=0, xend=axis1*2, yend=axis2*2), 
  #              colour="black", linewidth=0.7, arrow=arrow()) +
  geom_text_repel(data=bfish_species_SIMPER, # Adding in the SIMPER species labels
                  aes(x=axis1*2, y=axis2*2, label=labels),
                  colour="black", fontface="italic", size=4) +
  scale_fill_manual(values=colsmap) +
  theme_classic()
bfish_ordplot



### BENTHIC INVERT SPP

## OPTION 1: HELLINGER
binv_capmodel <- dbrda(binv_hel ~ Kelpdom+Tempave+exp_36+Depth_datum_m, data=preds_ord, comm=binv_wide, add=FALSE) 

## Quick ordiplot to visualize
binv_ordplot <- ordiplot(binv_capmodel, type="text", display="all")
ordisymbol(binv_ordplot, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(binv_capmodel, Kelpdom, col=colsmap, kind = "se", conf=0.95, label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)

## GGplot ordiplot to finalize

# Extracting the axes data from the model
binv_axis.long <- axis.long(binv_capmodel, choices=c(1,2))
binv_axis.long
# Extracting the locations of species from the model
binv_species.long <- species.long(binv_ordplot)
head(binv_species.long)
# Extracting the locations of sites from the model
binv_sites.long <- sites.long(binv_ordplot, env.data=preds_ord)
head(binv_sites.long)
# Extracting the locations of centroids from the sites.long output
binv_centroids.long <- centroids.long(binv_sites.long, grouping=preds_ord$Kelpdom, centroids.only=FALSE)
head(binv_centroids.long)
# Extracting env biplot correlations from the model
binv_envfit <- envfit(ord=binv_ordplot, env=preds_ord)
binv_envvectors <- c("Tempave", "Depth_datum_m", "exp_36") # Defining env vectors of interest
binv_vectorfit.long <- vectorfit.long(binv_envfit) %>%
  filter(vector %in% binv_envvectors)
head(binv_vectorfit.long)

# Extracting the ellipses from the model and ordiplot
binv_ordsimple <- ordiplot(binv_capmodel)
binv_ordellipses <- with(preds_ord, ordiellipse(binv_ordsimple, groups=Kelpdom, kind = "se", conf=0.95))
binv_ellipses <- ordiellipse.long(binv_ordellipses, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
binv_species.envfit <- envfit(binv_ordsimple, env=binv_hel) # The simple ordiplot above & Hellinger transformed spp data
binv_species.envfit.table <- data.frame(r=binv_species.envfit$vectors$r, p=binv_species.envfit$vectors$pvals)
binv_species.long.var <- species.long(binv_ordsimple, spec.data=binv_species.envfit.table)
binv_species.long.var
# Selecting for spp that explain > 50% of the observed variation among groups
binv_species.long.vargrt <- binv_species.long.var[binv_species.long.var$r >= 0.5, ]
binv_species.long.vargrt

# Selecting for SIMPER spp that explain > 70% of the observed variation among groups
binv_species_SIMPER <- binv_species.long.var %>%
  filter(labels %in% binv_simpSPP)
binv_species_SIMPER


# The ggplot
binv_ordplot <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(binv_axis.long[1, "label"]) +
  ylab(binv_axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  geom_mark_ellipse(data=binv_sites.long, # Ellipses hulls for all group points
                    aes(x=axis1, y=axis2, colour=Kelpdom,
                        fill=after_scale(alpha(colour, 0.2))),
                    expand=0, linewidth=0.2, show.legend=FALSE) +
  # geom_polygon(data=binv_ellipses, # Ellipses hulls for CIs
  #              aes(x=axis1, y=axis2, colour=Kelpdom,
  #                  fill=after_scale(alpha(colour, 0.2))),
  #              size=0.2, show.legend=FALSE) +
  geom_point(data=binv_sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, colour=Kelpdom),
             size=4) +
  geom_point(data=binv_centroids.long, # Adding in the centroids for groups
             aes(x=axis1c, y=axis2c, color=Centroid), size=4, shape=4, show.legend=F) +
  geom_segment(data=binv_vectorfit.long, # Adding in the biplot correlation arrows
               aes(x=0, y=0, xend=axis1, yend=axis2),
               colour="red", size=0.7, arrow=arrow()) +
  geom_text_repel(data=binv_vectorfit.long, # Adding in the biplot correlation labels
                  aes(x=axis1, y=axis2, label=vector),
                  colour="red", size=5) +
  # geom_segment(data=binv_species_SIMPER, # Adding in the SIMPER arrows
  #              aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
  #              colour="black", linewidth=0.7, arrow=arrow()) +
  geom_text_repel(data=binv_species_SIMPER, # Adding in the SIMPER species labels
                  aes(x=axis1*4, y=axis2*4, label=labels),
                  colour="black", fontface="italic", size=4) +
  scale_color_manual(values=colsmap) +
  scale_shape_manual(values=pchs) +
  theme_minimal()
binv_ordplot












### BY CLUSTER GROUPS

## OPTION 1: HELLINGER- Balanced emphasis for spp of all abundances (down weighting of highly abundant spp, no emphasis on rare spp)
capmodel <- capscale(sqrt(taxamat) ~ Cluster+Kelpdom, data=preds_ord, comm=taxawide, add=FALSE, dist="hellinger") 
# ## OPTION 2: BRAY/WISCONSIN- Emphasis on the influence of most abundant spp (lower weighting for rare species)
# capmodel <- capscale(wisconsin(taxarel^(1/4)) ~ Cluster+Kelpdom, data=preds_ord, comm=taxawide, add=FALSE, dist="bray")

## Quick ordiplot to visualize
ordplot <- ordiplot(capmodel, type="text", display="all")
ordisymbol(ordplot, preds_ord, "Cluster", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(capmodel, Cluster, col=colsmap, kind = "se", conf=0.95, label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Cluster), bty="n", pch=pchs, col = colsmap)

## Clean ggplots 
# https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html

# Extracting the axes data from the model
axis.long <- axis.long(capmodel, choices=c(1,2))
axis.long
# Extracting the locations of species from the model
species.long <- species.long(ordplot)
head(species.long)
# Extracting the locations of sites from the model
sites.long <- sites.long(ordplot, env.data=preds_ord)
head(sites.long)
# Extracting the locations of centroids from the sites.long output
centroids.long <- centroids.long(sites.long, grouping=preds_ord$Cluster, centroids.only=FALSE)
head(centroids.long)
# # Extracting env biplot correlations from the model
# envfit <- envfit(ord=ordplot, env=preds_ord)
# envvectors <- c("Tempave", "Depth_datum_m", "exp_36") # Defining env vectors of interest
# vectorfit.long <- vectorfit.long(envfit) %>%
#   filter(vector %in% envvectors)
# head(vectorfit.long)

# Extracting the CI ellipses from the model and ordiplot
ordsimple <- ordiplot(capmodel)
ordellipses <- with(preds_ord, ordiellipse(ordsimple, groups=Cluster, kind = "se", conf=0.95))
ellipses <- ordiellipse.long(ordellipses, grouping.name="Cluster")

# Extracting the most influential spp from the model and ordiplot
species.envfit <- envfit(ordsimple, env=taxahel) # The simple ordiplot above & Hellinger transformed spp data
species.envfit.table <- data.frame(r=species.envfit$vectors$r, p=species.envfit$vectors$pvals)
species.long.var <- species.long(ordsimple, spec.data=species.envfit.table)
species.long.var
# Selecting for spp that explain > 60% of the observed variation among groups
species.long.vargrt <- species.long.var[species.long.var$r >= 0.6, ]
species.long.vargrt


# The ggplot
plotgg <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long[1, "label"]) +
  ylab(axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-2.3,2)) + 
  scale_x_continuous(limits=c(-2,2)) + 
  geom_mark_ellipse(data=sites.long, # Ellipses hulls for all group points
                    aes(x=axis1, y=axis2, colour=Cluster,
                        fill=after_scale(alpha(colour, 0.2))),
                    expand=0, linewidth=0.2, show.legend=FALSE) +
  # geom_polygon(data=ellipses, # 95% confidence ellipses
  #              aes(x=axis1, y=axis2, colour=Cluster,
  #                  fill=after_scale(alpha(colour, 0.2))),
  #              size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long, # Adding in the main site data points
             aes(x=axis1, y=axis2, colour=Cluster, shape=Kelpdom),
             size=4) +
  geom_point(data=centroids.long, # Adding in the centroids for groups
             aes(x=axis1c, y=axis2c, color=Centroid), size=4, shape=4, show.legend=F) +
  # geom_segment(data=vectorfit.long, # Adding in the biplot correlation arrows
  #              aes(x=0, y=0, xend=axis1, yend=axis2), 
  #              colour="red", size=0.7, arrow=arrow()) +
  # geom_text_repel(data=vectorfit.long, # Adding in the biplot correlation labels
  #                 aes(x=axis1, y=axis2, label=vector),
  #                 colour="red") +
  geom_segment(data=species.long.vargrt, # Adding in the influential species data
               aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
               colour="black", linewidth=0.7, arrow=arrow()) +
  geom_text_repel(data=species.long.vargrt, # Adding in the species labels
                  aes(x=axis1*4, y=axis2*4, label=labels),
                  colour="black") +
  scale_color_manual(values=colsmap) +
  scale_shape_manual(values=pchs) +
  theme_minimal()
plotgg


#////////////////////////////#


### BY KELP SPP GROUPS

## OPTION 1: HELLINGER- Balanced emphasis for spp of all abundances (down weighting of highly abundant spp, no emphasis on rare spp)
capmodel2 <- capscale(taxahel ~ Kelpdom, data=preds_ord, comm=taxawide, add=FALSE) 
## OPTION 2: BRAY/WISCONSIN- Emphasis on the influence of most abundant spp (lower weighting for rare species)
capmodel2 <- capscale(taxabray ~ Kelpdom, data=preds_ord, comm=taxawide, add=FALSE)


## Quick ordiplot to visualize
ordplot2 <- ordiplot(capmodel2, type="text", display="all")
ordisymbol(ordplot2, preds_ord, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_ord, ordiellipse(capmodel2, Kelpdom, col=colsmap, kind = "ehull", label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_ord$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Clean ggplots 
# https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html

# Extracting the axes data from the model
axis.long2 <- axis.long(capmodel2, choices=c(1,2))
axis.long2
# Extracting the locations of species from the model
species.long2 <- species.long(ordplot2)
head(species.long2)
# Extracting the locations of sites from the model
sites.long2 <- sites.long(ordplot2, env.data=preds_ord)
head(sites.long2)
# Extracting the locations of centroids from the sites.long output
centroids.long2 <- centroids.long(sites.long2, grouping=preds_ord$Kelpdom, centroids.only=TRUE)
head(centroids.long2)
# # Extracting env biplot correlations from the model
# envfit2 <- envfit(ord=ordplot2, env=preds_ord)
# envvectors <- c("Tempave", "Depth_datum_m", "exp_36") # Defining env vectors of interest
# vectorfit.long2 <- vectorfit.long(envfit2) %>%
#   filter(vector %in% envvectors)
# head(vectorfit.long2)

# Extracting the CI ellipses from the ordiplot above
ordsimple2 <- ordiplot(capmodel2)
ordellipses2 <- with(preds_ord, ordiellipse(ordsimple2, groups=Kelpdom, kind = "se", conf=0.95))
ellipses2 <- ordiellipse.long(ordellipses2, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
species.envfit2 <- envfit(ordsimple2, env=taxahel) # The simple ordiplot above & Hellinger transformed spp data (see previous section)
species.envfit.table2 <- data.frame(r=species.envfit2$vectors$r, p=species.envfit2$vectors$pvals)
species.long.var2 <- species.long(ordsimple2, spec.data=species.envfit.table2)
species.long.var2
# Selecting for spp that explain > 60% of the observed variance among groups
species.long.vargrt2 <- species.long.var2[species.long.var2$r >= 0.6, ]
species.long.vargrt2

# The ggplot
plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long[1, "label"]) +
  ylab(axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-3,3)) + 
  scale_x_continuous(limits=c(-3,3)) + 
  # geom_mark_ellipse(data=sites.long2, # Ellipses hulls for all group points
  #                   aes(x=axis1, y=axis2, colour=Kelpdom,
  #                       fill=after_scale(alpha(colour, 0.2))),
  #                   expand=0, linewidth=0.2, show.legend=FALSE) +
  geom_polygon(data=ellipses2, # 95% confidence ellipses
               aes(x=axis1, y=axis2, colour=Kelpdom,
                   fill=after_scale(alpha(colour, 0.2))),
               size=0.2, show.legend=FALSE) +
  geom_point(data=sites.long2, # The main site data points
             aes(x=axis1, y=axis2, colour=Kelpdom, shape=Cluster),
             size=4) +
  geom_point(data=centroids.long2, # Adding in the centroids for groups
             aes(x=axis1c, y=axis2c, color=Centroid), size=4, shape=4, show.legend=F) +
  # geom_segment(data=vectorfit.long2, # Adding in the biplot correlation arrows
  #              aes(x=0, y=0, xend=axis1, yend=axis2), 
  #              colour="red", size=0.7, arrow=arrow()) +
  # geom_text_repel(data=vectorfit.long2, # Adding in the biplot correlation labels
  #                 aes(x=axis1, y=axis2, label=vector),
  #                 colour="red") +
  geom_segment(data=species.long.vargrt2, # Adding in the influential species data
               aes(x=0, y=0, xend=axis1, yend=axis2), 
               colour="black", size=0.7, arrow=arrow()) +
  geom_text_repel(data=species.long.vargrt2, # Adding in the species labels
                  aes(x=axis1, y=axis2, label=labels),
                  colour="black") +
  scale_shape_manual(values=pchs) +
  scale_color_manual(values=rev(colsmap)) +
  theme_minimal()
plotgg2





