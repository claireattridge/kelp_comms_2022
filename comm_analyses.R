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

# Slope or aspect (TBA)
slope <- read_csv()

# Kelp metrics 
kelps <- read_csv("./MSc_data/Data_new/kelpmetrics_2022.csv") %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  dplyr::select(-c(HeightSD, BiomassSD, DensitySD, MacroSD, NereoSD))
kelps <- kelps %>%
  mutate(MacroP = MacroM/DensityM) # Creating a column for proportion Macro at a given site


# Listing the variables of interest together
dfs1 <- list(kelps, temps, rei)

# Joining the variables into one dataframe
preds <- reduce(dfs1, dplyr::left_join, by="SiteName")
preds[is.na(preds)] <- 0 # Translating the NAs to zeros for clustering analysis

preds <- preds[-17,] # Removing second beach south as an outlier site 

# Polishing dataframe with environmental data to use for later use in CAP ordination
preds_cap <- preds %>%
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

#### Cluster analysis ----

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


#### Community data cleaning ----

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


#### Community data ordination ----

library(funrar)
library(ecodist)
library(BiodiversityR)
library(ggforce)
library(ggrepel)

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
  droplevels() %>%
  dplyr::select(-Method)

## For just spp. complexes of pelagic fish
taxap_grp <- taxap %>%
  mutate(CommonName = fct_collapse(CommonName,
                                   Rockfish = c("Black rockfish", "Copper rockfish", "Yellowtail rockfish"),
                                   Surfperch = c("Kelp perch", "Shiner perch", "Striped seaperch", "Pile perch"),
                                   Greenlings = c("Kelp greenling", "Whitespotted greenling", "Lingcod"),
                                   Bottomdwelling = c("Blackeye goby", "Longfin sculpin", "Red Irish lord", "Painted greenling"),
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

M2benthicfish <- c("Blackeye goby", "Crescent gunnel", "Kelp greenling", "Longfin sculpin", "Penpoint gunnel", "Red Irish lord", "Scalyhead sculpin", "Smoothhead sculpin", "Whitespotted greenling")

taxab_inv <- taxab %>%
  filter(!CommonName %in% M2benthicfish) %>%
  droplevels() 

## For just spp. of all benthic fish (benthic - inverts)
taxab_fish <- taxab %>%
  filter(CommonName %in% M2benthicfish) %>%
  filter(!CommonName %in% "Whitespotted greenling") %>%
  filter(!CommonName %in% "Kelp greenling") %>%
  droplevels() 

## For just spp. complexes of all inverts (benthic - demersal fish)

## For just spp. of all fish (pelagic + demersal fish)
  
## For just spp. complexes of all fish (pelagic + demersal fish)
  
## For just spp. complexes of all spp (pelagic + benthic)


### Spreading the data to wide format
taxawide <- taxap_grp %>% # Swap for any if the above community dataframes in this line, then run the lines below-
  spread(key = CommonName, value = TaxaAb)


# Making table of relative abundances
rownames(taxawide) <- taxawide$SiteName # Setting col as rownames
taxawide <- taxawide[,-1] # Removing the col used above
taxawide[is.na(taxawide)] <- 0 # (no NAs)

taxamat <- as.matrix(taxawide) # As matrix
taxarel <- make_relative(taxamat) # Relative abundance matrix of community data


#### CAP ordination (simplified dbRDA) and plots

### BY CLUSTER GROUPS

## OPTION 1: HELLINGER- Balanced emphasis for spp of all abundances (down weighting of highly abundant spp, no emphasis on rare spp)
capmodel <- capscale(sqrt(taxamat) ~ Cluster+exp_36+Tempave+Depth_datum_m, data=preds_cap, comm=taxawide, add=FALSE, dist="hellinger") 
# ## OPTION 2: BRAY/WISCONSIN- Emphasis on the influence of most abundant spp (lower weighting for rare species)
# capmodel <- capscale(wisconsin(taxarel^(1/4)) ~ Cluster, data=preds_cap, comm=taxawide, add=FALSE, dist="bray")

colsmap <- c("#225555", "#4477aa", "#e4632d", "#997700")
pchs <- c(15:18)

## Quick ordiplot to visualize
ordplot <- ordiplot(capmodel, type="text", display="all")
ordisymbol(ordplot, preds_cap, "Cluster", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_cap, ordiellipse(capmodel, Cluster, col=colsmap2, kind = "se", conf=0.95, label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_cap$Cluster), bty="n", pch=pchs, col = colsmap)

## Clean ggplots 
# https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html

# Extracting the axes data from the model
axis.long <- axis.long(capmodel, choices=c(1,2))
axis.long
# Extracting the locations of species from the model
species.long <- species.long(ordplot)
head(species.long)
# Extracting the locations of sites from the model
sites.long <- sites.long(ordplot, env.data=preds_cap)
head(sites.long)
# Extracting the locations of centroids from the sites.long output
centroids.long <- centroids.long(sites.long, grouping=preds_cap$Cluster, centroids.only=FALSE)
head(centroids.long)
# Extracting env biplot correlations from the model
envfit <- envfit(ord=ordplot, env=preds_cap)
envvectors <- c("Tempave", "Depth_datum_m", "exp_36") # Defining env vectors of interest
vectorfit.long <- vectorfit.long(envfit) %>%
  filter(vector %in% envvectors)
head(vectorfit.long)

# Extracting the CI ellipses from the model and ordiplot
ordsimple <- ordiplot(capmodel)
ordellipses <- with(preds_cap, ordiellipse(ordsimple, groups=Cluster, kind = "se", conf=0.95))
ellipses <- ordiellipse.long(ordellipses, grouping.name="Cluster")

# Extracting the most influential spp from the model and ordiplot
species.envfit <- envfit(ordsimple, env=taxahel) # The simple ordiplot above & Hellinger transformed spp data (see previous section)
species.envfit.table <- data.frame(r=species.envfit$vectors$r, p=species.envfit$vectors$pvals)
species.long.var <- species.long(ordsimple, spec.data=species.envfit.table)
species.long.var
# Selecting for spp that explain > 60% of the observed variance among groups
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
  # geom_mark_ellipse(data=sites.long, # Ellipses hulls for all group points
  #                   aes(x=axis1, y=axis2, colour=Cluster,
  #                       fill=after_scale(alpha(colour, 0.2))),
  #                   expand=0, linewidth=0.2, show.legend=FALSE) +
  geom_polygon(data=ellipses, # 95% confidence ellipses
               aes(x=axis1, y=axis2, colour=Cluster,
                   fill=after_scale(alpha(colour, 0.2))),
               size=0.2, show.legend=FALSE) +
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
               colour="black", size=0.7, arrow=arrow()) +
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
capmodel2 <- capscale(sqrt(taxamat) ~ Kelpdom+exp_36+Tempave+Depth_datum_m, data=preds_cap, comm=taxawide, add=FALSE, dist="hellinger") 
# ## OPTION 2: BRAY/WISCONSIN- Emphasis on the influence of most abundant spp (lower weighting for rare species)
# capmodel2 <- capscale(wisconsin(taxarel^(1/4)) ~ Kelpdom, data=preds_cap, comm=taxawide, add=FALSE, dist="bray")

colsmap <- c("#225555", "#4477aa", "#e4632d", "#997700")
pchs <- c(15:18)

## Quick ordiplot to visualize
ordplot2 <- ordiplot(capmodel2, type="text", display="all")
ordisymbol(ordplot2, preds_cap, "Kelpdom", legend=FALSE, colors=TRUE, col=colsmap)
with(preds_cap, ordiellipse(capmodel2, Kelpdom, col=colsmap2, kind = "se", conf=0.95, label=TRUE, lwd = 1.75))
legend("bottomleft", legend=levels(preds_cap$Kelpdom), bty="n", pch=pchs, col = colsmap)


## Clean ggplots 
# https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html

# Extracting the axes data from the model
axis.long2 <- axis.long(capmodel2, choices=c(1,2))
axis.long2
# Extracting the locations of species from the model
species.long2 <- species.long(ordplot2)
head(species.long2)
# Extracting the locations of sites from the model
sites.long2 <- sites.long(ordplot2, env.data=preds_cap)
head(sites.long2)
# Extracting the locations of centroids from the sites.long output
centroids.long2 <- centroids.long(sites.long2, grouping=preds_cap$Kelpdom, centroids.only=TRUE)
head(centroids.long2)
# Extracting env biplot correlations from the model
envfit2 <- envfit(ord=ordplot2, env=preds_cap)
envvectors <- c("Tempave", "Depth_datum_m", "exp_36") # Defining env vectors of interest
vectorfit.long2 <- vectorfit.long(envfit2) %>%
  filter(vector %in% envvectors)
head(vectorfit.long2)

# Extracting the CI ellipses from the ordiplot above
ordsimple2 <- ordiplot(capmodel2)
ordellipses2 <- with(preds_cap, ordiellipse(ordsimple2, groups=Kelpdom, kind = "se", conf=0.95))
ellipses2 <- ordiellipse.long(ordellipses2, grouping.name="Kelpdom")

# Extracting the most influential spp from the model and ordiplot
species.envfit2 <- envfit(ordsimple2, env=taxahel) # The simple ordiplot above & Hellinger transformed spp data (see previous section)
species.envfit.table2 <- data.frame(r=species.envfit2$vectors$r, p=species.envfit2$vectors$pvals)
species.long.var2 <- species.long(ordsimple2, spec.data=species.envfit.table2)
species.long.var2
# Selecting for spp that explain > 60% of the observed variance among groups
species.long.vargrt2 <- species.long.var2[species.long.var2$r >= 0.5, ]
species.long.vargrt2

# The ggplot
plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long[1, "label"]) +
  ylab(axis.long[2, "label"]) +  
  scale_y_continuous(limits=c(-2,2)) + 
  scale_x_continuous(limits=c(-2,2)) + 
  geom_mark_ellipse(data=sites.long2, # Ellipses hulls for all group points
                    aes(x=axis1, y=axis2, colour=Kelpdom,
                        fill=after_scale(alpha(colour, 0.2))),
                    expand=0, linewidth=0.2, show.legend=FALSE) +
  # geom_polygon(data=ellipses2, # 95% confidence ellipses
  #              aes(x=axis1, y=axis2, colour=Kelpdom,
  #                  fill=after_scale(alpha(colour, 0.2))),
  #              size=0.2, show.legend=FALSE) +
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
               aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
               colour="black", size=0.7, arrow=arrow()) +
  geom_text_repel(data=species.long.vargrt2, # Adding in the species labels
                  aes(x=axis1*4, y=axis2*4, label=labels),
                  colour="black") +
  scale_shape_manual(values=pchs) +
  scale_color_manual(values=rev(colsmap)) +
  theme_minimal()
plotgg2


#### Community data stats (SIMPER, Permanova) ----

library(pairwiseAdonis)

## Setting up the matrices

# Hellinger transformation of the raw species abundance matrix
taxahel <- decostand(sqrt(taxamat), method = "hellinger")
# Bray curtis transformation of the transformed relative abundance matrix
taxabray <- vegan::vegdist(wisconsin(taxarel^(1/4)), method="bray")


## Running PERMANOVA /or/ general anova-like permutation test (anova.cca())

## HELLINGER
# Cluster grouping
permh <- adonis2(sqrt(taxamat) ~ Cluster+exp_36+Tempave+Depth_datum_m, data=preds_cap, permutations = 999, method="hellinger")
permh_pair <- pairwise.adonis2(sqrt(taxamat) ~ Cluster+exp_36+Tempave+Depth_datum_m, data=preds_cap)
permh_anova <- anova.cca(capmodel, by="term") # Generalized anova-like permutation test

# Kelp spp grouping
permh2 <- adonis2(sqrt(taxamat) ~ Kelpdom+exp_36+Tempave+Depth_datum_m, data=preds_cap, permutations = 999, method="hellinger", by="term")
permh2_pair <- pairwise.adonis2(sqrt(taxamat) ~ Kelpdom, data=preds_cap)
permh2_anova <- anova.cca(capmodel2, by="term") # Generalized anova-like permutation test


# ## BRAY CURTIS
# # Cluster grouping
# permb <- adonis2(wisconsin(taxarel^(1/4)) ~ Cluster, data=preds_cap, permutations = 999, method="bray")
# permb_pair <- pairwise.adonis2(wisconsin(taxarel^(1/4)) ~ Cluster, data=preds_cap)
# permb_anova <- anova.cca(capmodel, by="term")
# 
# # Kelp spp grouping
# permb2 <- adonis2(wisconsin(taxarel^(1/4)) ~ Kelpdom, data=preds_cap, permutations = 999, method="bray", by="term")
# permb2_pair <- pairwise.adonis2(wisconsin(taxarel^(1/4)) ~ Kelpdom, data=preds_cap)
# permb2_anova <- anova.cca(capmodel2, by="term")


## Running SIMPER
simp <- simper()




