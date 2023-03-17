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

# Plotting out the relationship
mod <- lm(BiomassM ~ DensityM, data=preds)
modpred <- predict(mod)  
plot(BiomassM ~ DensityM, data=preds)
lines(preds$DensityM, modpred)


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
  mutate_at(c(2,3,4,5,6,8,9), funs(c(scale(.))))


# optimizing for k cluster values
rng <- 2:15 # K from 2 to 15
tries <- 100 #Run the K Means algorithm 100 times
avg.totw.ss <-integer(length(rng)) #Set up an empty vector to hold all of points

for(v in rng){ # For each value of the range variable
  v.totw.ss <- integer(tries) # Set up an empty vector to hold the 100 tries
  for(i in 1:tries){
    k.temp <- kmeans(na.omit(preds[,c(2,4,8,9)]), centers=v) # Run kmeans
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
k <- kmeans(na.omit(preds[,c(2,4,8,9)]), centers=4, nstart=15)
k$centers # looking at the centroid mean values
table(k$cluster) # looking at # sites per


cols <- c("#997700", "#225555", "#4477aa", "#e4632d")

tiff(file="./MSc_plots/Clusters_kmeans4.tiff", height = 4.5, width = 8, units = "in", res=400)

# plotting the clustered groups 
clusts <- fviz_cluster(k, data=na.omit(preds[,c(2,4,8,9)]), pointsize=2, repel=TRUE) +
  theme_classic() +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme(
    legend.position="none",
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


# # C1:
# Ed King East Inside, Ross Islet 2, Ross Islet Slug Island, Taylor Rock, Turf Island 2, Wizard Islet South
# # C2: 
# Between Scotts and Bradys, Danvers Danger Rock, Dixon Island Back (Bay), Dodger Channel 1, Flemming 112, Flemming 114, North Helby Rock, Tzartus 116
# # C3: 
# Cable Beach (Blow Hole), Less Dangerous Bay, Swiss Boy, Wizard Islet North
# # C4:
# Bordelais Island, Second Beach, Dodger Channel 2, Nanat Bay


# # dataframe of the center outcomes (to look for differences in metrics b/w the clusters)
# cframe <- as.data.frame(k$centers) %>%
#   mutate(Cluster = as.factor(c("C1", "C2", "C3", "C4", "C5"))) # Add as needed here


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

library(vegan)

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
  mutate(Species = as.factor(Species), common_name = as.factor(common_name), SiteName = as.factor(site_name),abundance = as.numeric(Total), Date = as.POSIXct(Date, format="%m/%d/%Y")) %>%
  dplyr::select(-c(site_name, Total)) 

# Filtering out the site redos at Ross Islet 2, Less Dangerous Bay, & Second Beach S in Sept 2022
# Also removing the acoustic baseline site of 'Sand Town'
comms <- comms %>%
  filter(Date <= as.POSIXct("2022-09-08") & SiteName != "Sand Town")
  

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
  summarize(observs = length(common_name)) %>%
  arrange(desc(observs)) %>%
  ungroup() %>%
  as.data.frame()

# Plot for species observations - Supp Fig
pdf(file="./MSc_plots/SuppFigs/SpeciesAbCutoff.pdf", height=10, width=7)

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
  geom_hline(yintercept = 55.5, linetype=2) # Adding in the cutoff line
obvsplot

dev.off()


## With this ^^ in mind: Removing uncommon species (< 5 site obvs) from the analysis frame (n = 37)
sppkeep <- obvs %>%
  filter(observs > 5)

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
library(vegan)
library(ecodist)

# Loading the RLS comms file
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  dplyr::select(SiteName, common_name, abundance, Method) %>%
  rename(TaxaAb = abundance) %>%
  rename(CommonName = common_name) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  filter(SiteName != "Second Beach South")

# For all species in community together
taxa_all <- taxa %>%
  dplyr::select(-Method) %>% # Removing RLS method col when using all species together
  group_by(SiteName, CommonName) %>%
  summarise(TaxaAb = sum(TaxaAb)) %>%
  ungroup() %>%
  as.data.frame()
# For just RLS Method 1 species (pelagic)
taxap <- taxa %>%
  filter(Method == 1)
# For just RLS Method 2 species (benthic)
taxab <- taxa %>%
  filter(Method == 2)

# Spreading the data to wide format (CHANGE TO PELAGIC OR BENTHIC HERE)
taxawide <- taxa_all %>%
  spread(key = CommonName, value = TaxaAb)

# Making table of relative abundances
rownames(taxawide) <- taxawide$SiteName
taxamat <- as.matrix(taxawide[,-1])
taxamat[is.na(taxamat)] <- 0

# Abundance matrix to a relative abundance matrix
taxarel <- make_relative(taxamat)

# Bray-curtis dissimilarity matrix calculations (with double wisconsin transform first)
taxabray <- vegan::vegdist(wisconsin(taxarel^(1/4)), method = "bray")



##PCOA??
pcoa.meio.bray <- cmdscale(taxabray, k = 2, eig = T)

pcoa.meio.bray.plotting <- as.data.frame(pcoa.meio.bray$points)
colnames(pcoa.meio.bray.plotting) <- c("axis_1", "axis_2")
pcoa.meio.bray.plotting$site <- rownames(pcoa.meio.bray.plotting)

# Adding on a column to identify the cluster groups (see 'comm_analyses.R' for details)
pcoa.meio.bray.plotting <- pcoa.meio.bray.plotting %>%
  mutate(Cluster = case_when(site == "Between Scotts and Bradys" | site == "Danvers Danger Rock" | site == "Dixon Island Back (Bay)" | site == "Dodger Channel 1" | site == "Flemming 112" | site == "Flemming 114" | site == "North Helby Rock" | site == "Tzartus 116" ~ "C2",
                             site == "Ed King East Inside" | site == "Ross Islet 2" | site == "Ross Islet Slug Island" | site == "Taylor Rock" | site == "Turf Island 2" | site == "Wizard Islet South" ~ "C1",
                             site == "Bordelais Island" | site == "Second Beach" | site == "Dodger Channel 2" | site == "Nanat Bay" ~ "C4",
                             site == "Cable Beach (Blow Hole)" | site == "Less Dangerous Bay" | site == "Swiss Boy" | site == "Wizard Islet North" ~ "C3",
                             TRUE ~ "Aux"))


pcoa.meio.bray$eig[1]/(sum(pcoa.meio.bray$eig))
pcoa.meio.bray$eig[2]/(sum(pcoa.meio.bray$eig))


colsmap <- c("#225555", "#4477aa", "#997700", "#e4632d")

pcoa.meio.bray.plot <- ggplot(pcoa.meio.bray.plotting, aes(x = axis_1, y = axis_2, color = Cluster, shape = Cluster, label=site)) +
  geom_point(size = 3) +
  geom_text(hjust=0, vjust=0) +
  scale_color_manual(values=colsmap, name="") +
  theme_bw() + 
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1) 
  
pcoa.meio.bray.plot


##CAP??

library(BiodiversityR)

preds_cap <- preds %>%
  mutate(Cluster = case_when(SiteName == "Between Scotts and Bradys" | SiteName == "Danvers Danger Rock" | SiteName == "Dixon Island Back (Bay)" | SiteName == "Dodger Channel 1" | SiteName == "Flemming 112" | SiteName == "Flemming 114" | SiteName == "North Helby Rock" | SiteName == "Tzartus 116" ~ "C2",
                             SiteName == "Ed King East Inside" | SiteName == "Ross Islet 2" | SiteName == "Ross Islet Slug Island" | SiteName == "Taylor Rock" | SiteName == "Turf Island 2" | SiteName == "Wizard Islet South" ~ "C1",
                             SiteName == "Bordelais Island" | SiteName == "Second Beach" | SiteName == "Dodger Channel 2" | SiteName == "Nanat Bay" ~ "C4",
                             SiteName == "Cable Beach (Blow Hole)" | SiteName == "Less Dangerous Bay" | SiteName == "Swiss Boy" | SiteName == "Wizard Islet North" ~ "C3",
                                              TRUE ~ "Aux")) %>%
  mutate(Cluster = as.factor(Cluster), Composition = as.factor(Composition))


# Inputting my relative abundance (transformed) dataframe, to then apply bray-curtis distances within the CAPdiscrim() function
capmodel <- CAPdiscrim((wisconsin(taxarel^(1/4)))~Cluster, data=preds_cap, axes=2, m=0, dist="bray", add=FALSE)

plot1 <- ordiplot(capmodel, type="none")
ordisymbol(plot1, preds_cap, "Cluster", legend=TRUE)

















