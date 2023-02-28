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

# RLS taxa
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  rename(SpeciesAb = Total) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName))

dfs1 <- list(taxa, kelps, temps, rei)
dfs2 <- list(kelps, temps, rei)
  
comms <- reduce(dfs1, dplyr::left_join, by="SiteName")
preds <- reduce(dfs2, dplyr::left_join, by="SiteName")
preds[is.na(preds)] <- 0 # Translating the NAs to zeros for clustering analysis
  

#### Cluster analysis ----

preds <- preds[-17,] # Removing second beach south as an outlier site 

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
k <- kmeans(na.omit(preds[,c(2,4,8,9)]), centers=5, nstart=15)
k$centers # looking at the centroid mean values
table(k$cluster) # looking at # sites per cluhttp://127.0.0.1:22719/graphics/plot_zoom_png?width=1236&height=681ster


cols <- c("#e4632d", "#994455", "#225555", "#4477aa", "#997700")

tiff(file="./MSc_plots/Clusters_kmeans.tiff", height = 4.5, width = 8, units = "in", res=400)

# plotting the clustered groups 
clusts <- fviz_cluster(k, data=na.omit(preds[,c(2,4,8,9)]), pointsize=2, repel=TRUE) +
  theme_classic() +
  scale_color_manual(values=cols) +
  theme(
    legend.position="none",
    axis.text = element_text(color="black", size="12"),
    axis.title = element_text(color="black", size="13", vjust=1),
    plot.margin = unit(c(0,1,0.5,1.2), "cm"), # Selecting margins
    # panel.border = element_rect(colour="black", fill=NA, linewidth=1), # Border option
    axis.line = element_line(color="black", linewidth=0.5)) + 
  ggtitle("") +
  annotate("text", label = "C1", size=6, fontface=2, x = 1.5, y = 1.6) +
  annotate("text", label = "C2", size=6, fontface=2, x = -1.2, y = -2.1) +
  annotate("text", label = "C3", size=6, fontface=2, x = -0.2, y = 1.4) +
  annotate("text", label = "C4", size=6, fontface=2, x = -2.2, y = 0.9) +
  annotate("text", label = "C5", size=6, fontface=2, x = 1.8, y = -2.0) 
# Dim1 - Dens & prop Macro (reversed)
# Dim2 - Height & area (reversed)
clusts

dev.off()


# # C5: 
# Cable Beach (Blow Hole), Less Dangerous Bay, Swiss Boy, Wizard Islet North
# # C2:
# Dodger Channel 2, Nanat Bay
# # C3:
# Ed King East Inside, Ross Islet 2, Ross Islet Slug Island, Taylor Rock, Turf Island 2, Wizard Islet South
# # C1: 
# Between Scotts and Bradys, Danvers Danger Rock, Dixon Island Back (Bay), Dodger Channel 1, Flemming 112, Flemming 114, North Helby Rock, Tzartus 116
# # C4:
# Bordelais Island, Second Beach


# dataframe of the center outcomes (to look for differences in metrics b/w the clusters)
cframe <- as.data.frame(k$centers) %>%
  mutate(Cluster = as.factor(c("C1", "C2", "C3", "C4", "C5"))) # Add as needed here


## IMAGE PLOTS FOR CLUSTERS

library(imager)
library(magick)
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html

#C1
img1 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251965.JPG')
img1bd <- img1 %>%
  magick::image_border(color="#e4632d", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C2
img2 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251967.JPG')
img2bd <- img2 %>%
  magick::image_border(color="#994455", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C3
img3 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251975.JPG')
img3bd <- img3 %>%
  magick::image_border(color="#015f60", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C4
img4 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251989.JPG')
img4bd <- img4 %>%
  magick::image_border(color="#4477aa", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

#C5
img5 <- image_read('C:/Users/Claire/Desktop/MSc/Media/BMSC2022/CMA/P7251994.JPG')
img5bd <- img5 %>%
  magick::image_border(color="#997700", "150x150") %>% # Adding cluster col border
  magick::image_scale("1200") # Scaling image to smaller size (1200px width)

imgs <- c(img1bd,img2bd,img3bd,img4bd,img5bd) # Combining all imgs to magick vector

imgstk <- magick::image_append(image_scale(imgs, "480"), stack=TRUE) # Stacking the images vertically


# # optional formats to shift to
# 
# imgras <- as.raster(imgstk) 
# 
# imggrb <- rasterGrob(imgras)





