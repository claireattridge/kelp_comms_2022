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

#### Loading all data sheets ----

# wave exp
rei <- read_csv("./MSc_data/Data_new/REI_2022.csv") %>%
  dplyr::select(-c(x,y)) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName), exp_scaled = scale(exp_36)) %>% # scaling REI
  mutate(exp_scaled = as.numeric(exp_scaled))

studysites <- as.vector(rei$SiteName) # pulling a vector of study site names

# temperature
temps <- read_csv("./MSc_data/Data_new/temps_2022.csv") %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  ungroup() %>%
  filter(SiteName %in% studysites) %>%
  droplevels() %>%
  dplyr::select(-SD_Tempave)

# slope or aspect (TBA)
slope <- read_csv()

# kelp metrics 
kelps <- read_csv("./MSc_data/Data_new/kelpmetrics_2022.csv") %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName)) %>%
  dplyr::select(-c(HeightSD, BiomassSD, DensitySD, MacroSD, NereoSD))
kelps <- kelps %>%
  mutate(MacroP = MacroM/DensityM)

# RLS taxa
taxa <- read_csv("./MSc_data/Data_new/RLS_2022.csv") %>%
  rename(SpeciesAb = Total) %>%
  as.data.frame() %>%
  mutate(SiteName = as.factor(SiteName))

dfs1 <- list(taxa, kelps, temps, rei)
dfs2 <- list(kelps, temps, rei)
  
comms <- reduce(dfs1, dplyr::left_join, by="SiteName")
preds <- reduce(dfs2, dplyr::left_join, by="SiteName")
preds[is.na(preds)] <- 0
  

#### Cluster analysis ----

mod <- lm(SpeciesAb ~ BiomassM + DensityM + HeightM + Area_m2 + MacroP, data=comms)
vif(mod)
# Biomass is the highly correlated with other kelp predictors in my dataset
# Once removed the VIF scores of remaining predictors are within a reasonable range

plot(BiomassM ~ DensityM, data=preds) # Most strongly correlated to forest density
plot(BiomassM ~ HeightM, data=preds)
plot(BiomassM ~ Area_m2, data=preds)
plot(BiomassM ~ MacroP, data=preds)


## PCA (for density & biomass)
pca <- prcomp(na.omit(preds[,c(3:4)]), center = TRUE, scale = TRUE)
fviz_eig(pca) # Scree plot: PC1 explains vast majority of variation (83.7%)
fviz_pca_biplot(pca) # Biplot

pc1 <- as.vector(pca$x[,1]) # Extracting the pc1 values for further use

preds$PC1 <- pc1 # Adding pc1 values to the working dataframe
  

## K MEANS

preds <- preds[-17,] # Removing second beach south as a strong outlier site

preds <- preds %>% # Scaling all predictors prior to clustering (Euc dists are sens to scale)
  mutate_at(c(2,3,4,5,6,8,9), funs(c(scale(.))))

# only using the kelp forest attributes for clustering sites
k <- kmeans(na.omit(preds[,c(4,8,9,15)]), centers=4, nstart=15) # use some guess of clusters
k$centers
table(k$cluster)


# plotting the clustered groups 
fviz_cluster(k, data=na.omit(preds[,c(4,8,9,15)]))
# Dim1 - Dens & PC1 (biomass/canopy height)
# Dim2 - Area & Macro prop.


# checking other k cluster values
rng <- 2:15 # K from 2 to 15
tries <- 100 #Run the K Means algorithm 100 times
avg.totw.ss <-integer(length(rng)) #Set up an empty vector to hold all of points

for(v in rng){ # For each value of the range variable
  v.totw.ss <- integer(tries) # Set up an empty vector to hold the 100 tries
  for(i in 1:tries){
    k.temp <- kmeans(na.omit(preds[,c(2:4,8,9)]), centers=v) # Run kmeans
    v.totw.ss[i] <-k.temp$tot.withinss # Store the total withinss
  }
  avg.totw.ss[v-1] <-mean(v.totw.ss) #Average the 100 total withinss
}

# plot to show ave within SS as K increases (when is within SS reduced as we ^ clusters)
plot(rng,avg.totw.ss,type="b", main="Total Within SS by Various K",
     ylab="Average Total Within Sum of Squares",
     xlab="Value of K")
# looks like > 5 clusters will not substantially minimize error within clusters


# dataframe out the center outcomes (to look for differences in metrics b/w the clusters)
cframe <- as.data.frame(k$centers) %>%
  mutate(Cluster = as.factor(c("C1", "C2", "C3", "C4"))) # Add as needed here


# library(FSA)
# 
# kruskal.test(HeightM ~ Cluster, data = cframe) # different
# dunnTest(HeightM ~ Cluster, data=cframe,
#          method="bh") # correction for multiple comparisons
