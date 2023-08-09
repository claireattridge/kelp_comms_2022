## Barkley Sound RLS photo quadrats - substrate processing
## Author: Claire Attridge
## Origin date: Dec 2022

# Load base packages
library(tidyverse)
library(MetBrewer)
library(cowplot)


#### Similarity analysis ----

### Cleaning up the data

# Reading in the raw data
rawsim <- read.csv("./MSc_data/Data_new/BenthicCover_TestingSheet_CMA.csv", header=F)
names(rawsim) <- rawsim[2,] # Making second row the header names
rawsim <- rawsim[-c(1:2),] # Removing first two rows

# Cleaning up column types
rawsim <- rawsim %>%
  mutate(across(.cols=9:16, .fns=as.numeric)) %>%
  mutate(SiteName = as.factor(SiteName), TotalPts = as.factor(TotalPts), Depth_m = as.numeric(Depth_m), Date_photos = as.POSIXct(Date_photos, format="%m/%d/%Y")) %>%
  dplyr::select(-c(Processor, ImageID, ImageWidth_cm, ImageLength_cm, AreaPerPt))

# Calculating proportion cover types
propsim <- rawsim %>%
  rowwise() %>%
  mutate(Pcanopy = (Allcanopy/Total),
         Punder = (Allunderstory/Total),
         Phard = (Hardbottom/Total),
         Psand = (Sand/Total),
         Pbiol = (Biological/Total),
         Pother = (Other/Total)) %>%
  ungroup() 

# Calculating averages by pt number
avesim <- propsim %>%
  group_by(TotalPts, SiteName) %>%
  summarise(Pcanopy = mean(Pcanopy, na.rm=T),
            Punder = mean(Punder, na.rm=T),
            Phard = mean(Phard, na.rm=T),
            Psand = mean(Psand, na.rm=T),
            Pbiol = mean(Pbiol, na.rm=T),
            Pother = mean(Pother, na.rm=T)) %>%
  ungroup() %>%
  as.data.frame()
  # dplyr::select(-SiteName) # Removing the site names for further manipulations

avesim[avesim == "NaN"] <- NA # Replacing NaNs with NAs

# Flipping so that dataframe runs from more to fewer pts in row order
avesim <- avesim %>%
  mutate(TotalPts = as.numeric(as.character(TotalPts))) %>%
  arrange(desc(TotalPts)) %>%
  mutate(TotalPts = as.factor(TotalPts)) # Turning back into a factor column

# Splitting dataframe by substrate categories & spreading wider 
canopy <- avesim[,c(1,2,3)] %>%
  spread(TotalPts, Pcanopy)
under <- avesim[,c(1,2,4)] %>%
  spread(TotalPts, Punder)
hard <- avesim[,c(1,2,5)] %>%
  spread(TotalPts, Phard)
sand <- avesim[,c(1,2,6)] %>%
  spread(TotalPts, Psand)
biol <- avesim[,c(1,2,7)] %>%
  spread(TotalPts, Pbiol)
other <- avesim[,c(1,2,8)] %>%
  spread(TotalPts, Pother)

# Applying the ccf() function for correlation among columns / point categories

# canopy
ccfcanopy <- sapply(canopy, function(x) ccf(canopy$'99', x, lag.max = 0, type="correlation", plot=FALSE, na.action=na.exclude)$acf) 
ccfcanopy <- ccfcanopy[-1] %>% # Removing the irrelevant comparison to 'SiteName'
  as.data.frame() %>% setNames("Canopy_algae")
ccfcanopy <- ccfcanopy %>%
  rownames_to_column("Pts") %>%
  mutate(Pts = factor(Pts, levels=unique(Pts), ordered=TRUE))
  
# understory
ccfunder <- sapply(under, function(x) ccf(under$'99', x, lag.max = 0, type="correlation", plot=FALSE, na.action=na.exclude)$acf) 
ccfunder <- ccfunder[-1] %>% # Removing the irrelevant comparison to 'SiteName'
  as.data.frame() %>% setNames("Understory_algae")
ccfunder <- ccfunder %>%
  rownames_to_column("Pts") %>%
  mutate(Pts = factor(Pts, levels=unique(Pts), ordered=TRUE))

# hard bottom
ccfhard <- sapply(hard, function(x) ccf(hard$'99', x, lag.max = 0, type="correlation", plot=FALSE, na.action=na.exclude)$acf) 
ccfhard <- ccfhard[-1] %>% # Removing the irrelevant comparison to 'SiteName'
  as.data.frame() %>% setNames("Hard_bottom")
ccfhard <- ccfhard %>%
  rownames_to_column("Pts") %>%
  mutate(Pts = factor(Pts, levels=unique(Pts), ordered=TRUE))

# sandy bottom
ccfsand <- sapply(sand, function(x) ccf(sand$'99', x, lag.max = 0, type="correlation", plot=FALSE, na.action=na.exclude)$acf) 
ccfsand <- ccfsand[-1] %>% # Removing the irrelevant comparison to 'SiteName'
  as.data.frame() %>% setNames("Sand_bottom")
ccfsand <- ccfsand %>%
  rownames_to_column("Pts") %>%
  mutate(Pts = factor(Pts, levels=unique(Pts), ordered=TRUE))

# biological
ccfbiol <- sapply(biol, function(x) ccf(biol$'99', x, lag.max = 0, type="correlation", plot=FALSE, na.action=na.exclude)$acf) 
ccfbiol <- ccfbiol[-1] %>% # Removing the irrelevant comparison to 'SiteName'
  as.data.frame() %>% setNames("Biological")
ccfbiol <- ccfbiol %>%
  rownames_to_column("Pts") %>%
  mutate(Pts = factor(Pts, levels=unique(Pts), ordered=TRUE))

# other
ccfother <- sapply(other, function(x) ccf(other$'99', x, lag.max = 0, type="correlation", plot=FALSE, na.action=na.exclude)$acf) 
ccfother <- ccfother[-1] %>% # Removing the irrelevant comparison to 'SiteName'
  as.data.frame() %>% setNames("Other")
ccfother <- ccfother %>%
  rownames_to_column("Pts") %>%
  mutate(Pts = factor(Pts, levels=unique(Pts), ordered=TRUE))


### Plotting the similarity data

# understory algae
uplot <- ggplot() +
  geom_point(data=ccfunder, aes(x=Pts, y=Understory_algae), size=3, shape=21, fill="white", color="black") +
  scale_y_continuous(limits=c(0.6,1.0)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color="black", size="10"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(color="black", size="12")) +
  xlab("Points") + ylab(expression("Similarity")) +
  geom_hline(yintercept = 0.95, color="black") +
  geom_vline(xintercept = 4, color="red") +
  ggtitle("Understory algae")
uplot

# canopy
cplot <- ggplot() +
  geom_point(data=ccfcanopy, aes(x=Pts, y=Canopy_algae), size=3, shape=21, fill="white", color="black") +
  scale_y_continuous(limits=c(0.6,1.0)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color="black", size="10"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(color="black", size="12")) +
  xlab("Points") + ylab(expression("")) +
  geom_hline(yintercept = 0.95, color="black") +
  geom_vline(xintercept = 1, color="red") +
  ggtitle("Canopy algae")
cplot

# hard bottom
hplot <- ggplot() +
  geom_point(data=ccfhard, aes(x=Pts, y=Hard_bottom), size=3, shape=21, fill="white", color="black") +
  scale_y_continuous(limits=c(0.6,1.0)) +
  theme_classic() +
  theme(
  axis.text = element_text(color="black", size="10"),
  axis.title.x = element_text(color="black", size="12"), 
  axis.title.y = element_text(color="black", size="12")) +
  xlab("Points") + ylab(expression("Similarity")) +
  geom_hline(yintercept = 0.95, color="black") +
  geom_vline(xintercept = 6, color="red") +
  ggtitle("Hard bottom")
hplot

# sand bottom
splot <- ggplot() +
  geom_point(data=ccfsand, aes(x=Pts, y=Sand_bottom), size=3, shape=21, fill="white", color="black") +
  scale_y_continuous(limits=c(0.6,1.0)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black", size="10"),
    axis.title.x = element_text(color="black", size="12"), 
    axis.title.y = element_text(color="black", size="12")) +
  xlab("Points") + ylab(expression("")) +
  geom_hline(yintercept = 0.95, color="black") +
  geom_vline(xintercept = 4, color="red") +
  ggtitle("Sand bottom")
splot


tiff(file="./MSc_plots/SuppFigs/SubstrateSimAnalysis.tiff", height = 6, width = 7.5, units = "in", res=300)

# Putting the four plots together
togplot <- plot_grid(uplot, cplot, hplot, splot, ncol=2)
togplot

dev.off()


#### Percent cover processing ----

# Reading in the raw data
raw <- read.csv("./MSc_data/Data_new/BenthicCover_RawSheet.csv", header=F)
names(raw) <- raw[2,] # Making second row the header names
raw <- raw[-c(1:2),] # Removing first two rows

# Cleaning up column types
raw <- raw %>%
  mutate(across(.cols=10:47, .fns=as.numeric)) %>%
  mutate(SiteName = as.factor(SiteName), TotalPts = as.factor(TotalPts), Depth = as.numeric(Depth), Date_photos = as.POSIXct(Date_photos, format="%m/%d/%Y")) %>%
  dplyr::select(-c(Processor, ImageID, ImageWidth_cm, ImageLength_cm, AreaPerPt, Notes))

# sorting data categories into broader groupings
rawgrp <- raw %>%
  rowwise() %>%
  mutate(TotalPoints = sum(raw[,c(10:41)]),
         Allbrown = DesmerestiaSpp + SargassumMuticum + EgregiaMenziesii + LaminariaSetchellii + AlariaMarginata + AgarumSpp + SaccharinaLatissima + CostariaCostata + PterygophoraCalifornica + PleurophycusGardneri + FucusSpp + HaplogloiaAndersonii + AnalipusJaponicus + UnknownBrown,
         Allred = ChondracanthusExasperatus + MazzaellaSplendens + UnknownRed,
         Allgreen = UlvaLactuca + UnknownGreen,
         Allunderstory = DesmerestiaSpp + SargassumMuticum + EgregiaMenziesii + LaminariaSetchellii + AlariaMarginata + AgarumSpp + SaccharinaLatissima + CostariaCostata + PterygophoraCalifornica + PleurophycusGardneri + FucusSpp + HaplogloiaAndersonii + AnalipusJaponicus + UnknownBrown + ChondracanthusExasperatus + MazzaellaSplendens + UnknownRed + UlvaLactuca + UnknownGreen,
         Allcanopy = MacrocystisPyrifera + NereocystisLuetkeana,
         Allcoralline = CorallineStanding + CorallineStanding,
         Hardbottom = Bedrock + Boulder + Cobble + CorallineCrustose,
         Softbottom = Sand + Silt,
         Sandbottom = Sand,
         Siltbottom = Silt,
         Invasive = DesmerestiaSpp,
         Turf = TurfAlgaeMix,
         AnalyzedPoints = Total - Unknown - Object) %>%
  ungroup()

# Calculating proportion cover types
rawprop <- rawgrp %>%
  rowwise() %>%
  mutate(Pbrown = (Allbrown/AnalyzedPoints),
         Pred = (Allred/AnalyzedPoints),
         Pgreen = (Allgreen/AnalyzedPoints),
         Punderstory = (Allunderstory/AnalyzedPoints),
         Pcanopy = (Allcanopy/AnalyzedPoints),
         Pcoralline = (Allcoralline/AnalyzedPoints),
         Pharbottom = (Hardbottom/AnalyzedPoints),
         Psoftbottom = (Softbottom/AnalyzedPoints),
         Psandbottom = (Sandbottom/AnalyzedPoints),
         Psiltbottom = (Siltbottom/AnalyzedPoints),
         Panimal = (AnimalMacro/AnalyzedPoints),
         Pinvasive = (DesmerestiaSpp/AnalyzedPoints),
         Pturf = (TurfAlgaeMix/AnalyzedPoints)) %>%
  ungroup() 

# Calculating site averages by proportion (to percent) cover
rawsite <- rawprop %>%
  group_by(SiteName) %>%
  summarise(Pbrown = mean(Pbrown*100, na.rm=T),
            Pred = mean(Pred*100, na.rm=T),
            Pgreen = mean(Pgreen*100, na.rm=T),
            Punderstory = mean(Punderstory*100, na.rm=T),
            Pcanopy = mean(Pcanopy*100, na.rm=T),
            Pcoralline = mean(Pcoralline*100, na.rm=T),
            Phardbottom = mean(Pharbottom*100, na.rm=T),
            Psoftbottom = mean(Psoftbottom*100, na.rm=T),
            Psandbottom = mean(Psandbottom*100, na.rm=T),
            Psiltbottom = mean(Psiltbottom*100, na.rm=T),
            Panimal = mean(Panimal*100, na.rm=T),
            Pinvasive = mean(Pinvasive*100, na.rm=T),
            Pturf = mean(Pturf*100, na.rm=T))


# saving a .csv file of the substrate types (%s) at all sites
write.csv(rawsite, "./MSc_data/Data_new/Substrates_2022.csv", row.names=F)


# Melting data into long format
library(reshape2)
rawplot <- rawsite %>%
  melt(id.vars="SiteName")

# Pulling out key groups of interest
subs <- c("Pturf", "Punderstory", "Pcanopy", "Phardbottom", "Psoftbottom", "Panimal")
rawplot_subs <- rawplot %>%
  filter(variable %in% subs)

  
### Plotting by substrate types
  
tiff(file="./MSc_plots/SuppFigs/SubstratePercentages.tiff", height = 6, width = 8, units = "in", res=300)

# All bottom types
plotsubs <- ggplot(rawplot_subs) +
  geom_col(aes(x = SiteName, 
               y = value,
               fill = variable)) +
  ylab("Mean percent cover") +
  theme_classic() +
  scale_fill_manual(values=met.brewer("Redon"),
                      labels = c("Understory algae", "Canopy kelps", "Hard bottom", "Soft bottom", "Animal", "Turf algae"), name = NULL) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = 10),
        legend.spacing.y = unit(0.25, "cm"),
        legend.key.size = unit(0.9, "cm"),
        axis.text.y = element_text(size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12))
plotsubs

dev.off()







