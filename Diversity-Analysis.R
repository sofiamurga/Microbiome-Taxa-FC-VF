# Diversity analysis metaphlan3 (through HUMAnN3) and virulence factors (through abricate-VFDB)
# HEX colors NW "#D3C331", OWOB "#0AA0A0"

library(phyloseq)
library(vegan)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(scales)

# Import feature data
feature <- read_csv("~/FeatureDataframe.csv")
# Import metadata
meta <- read_csv("~/metadata.csv")

# order data
row.names(feature) <- feature$SampleID
row.names(meta) <- meta$SampleID
feature = feature[order(row.names(feature)),]
meta = meta[order(row.names(meta)),]

# remove columns with all 0
feature <- feature[, !apply(feature == 0, 2, all)]

# Proportions (change columns)
featureProp <- sweep(feature[,2:ncol(feature)],1,rowSums(feature[,2:ncol(feature)]),"/")

# remove rare taxa (>0.001%)
## apply max function to columns and return those greater than 0.0001 
feature.abund=featureProp[,which(apply(featureProp,2,max)>0.0001)]

# Graph taxa into barplots
featureRelAb <- cbind(feature[,"SampleID",drop=F], featureProp)
featureRelAb_melt <- melt(featureRelAb, id.vars = "SampleID")
# plot Phylum stacked relative abundance
ggplot(featureRelAb_melt, aes(x = SampleID, y = value, fill = variable)) + 
  geom_bar(position = "fill", stat = "identity") + scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Relative Abundance") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 19), axis.text.x = element_text(size = 12))  

# Set seed
set.seed(8765)

# Alpha diversity
Phlan.physeq = otu_table(as.matrix(feature.abund[,1:ncol(feature.abund)]), taxa_are_rows = FALSE)
meta.physeq = sample_data(meta)

physeq.alpha = phyloseq(Phlan.physeq, meta.physeq)
sample_data(physeq.alpha)$shannon.physeq_feature <- estimate_richness(physeq.alpha, measures = "Shannon")
plot_richness(physeq.alpha, "percentil", measures = "Shannon")
plot_richness(physeq.alpha, "percentil", measures = "Simpson")

meta$shannon.vegan_feature <- diversity(feature.abund[,1:ncol(feature.abund)], index = "shannon")
meta$simpson.vegan_feature <- diversity(feature.abund[,1:ncol(feature.abund)], index = "simpson")
meta$invsimpson.vegan_feature <- diversity(feature.abund[,1:ncol(feature.abund)], index = "invsimpson")

hist(meta$shannon.vegan_feature, main="Shannon diversity", xlab="", breaks=10)
hist(meta$simpson.vegan_feature, main="Simpson diversity", xlab="", breaks=10)
hist(meta$invsimpson.vegan_feature, main="InvSimpson diversity", xlab="", breaks=10)

# Normalcy test 
shapiro.test(meta$shannon.vegan_feature)
shapiro.test(meta$simpson.vegan_feature)
shapiro.test(meta$invsimpson.vegan_feature)

pairwise.wilcox.test(meta$shannon.vegan_feature, meta$percentil, p.adjust.method = "fdr") #species and genera
pairwise.wilcox.test(meta$simpson.vegan_feature, meta$percentil, p.adjust.method = "fdr") #species and genera
pairwise.t.test(meta$invsimpson.vegan_feature, meta$percentil) #species
pairwise.wilcox.test(meta$invsimpson.vegan_feature, meta$percentil, p.adjust.method = "fdr") #genera


### PLOTS alpha 
meta %>%
  filter(percentil %in% c("NW", "OWOB")) %>%
  ggplot(aes(x = percentil, y = shannon.vegan_feature, fill = factor(percentil))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D", "#01BFC4")) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), na.rm = TRUE) +
  labs(y = "shannon alpha diversity") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 19), axis.text.x = element_text(size = 22)) 
meta %>%
  filter(percentil %in% c("NW", "OWOB")) %>%
  ggplot(aes(x = percentil, y = invsimpson.vegan_feature, fill = factor(percentil))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D", "#01BFC4")) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), na.rm = TRUE) +
  labs(y = "simpson inversed alpha diversity") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 19), axis.text.x = element_text(size = 22)) 

##################
#####################
#########################
############################
#################################
####################################
########################################
##############################################
## BETA diversity

#### Bray-Curtis: stress values <0.2
BC.nmds = metaMDS(feature.abund[,1:ncol(feature.abund)], distance="bray", k=2, trymax=10000)#### PERMANOVA 
BC.dist = vegdist(feature.abund[,1:ncol(feature.abund)], distance = "bray")
adonis(BC.dist ~ meta$percentil, data = meta, permutations = 10000)
# Plot nMDS
plot(BC.nmds, type="n", main="Bray-Curtis")
ordiellipse(BC.nmds, groups=meta$percentil, display="sites", kind="se", label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("NW"), border=FALSE)
ordiellipse(BC.nmds, groups=meta$percentil, display="sites", kind="se", label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("OWOB"), border=FALSE)


##### Jaccard: stress values <0.2
J.nmds = metaMDS(feature.abund[,1:ncol(feature.abund)], distance="jaccard", k=2, trymax=10000)
##### PERMANOVA 
J.dist = vegdist(feature.abund[,1:107], distance="jaccard")
adonis(J.dist ~ meta$percentil, data = meta, permutations = 10000)
# Plot nMDS
plot(J.nmds, type="n", main="Jaccard")
ordiellipse(J.nmds, groups=meta$percentil, display="sites", kind="se", label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("NW"), border=FALSE)
ordiellipse(J.nmds, groups=meta$percentil, display="sites", kind="se", label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("OWOB"), border=FALSE)


##### Weighted UniFrac: stress values <0.2
wUF.dist <- read_csv("~/unifrac_metaphlanbugs.csv")
# take out column ID and add row.names
wUF.distM <- wUF.dist[,-1]
row.names(wUF.distM) <- wUF.dist$SampleID
##### PERMANOVA 
adonis(wUF.distM ~ meta$percentil, permutations = 1000)
# Plot nMDS
wUF.nmds = metaMDS(wUF.distM, k=2, trymax=10000)
plot(wUF.nmds, type="n", main="Weighted UniFrac")
ordiellipse(wUF.nmds, groups=meta$percentil, display="sites", kind="se", label=FALSE, col="pink", draw="polygon", alpha=200, show.groups = c("NW"), border=FALSE)
ordiellipse(wUF.nmds, groups=meta$percentil, display="sites", kind="se", label=FALSE, col="cyan", draw="polygon", alpha=200, show.groups = c("OWOB"), border=FALSE)

### Beta dispersion
disp.imc = betadisper(BC.dist, meta$percentil)
permutest(disp.imc, pairwise=TRUE, permutations = 10000)


### SIMPER’s output is a list of OTUs which cumulatively explain 70%+ of the variation between each comparison. The numbers below the OTUs are cumulative, so to get each OTU’s contribution, you must subtract the previous OTU’s value.
simper(feature.abund[,1:ncol(feature.abund)], meta$percentil, permutations=999)
#########################################################################
############## SIMPER's output pairwise comparison per NW vs. OWOB
pairwise.wilcox.test(feature[,2:ncol(feature)]$"", meta$percentil, p.adjust.method="fdr") #add taxa
boxplot(feature[,2:ncol(feature)]$"" ~ meta$percentil, ylab= "Relative abundance", main="Taxa", xlab= "", col =  c("#F8766D", "#01BFC4")) #add taxa
meta$Eubacterium <- feature$`k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium`


### Vectors ##############  
# Bray-Curtis beta diversity
fit.BC = envfit(BC.nmds, meta)
fit.BC
# weighted UniFrac beta diversity
fit.wUF = envfit(wUF.nmds, meta)
fit.wUF

### END