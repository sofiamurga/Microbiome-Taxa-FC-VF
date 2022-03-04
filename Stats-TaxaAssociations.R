# taxonomy analysis metaphlan3 (through HUMAnN3)
# HEX colors NW "#D3C331", OWOB "#0AA0A0"
## https://htmlcolorcodes.com/es/

### 

library(tidyverse)
library(gplots)
library(vegan)
library(phyloseq)
library(metagenomeSeq)
library(Maaslin2)


#####
setwd("~/WorkingDirectory/") # Change the current working directory 

##################

# import METAdata 
input_metadata = read.table(file = "~/metadata.tsv", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# import Taxa, species
input_Speciesdata = read.table(file = "~/MetaPhlan3_species.tsv", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

##############
# order data
input_metadata = input_metadata[order(row.names(input_metadata)),]
input_Speciesdata = input_Speciesdata[order(row.names(input_Speciesdata)),]

# Set seed
set.seed(8765)

##################
#####################
#########################
############################
#################################
####################################
########################################
##############################################
## Taxa BETA diversity
#### PERMANOVA ANOSIM
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
BC.dist = vegdist(input_Speciesdata[,1:247], distance = "bray")
adonis(BC.dist ~ sexo+age+ahfob_sobrep_padres+tiemtotaflib+zBMIfa, data = input_metadata, permutations = 10000)
#                     Df SumsOfSqs MeanSqs F.Model  R2    Pr(>F)  
#sexo                 1    0.4605 0.46050 1.82499 0.04000 0.0318 *
#age                  1    0.4802 0.48021 1.90307 0.04171 0.0220 *


###################
############################
################################### Maaslin2
#######################################
###############################
##### Fixed effects
# Anthropometric metadata > "ind_cint_cadera", "zBMIfa", "percentil"
# Serum metadata > "glucmgdl", "trigmgdl", "coltotalmgdl"
# Total Energy metadata > "energ_total"
# Macronutrients % metadata > "pct_en_ch", "pct_en_lip", "pct_en_prot", "pct_azucar", "pct_fibra", "pct_agsat", "pct_agmono", "pct_agpoli", "pct_agtrans"
# Pattern 1 and Pattern 2 metadata > "Pattern1", "Pattern2"

### TAXA -adjust by age and sex for all-
# CPLM (Running selected normalization method: TSS / Min samples required with min abundance for a feature not to be filtered: 4.500000 (Total samples in data: 45 -10%) / Applying z-score to standardize continuous metadata / transform method: LOG / selected analysis method: CPLM)
fit_data = Maaslin2(
  input_data = input_Speciesdata, 
  input_metadata = input_metadata, 
  output = "MaAsLin2_TaxaAssociations", 
  fixed_effects = c("age", "sexo"), # add metadata fixed effects
  analysis_method = "CPLM")

# Preprocess MaAsLin2 output to make heatmaps 
#######################################
###############################
# adjust p-values for multiple comparisons for individual MaAsLin2 association analisis: Taxa vs. Glu, TG, and CT -all adjusted by age & sex-  
#####
######## Serum
SerumGlu <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
SerumTG <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
SerumCT <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
SepSerum <- rbind(SerumGlu, SerumTG, SerumCT)

FDR <-p.adjust(SepSerum$pval, method = "fdr", n = length(SepSerum$pval))
SepSerum$FDR <- FDR
#write_csv(SepSerum, "SepSerum-FDR.csv")
#######
# Join dataframes (AssociatonTaxa-Dataframe.csv) form separated association analysis and add column X=(−log(qval)*sign(coeff)) to use it further for heatmaps comparison 
######################################
#### HEATMAPs
library(ggplot2)
library(dplyr) # easier data wrangling 
library(lubridate) # for easy date manipulation
library(ggExtra) # because remembering ggplot theme options is beyond me
library(tidyr) 
library(reshape2)

# input data frames and select columns
# X = −log(qval)*sign(coeff)
TaxaAssociationDataframe <- read.csv("~/AssociatonTaxa-Dataframe.csv")

TaxaAssociationDataframe <- TaxaAssociationDataframe %>%
  select(feature.Taxa, Effect, X)

# reshape the data frame
# variable = −log(qval)*sign(coeff)

TaxaAssociationMatrix <- melt(TaxaAssociationDataframe)
TaxaAssociationMatrix <- dcast(TaxaAssociationMatrix, feature.Taxa ~ Effect + variable)
TaxaAssociationMatrix[is.na(TaxaAssociationMatrix)] <- 0
rownames(TaxaAssociationMatrix) <- TaxaAssociationMatrix$feature.Taxa
TaxaAssociationMatrix$feature.Taxa <- NULL

############ HEATMAP ggplot
# https://www.r-graph-gallery.com/283-the-hourly-heatmap.html

ggplot(TaxaAssociationDataframe,aes(Effect,feature.Taxa,fill=X)) +
  geom_tile(color= "white",size=0.1) + 
  scale_fill_gradient2(low = ("#03B089"),
                       mid = "white",
                       high = ("#C7048C"),
                       midpoint = 0) +
  theme(panel.background = element_rect(fill = "#F9F8F8"),
    panel.grid.major = element_line(size = 1, linetype = 'solid',
                                    colour = "#F1F1F1"))

# convert back to a matrix to use Heatmap
TaxaAssociationMatrixHeatmap <- as.matrix(TaxaAssociationMatrix)

###### HEATMAP
library(ComplexHeatmap)
library(RColorBrewer)

GnWhPk <- colorpanel(250, "#03B089", "white", "#C7048C")

# variable = −log(qval)*sign(coeff)
heatmap(TaxaAssociationMatrixHeatmap, Colv = NA, scale = 'none', col = GnWhPk)


#############################################
######################################################
############################################
# Diet (macronutrients %)
#############################################
# Preprocess MaAsLin2 output to make heatmaps
#######################################
###############################
# adjust p-values for multiple comparisons for individual MaAsLin2 association analisis: Taxa vs. each macronutrient % -all adjusted by age & sex- 
#####
#####
######## Diet %
DietChPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietFibPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietLipPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietMFPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietPFPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietProtPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietSFPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietSugPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DietTFPct <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

DietPct <- rbind(DietChPct, DietFibPct, DietLipPct, DietMFPct, DietPFPct, DietProtPct, DietSFPct, DietSugPct, DietTFPct)

FDR <-p.adjust(DietPct$pval, method = "fdr", n = length(DietPct$pval))
DietPct$FDR <- FDR
#write_csv(DietPct, "DietPct-FDR.csv")
##########
# Join dataframes (AssociatonTaxaDiet-Dataframe.csv) form separated association analysis and add column X=(−log(qval)*sign(coeff)) to use it further for heatmaps comparison 
######################################
#### HEATMAPs for Diet

# input data frames and select columns
# X = −log(qval)*sign(coeff)
DietDataframe <- read.csv("AssociatonTaxaDiet-Dataframe.csv")

DietDataframe <- DietDataframe %>%
  select(feature.Taxa, Effect, X)

# reshape the data frame
# variable = −log(qval)*sign(coeff)
DietMatrix <- melt(DietDataframe)
DietMatrix <- dcast(DietMatrix, feature.Taxa ~ Effect + variable)
DietMatrix[is.na(DietMatrix)] <- 0
rownames(DietMatrix) <- DietMatrix$feature.Taxa
DietMatrix$feature.Taxa <- NULL

############ HEATMAP ggplot
# https://www.r-graph-gallery.com/283-the-hourly-heatmap.html

ggplot(DietDataframe,aes(Effect,feature.Taxa,fill=X)) +
  geom_tile(color= "white",size=0.1) + 
  scale_fill_gradient2(low = ("#03B089"),
                       mid = "white",
                       high = ("#C7048C"),
                       midpoint = 0) +
  theme(panel.background = element_rect(fill = "#F9F8F8"),
        panel.grid.major = element_line(size = 1, linetype = 'solid',
                                        colour = "#F1F1F1"))

# convert back to a matrix to use Heatmap
DietMatrixHeatmap <- as.matrix(DietMatrix)

###### HEATMAP
library(ComplexHeatmap)
library(RColorBrewer)

GnWhPk <- colorpanel(250, "#03B089", "white", "#C7048C")
# variable = −log(qval)*sign(coeff)
heatmap(DietMatrixHeatmap, Colv = NA, scale = 'none', col = GnWhPk)

### END
