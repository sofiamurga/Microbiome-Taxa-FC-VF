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

# import virulence factor database
input_Virulencedata <- read.delim(file = "~/virulencefactors.tsv", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

##############
# order data
input_metadata = input_metadata[order(row.names(input_metadata)),]
input_Virulencedata = input_Virulencedata[order(row.names(input_Virulencedata)),]

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
## Virulence BETA diversity (counts 0s)
#### PERMANOVA ana ANOSIM
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
BC.distPath = vegdist(input_Virulencedata[,1:53], distance = "bray")
adonis(BC.distPath ~ sexo+age+ahfob_sobrep_padres+tiemtotaflib+zBMIfa, data = input_metadata, permutations = 10000)
#VF VFID clustered   Df SumsOfSqs MeanSqs F.Model  R2    Pr(>F)  
#      NONE

##################
#####################
#########################
############################
#################################
####################################
########################################
##############################################
######
###################
############################
################################### Maaslin2
##### Fixed effects
# Anthropometric metadata > "ind_cint_cadera", "zBMIfa", "percentil"
# Serum metadata > "glucmgdl", "trigmgdl", "coltotalmgdl"
# Total Energy metadata > "energ_total"
# Macronutrients % metadata > "pct_en_ch", "pct_en_lip", "pct_en_prot", "pct_azucar", "pct_fibra", "pct_agsat", "pct_agmono", "pct_agpoli", "pct_agtrans"
# Pattern 1 and Pattern 2 metadata > "Pattern1", "Pattern2"

## CPLM or LM (normalization TSS/Min samples required with min abundance for a feature not to be filtered: 4.500000/z-scores to standarize continuous metadata/ transform method LOG, analysis method LM)
fit_func = Maaslin2(
  input_data = input_Virulencedata, 
  input_metadata = input_metadata, 
  output = "MaAsLin2_VirulenceFactorsAssociations", 
  fixed_effects = c(""),  # add metadata fixed effects
  analysis_method = "CPLM",
  min_abundance = 0.0000001,
  min_prevalence = 0.00001)
#  
#############################################
######################################################
# Preprocess MaAsLin2 output to make heatmaps
#######################################
############################################
# adjust p-values for multiple comparisons for individual MaAsLin2 association analisis: Taxa vs. each macronutrient %  
########### Macronutrients 
pct_agtrans <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_agpoli <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_agmono <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_agsat <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_fibra <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_en_lip <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_en_ch <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pct_azucar <- read.table(file = "~/all_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

VFIDMacronuts <- rbind(pct_agtrans, pct_agpoli, pct_agmono, pct_agsat, pct_fibra, pct_en_lip, pct_en_ch, pct_azucar)
FDR <-p.adjust(VFIDMacronuts$pval, method = "fdr", n = length(VFIDMacronuts$pval))
VFIDMacronuts$FDR <- FDR
#write_csv(VFIDMacronuts, "~/VFIDMacronuts-FDR.csv")

# Join dataframes (AssociatonVirulenceFactors-Dataframe.csv) form separated association analysis and add column X=(−log(qval)*sign(coeff)) to use it further for heatmaps comparison 

######################################
#### HEATMAPs for VF
library(ggplot2)
library(dplyr) # easier data wrangling 
library(lubridate) # for easy date manipulation
library(ggExtra) # because remembering ggplot theme options is beyond me
library(tidyr) 
library(reshape2)

#########Virulence Factors
#### input data frames and select columns
#### X = −log(qval)*sign(coeff)

VFIDCDataframe <- read.csv("AssociatonVirulenceFactors-Dataframe.csv")
VFIDCDataframe <- VFIDCDataframe %>%
  select(feature.Taxa, Effect, X)

VFIDCDataframe$Effect <- factor(VFIDCDataframe$Effect, levels = c("zBMIfa", "ind_cint_cadera", "OWOB",  "glucmgdl", "trigmgdl", "coltotalmgdl",  "energ_total", "Pattern1", "Pattern2", "pct_en_ch", "pct_azucar", "pct_en_lip", "pct_agsat", "pct_agtrans", "pct_agmono"))
VFIDCDataframe$feature.Taxa <- factor(VFIDCDataframe$feature.Taxa, levels = c("CNF-1 (VF0240) ", "CT (VF0128) ", "EAST1 (VF0216) ", "Heat-stable toxin (ST) (VF0211) ", "iota-toxin (VF0381) ", "HSI-I (VF0334) ", "Map (VF0195) ", "Tsh (VF0233) ", "Tir (VF0193) ", "Paa (VF0194) ", "S fimbriae (VF0222) ", "F1C fimbriae (VF0224) ", "TcpC (VF0413) ", "IroN (VF0230) "))

############ HEATMAP ggplot
# https://www.r-graph-gallery.com/283-the-hourly-heatmap.html
ggplot(VFIDCDataframe,aes(Effect,feature.Taxa,fill=X)) +
  geom_tile(color= "white",size=0.1) + 
  scale_fill_gradient2(low = ("#03B089"),
                       mid = "white",
                       high = ("#C7048C"),
                       midpoint = 0) +
  theme(panel.background = element_rect(fill = "#F9F8F8"),
        panel.grid.major = element_line(size = 1, linetype = 'solid',
                                        colour = "#F1F1F1"))

# convert back to a matrix to use Heatmap
VFIDAssociationM <- melt(VFIDCDataframe)
VFIDAssociationM <- dcast(VFIDAssociationM, feature.Taxa ~ Effect + variable)
VFIDAssociationM[is.na(VFIDAssociationM)] <- 0
rownames(VFIDAssociationM) <- VFIDAssociationM$feature.Taxa
VFIDAssociationM$feature.Taxa <- NULL
VFIDssociationHP <- as.matrix(VFIDAssociationM)

###### HEATMAP
library(ComplexHeatmap)
library(RColorBrewer)

GnWhPk <- colorpanel(250, "#03B089", "white", "#C7048C")

# variable = −log(qval)*sign(coeff)
#heatmap(AnthropMatrixHP, Colv = NA, col = coul)
heatmap(VFIDssociationHP, Colv = NA, scale = 'none', col = GnWhPk)

### END