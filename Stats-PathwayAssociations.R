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

# import PathwayAbundance (PathwayAbundance.tsv or genefamilies-Level4EC-CPM.tsv)
input_Pathwaydata = read.delim(file = "~/PathwayAbundance.tsv", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
# Unstratify pathways if needed:
#This can also be done with with the HUMAnN 3 untiliy `humann_split_stratified_table`
#unstrat_pathways <-function(dat_path){
#  temp = dat_path[!grepl("\\|",rownames(dat_path)),]
#  return(temp)
#}
#
##############
# order data
input_metadata = input_metadata[order(row.names(input_metadata)),]
input_Pathwaydata = input_Pathwaydata[order(row.names(input_Pathwaydata)),]

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
############################################
## Pathway BETA diversity (PathwayAbundance.tsv)
#### PERMANOVA ANOSIM
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
BC.distPath = vegdist(input_Pathwaydata[,1:436], distance = "bray")
adonis(BC.distPath ~ sexo+age+ahfob_sobrep_padres+tiemtotaflib+zBMIfa, data = input_metadata, permutations = 10000)
#                     Df SumsOfSqs MeanSqs F.Model  R2    Pr(>F)  
#tiemtotaflib         1    0.4326 0.43264  5.3047 0.11067 0.008599 **

## Pathway BETA diversity (genefamilies-Level4EC-CPM.tsv)
#### PERMANOVA ANOSIM
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
BC.distPath = vegdist(input_Pathwaydata[,1:2016], distance = "bray")
adonis(BC.distPath ~ sexo+age+ahfob_sobrep_padres+tiemtotaflib+zBMIfa, data = input_metadata, permutations = 10000)
#                     Df SumsOfSqs MeanSqs F.Model  R2    Pr(>F)  
#age                  1   0.03516 0.035164 1.78748 0.04030 0.08289 .

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
#####################
##### Fixed effects
# Anthropometric metadata > "ind_cint_cadera", "zBMIfa", "percentil"
# Serum metadata > "glucmgdl", "trigmgdl", "coltotalmgdl"
# Total Energy metadata > "energ_total"
# Macronutrients % metadata > "pct_en_ch", "pct_en_lip", "pct_en_prot", "pct_azucar", "pct_fibra", "pct_agsat", "pct_agmono", "pct_agpoli", "pct_agtrans"
# Pattern 1 and Pattern 2 metadata > "Pattern1", "Pattern2"

### Pathways PathwayAbundance.tsv -adjust all by tiemtotaflib- 
### Pathways genefamilies-Level4EC-CPM.tsv -adjust all by age- 
## CPLM or LM (normalization TSS/Min samples required with min abundance for a feature not to be filtered: 4.500000/z-scores to standarize continuous metadata/ transform method LOG, analysis method LM)
fit_func = Maaslin2(
  input_data = input_Pathwaydata, 
  input_metadata = input_metadata, 
  output = "MaAsLin2_PathwayAssociations", 
  fixed_effects = c("tiemtotaflib"), # add metadata fixed effects
  min_abundance = 0.0001,
  min_prevalence = 0.1)

### END