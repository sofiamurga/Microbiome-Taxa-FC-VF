# Virulence factors of gut microbiota are associated with BMI and biochemical blood markers in children with obesity

## REPOSITORY DESCRIPTION

This repository contains analysis scripts and data for the study titled "Virulence factors of gut microbiota are associated with BMI and biochemical blood markers in children with obesity".

### CONTENTS

This repository contains two folders:

1. **CODE** :
  This foldder contains one shell scipt file with the code followed for the HMP Unified Metabolic Analysis Network (HUMAnN version 3.0.0, http://huttenhower.sph.harvard.edu/humann) pipeline and five essential R scripts for our analysis.
  
  - Humann_docker.sh
  - Diversity-Analysis.R
  - Stats-TaxaAssociations.R
  - Stats-PathwayAssociations.R
  - Stats-VirulenceFactorsAssociations.R
  - EdgeR-DAG.R

2. **DATASETS** : 
  This folder contians datasets resulting from the collection of participants metadata, HUMAnN 3.0 pipeline followed for annotation of high quality microbial sequences using for alignment ChocoPhlAn and UniRef90 full databases, and virulence factor screening using ABRicate v1.0.1 (https://github.com/tseemann/abricate)
  
  - MetaData.tsv
  - MetaPhlan3_generaWAL01percent.tsv
  - MetaPhlan3_speciesWAL01percent.tsv
  - PathwayAbundance.tsv
  - unifrac_metaphlanbugs.csv
  - VF_VFIDclustered_counts0.tsv
 
## CONTACT
   - Sofia M. Murga-Garrido - sofiamurgaga@gmail.com
