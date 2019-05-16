# Explore KIRP Data
# Gabriel Odom
# 2018-11-14
# Last edit: 2019-05-14

library(tidyverse)
library(pathwayPCA)


######  Import Data  ##########################################################
# This data is difficult to find directly. Here is the search process I used:
# Go to the Xena "Hub" website directly (I was not able to find a link from
#   the Xena Browser home page):
# https://xenabrowser.net/hub/
#   Click "TCGA hub"; this will take you to
# https://xenabrowser.net/datapages/?host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#   Click "TCGA Kidney Papillary Cell Carcinoma (KIRP)"; this will take you to
# https://xenabrowser.net/datapages/?cohort=TCGA%20Kidney%20Papillary%20Cell%20Carcinoma%20(KIRP)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443



###  Assay  ###
# For RNAseq, scroll to the "gene expression RNAseq" section; click on the
#   "IlluminaHiSeq" link to take you to the page for log2(x+1) transformed RSEM
#   normalised counts. This will take you to
# https://xenabrowser.net/datapages/?dataset=TCGA.KIRP.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#   From here, click the link next to "download", but note that you can only
#   download the data once (the link will disappear after you click it). Just
#   in case, we have included the direct download link below. Unzip the data,
#   and follow the import steps below.
# Direct Link:
# https://tcga.xenahubs.net/download/TCGA.KIRP.sampleMap/HiSeqV2.gz

kidneyAssay_df <- read_delim(
  "inst/extdata/HiSeqV2",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
# This has class "spec_tbl_df" in addition to the tibble classes. Remove this
#   extra class via subsetting
kidneyAssay_df <- kidneyAssay_df[]

kirpRNAseq_df <- TransposeAssay(kidneyAssay_df)
anyNA(kirpRNAseq_df)
rm(kidneyAssay_df)

# usethis::use_data(kirpRNAseq_df)



###  Pheno  ###
# For clinical data, scroll to the "phenotype" section; click on the
#   "Phenotypes" link to take you to
# https://xenabrowser.net/datapages/?dataset=TCGA.KIRP.sampleMap%2FKIRP_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#   From here, click the link next to "download". We have included the direct
#   download link below as well.
# Direct Link:
# https://tcga.xenahubs.net/download/TCGA.KIRP.sampleMap/KIRP_clinicalMatrix.gz

kidneyPheno_df <- read_delim(
  "inst/extdata/KIRP_clinicalMatrix",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
kidneyPheno_df %>% select(gender) %>% unique

kirpPheno_df <-
  kidneyPheno_df %>%
  select(
    Sample = sampleID,
    time = OS.time,
    status = OS,
    gender
  ) %>%
  mutate(male = (gender == "MALE")) %>%
  select(-gender) %>%
  na.omit()
rm(kidneyPheno_df)

# usethis::use_data(kirpPheno_df)
