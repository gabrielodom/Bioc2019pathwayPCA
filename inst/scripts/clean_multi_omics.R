# LinkedOmics Firehose Check
# Gabriel Odom
# 2019-05-13

library(tidyverse)
library(pathwayPCA)

# Download the clinical annotation, PNNL proteome data for tumor samples
#   (log-ratio normalized), and gene-level copy-number change (GISTIC2 log
#   ratio) from:
# http://linkedomics.org/data_download/TCGA-OV/



######  Clinical Data  ########################################################
# This file is named a much longer "firehose" name:
# Human__TCGA_OV__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi,
# but I saved it under Human_TCGA_OV_Clinical_Annotation.tsi instead.

linkedOmics_TCGAOVPheno_df <- read_delim(
  "inst/extdata/Human_TCGA_OV_Clinical_Annotation.tsi",
  "\t", escape_double = FALSE, trim_ws = TRUE
)

LOPhenoT_df <- TransposeAssay(linkedOmics_TCGAOVPheno_df) %>%
  # Replace "." with "-"
  mutate(Sample = str_replace_all(Sample, "\\.", "-")) %>%
  # character to numeric and change case to lower
  mutate(years_to_birth = as.integer(years_to_birth)) %>%
  rename(tumor_purity = Tumor_purity) %>%
  mutate(tumor_purity = as.numeric(tumor_purity)) %>%
  rename(platinum_status = Platinum_status) %>%
  # Recode radiation to logical
  mutate(radiation_therapy = radiation_therapy == "yes") %>%
  # survival time and status
  rename(OS_time = overall_survival) %>%
  mutate(OS_time = as.numeric(OS_time)) %>%
  rename(OS_status = status) %>%
  mutate(OS_status = as.integer(OS_status)) %>%
  select(-overallsurvival)
# 591 obs x 10 features

LOsurv_df <- LOPhenoT_df %>%
  select(Sample, OS_time, OS_status)

ovPheno_df <- LOsurv_df[complete.cases(LOsurv_df), ]
# 565 obs x 3 features

rm(linkedOmics_TCGAOVPheno_df, LOPhenoT_df, LOsurv_df)

# usethis::use_data(ovPheno_df)



######  Proteomics  ###########################################################
# The "firehose" name for this file is:
# Human__TCGA_OV__PNNL__Proteome__Velos___QExact__01_28_2016__PNNL__Gene__CDAP_iTRAQ_UnsharedLogRatio_r2.cct
# I named it Human_TCGA_OV_PNNL_Proteome_logRatio.cct.
linkedOmics_TCGAOV_Proteomics_df <- read_delim(
  "inst/extdata/Human_TCGA_OV_PNNL_Proteome_logRatio.cct",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
# For some bizarre reason, this object is also a "spec_tbl_df"?

LOProteomicsT_df <- TransposeAssay(linkedOmics_TCGAOV_Proteomics_df) %>%
  # Replace "." with "-"
  mutate(Sample = str_replace_all(Sample, "\\.", "-"))
rm(linkedOmics_TCGAOV_Proteomics_df)



###  Check Missingness  ###
LOProtMissing_num <- sapply(LOProteomicsT_df, function(x) mean(is.na(x)))
summary(LOProtMissing_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.0000  0.0000  0.1225  0.2143  0.5952

LOProtSampMissing_num <- apply(
  LOProteomicsT_df,
  MARGIN = 1,
  function(x) mean(is.na(x))
)
summary(LOProtSampMissing_num)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04046 0.06841 0.11736 0.12251 0.15936 0.23873

# Since we have some features with a ton of missing data, we are first going to
#   match the clinical data to the proteomics data to find if any subjects are
#   also missing survival time / outcome.



###  Subset by Recorded Clinical Outcome  ###
keepProtSamples_char <- intersect(
  ovPheno_df$Sample,
  LOProteomicsT_df$Sample
)
LOProteomicsT_df <- LOProteomicsT_df %>%
  filter(Sample %in% keepProtSamples_char)
# We lose 1

# Now, we can cut any proteins with greater than 20% missingness
LOProtMissing_num <- sapply(LOProteomicsT_df, function(x) mean(is.na(x)))
summary(LOProtMissing_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.0000  0.0000  0.1227  0.2169  0.6024

LOProtTrimT_df <- LOProteomicsT_df[, which(LOProtMissing_num < 0.2)]
# Subsetting the LOProteomicsT_df object removes the "spec_tbl_df" class?

rm(
  keepProtSamples_char,
  LOProtMissing_num,
  LOProtSampMissing_num,
  LOProteomicsT_df
)
mean(is.na(LOProtTrimT_df[, -1]))
# 2.85% missingness

###  Impute Missing Data  ###
# Use kNN from Bioconductor's package impute:: with default settings
LOProtTrim_mat <-
  LOProtTrimT_df %>%
  dplyr::select(-Sample) %>%
  t
colnames(LOProtTrim_mat) <- LOProtTrimT_df$Sample

LOProtTrim_ls <- impute::impute.knn(LOProtTrim_mat)

ovProteome_df <- TransposeAssay(
  as_tibble(LOProtTrim_ls$data, rownames = "Proteins"),
  omeNames = "firstCol"
)

rm(LOProtTrimT_df, LOProtTrim_mat, LOProtTrim_ls)

# usethis::use_data(ovProteome_df)



######  Copy Number Data  #####################################################
# The "firehose" name for this file is:
# Human__TCGA_OV__BI__SCNA__SNP_6.0__01_28_2016__BI__Gene__Firehose_GISTIC2.cct
# I named it Human_TCGA_OV_SCNV_GISTIC2.cct.
linkedOmics_TCGAOV_CNV_df <- read_delim(
  "inst/extdata/Human_TCGA_OV_SCNV_GISTIC2.cct",
  "\t", escape_double = FALSE, trim_ws = TRUE
)
# This also has class "spec_tbl_df". I think it's from readr' import "specs"
ovCNV_df <- TransposeAssay(linkedOmics_TCGAOV_CNV_df) %>%
  # Replace "." with "-"
  mutate(Sample = str_replace_all(Sample, "\\.", "-"))
anyNA(ovCNV_df)

rm(linkedOmics_TCGAOV_CNV_df)

# When I subsetted the proteome data frame, I lose that "spec" class. From the
#   readr 1.3.0 NEWS, subsetting removes this class and returns a regular
#   tibble. See https://cran.r-project.org/web/packages/readr/news/news.html.
ovCNV_df <- ovCNV_df[]

# usethis::use_data(ovCNV_df)




