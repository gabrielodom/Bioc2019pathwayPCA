#' Ovarian Cancer Proteomics Data
#'
#' @description Clinical Proteomic Tumor Analysis Consortium's log-ratio
#'   normalized protein abundance expression, as processed by the Pacific
#'   Northwest National Laboratory. We removed one subject due to missing
#'   survival outcome. We removed 1712 protein features because they had greater
#'   than 20%% missingness. The remaining 2.85%% missing values in this protein
#'   expression data were imputed using the Bioconductor package `impute` under
#'   default settings.
#'
#' @format A data frame containing 4763 protein expression values for 83 samples.
#'   The first column is the TCGA sample ID (\code{Sample}).
#'
#' @source Proteome (PNNL, Gene level)
#'   \url{http://linkedomics.org/data_download/TCGA-OV/}
#'
#'   Cleaning script: \code{inst/scripts/clean_multi_omics.R}
"ovProteome_df"
