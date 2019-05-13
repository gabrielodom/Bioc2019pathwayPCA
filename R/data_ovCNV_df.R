#' Ovarian Cancer Copy-number Variation Data
#'
#' @description Clinical Proteomic Tumor Analysis Consortium's copy-number
#'   change at the gene level, in GISTIC2 log ratio format.
#'
#' @format A data frame containing 24776 copy-number values for 579 samples.
#'   The first column is the TCGA sample IDs (\code{Sample}).
#'
#' @source SCNV (Gene level, log-ratio)
#'   \url{http://linkedomics.org/data_download/TCGA-OV/}
#'
#'   Cleaning script: \code{inst/scripts/clean_multi_omics.R}
"ovCNV_df"
