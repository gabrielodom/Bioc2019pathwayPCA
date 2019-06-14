#' Ovarian Cancer Gene Expression Data
#'
#' @description Clinical Proteomic Tumor Analysis Consortium's gene expression  
#'    change at the mRNA level, in GISTIC2 log2(fpkm-uq+1) format.
#'
#' @format A data frame containing 23316 gene expression values for 381 samples.
#'   The first column is the TCGA sample IDs (\code{Sample}).
#'
#' @source SCNV (Gene level, log-ratio)
#'   \url{http://linkedomics.org/data_download/TCGA-OV/}
#'
#'   Cleaning script: \code{inst/scripts/clean_multi_omics.R}
"ovmRNA_df"