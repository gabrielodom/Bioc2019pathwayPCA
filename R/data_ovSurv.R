#' Ovarian Cancer -Omics Data
#'
#' @description Clinical Proteomic Tumor Analysis Consortium's normalized protein
#'   abundance expression, as processed by the Pacific Northwest National
#'   Laboratory and saved under “Proteome (PNNL, Gene level)”. We removed one
#'   subject due to missing survival outcome. Missing values in this protein
#'   expression data were imputed using the Bioconductor package `impute` under
#'   default settings.
#'
#' @format A data frame containing 5162 protein expression values for 83 samples.
#'   The first three columns are the TCGA sample ID (\code{Sample}), Overall
#'   Survival time (\code{OS_time}), and event/death indicator (\code{OS_event}).
#'
#' @import pathwayPCA
#'
#' @source TCGA_Ovarian_PNNL_Proteome; Study ID: S020-3
#'   \url{https://cptac-data-portal.georgetown.edu/cptac/s/S020}
"ovSurv_df"
