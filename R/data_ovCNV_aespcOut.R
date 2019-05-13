#' Ovarian Cancer Copy-number Variation AES-PCA Output
#'
#' @description The AES-PCA analysis of the TCGA ovarian cancer data set
#'   included in this package
#'
#' @format An object of class \code{aespcOut}. Access the principal components
#'   and corresponding loadings for each pathway from this analysis with the
#'   function \code{\link[pathwayPCA]{getPathPCLs}}. Access the table of pathway
#'   significance levels with the function
#'   \code{\link[pathwayPCA]{getPathpVals}}.
#'
#' @import pathwayPCA
#'
#' @source calculated via the package vignette
"ovCNV_aespcOut"
