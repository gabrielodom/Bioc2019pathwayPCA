#' @name Bioc2019pathwayPCA
#'
#' @title A Workshop on the pathwayPCA Package
#'
#' @description With the advance in high-throughput technology for molecular
#'   assays, multi-omics datasets have become increasingly available. In this
#'   workshop, we will demonstrate using the pathwayPCA package to perform
#'   integrative pathway-based analyses of multi-omics datasets. In particular,
#'   we will demonstrate through three case studies the capabilities of
#'   \code{pathwayPCA}
#'   \enumerate{
#'     \item perform pathway analysis with gene selection,
#'     \item integrate multi-omics datasets to identify driver genes,
#'     \item estimate and visualize sample-specific pathway activities in
#'       ovarian cancer, and
#'     \item identify pathways with sex-specific effects in kidney cancer.
#'   }
#'
#' @section \code{pathwayPCA} data functions:
#'   \describe{
#'     \item{\code{\link[pathwayPCA]{read_gmt}} - }{imports a \code{.gmt} file
#'       as a pathway collection}
#'     \item{\code{\link[pathwayPCA]{SE2Tidy}} - }{extracts an assay from a
#'       SummarizedExperiment object
#'       (\url{https://doi.org/10.18129/B9.bioc.SummarizedExperiment})
#'       and turns it into a ``tidy'' data frame}
#'     \item{\code{\link[pathwayPCA]{TransposeAssay}} - }{is a variant of the
#'       base \code{\link[base]{t}} function designed specifcially for data
#'       frames and tibbles. It preserves row and column names after
#'       transposition.}
#'   }
#'
#' @section \code{pathwayPCA -Omics} functions:
#'   \describe{
#'     \item{\code{\link[pathwayPCA]{CreateOmics}} - }{takes in a collection of
#'       pathways, a single -omics assay, and a clinical response data frame
#'       and creates a data object of class \code{Omics*}}
#'     \item{\code{\link[pathwayPCA]{SubsetPathwayData}} - }{can extract the
#'       pathway-specific assay values and responses for a given pathway from
#'       an \code{Omics*} object}
#'   }
#'
#' @section \code{pathwayPCA} methods:
#'   \describe{
#'     \item{\code{\link[pathwayPCA]{AESPCA_pVals}} - }{takes in an \code{Omics*}
#'       object and calculates pathway \eqn{p}-values (parametrically or non-
#'       parametrically), principal components, and loadings via AESPCA. This
#'       returns an object of class \code{aespcOut}.}
#'     \item{\code{\link[pathwayPCA]{SuperPCA_pVals}} - }{takes in an
#'       \code{Omics*} object with valid response information and calculates
#'       pathway parametric \eqn{p}-values, principal components, and loadings
#'       via SuperPCA. This returns an object of class \code{superpcOut}.}
#'   }
#'
#' @section \code{pathwayPCA} results functions:
#'   \describe{
#'     \item{\code{\link[pathwayPCA]{getPathPCLs}} - }{takes in an object of
#'       class \code{aespcOut} or \code{superpcOut} and the \code{TERMS} name
#'       of a pathway. This function extracts 1) the data frame of principal
#'       components and subject IDs for the given pathway, and 2) a data frame
#'       of sparse loadings and feature names for the given pathway.}
#'     \item{\code{\link[pathwayPCA]{getPathpVals}} - }{takes in an object of
#'       class \code{aespcOut} or \code{superpcOut} and returns a table of the
#'       \eqn{p}-values and false discovery rates for each pathway}
#'   }
#'
#' @import pathwayPCA
#' @import tidyverse
#'
#' @docType package
#' @name mvMonitoring

NULL
