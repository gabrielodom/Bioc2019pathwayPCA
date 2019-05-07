# Integrative Pathway Analysis with pathwayPCA
#### 15 March 2019

</br>

#### Gabriel Odom (<Gabriel.odom@med.miami.edu>)

#### James Ban (<Yuguang.ban@med.miami.edu>)

#### Lily Wang (<lily.wang@med.miami.edu>)

#### Xi Chen (<steven.chen@med.miami.edu>)


# Workshop Description

With the advance in high-throughput technology for molecular assays, multi-omics datasets have become increasingly available. In this workshop, we will demonstrate using the pathwayPCA package to perform integrative pathway-based analyses of multi-omics datasets. In particular, we will demonstrate through three case studies the capabilities of pathwayPCA to

- perform pathway analysis with gene selection,
- integrate multi-omics datasets to identify driver genes,
- estimate and visualize sample-specific pathway activities in ovarian cancer, and
- identify pathways with sex-specific effects in kidney cancer.


## Pre-requisites

* Basic knowledge of R syntax
* Familiarity with RStudio
* Familiarity with pathway analysis 
* Knowledge of [*Principal Components Analysis*](https://en.wikipedia.org/wiki/Principal_component_analysis) (PCA) is helpful but not required

## Workshop Participation

This will be a hands-on workshop, students will live-code with us. It would be helpful for participants to bring a laptop with RStudio and pathwayPCA package installed (https://gabrielodom.github.io/pathwayPCA/index.html#installing-the-package). 

## _R_ / _Bioconductor_ packages used

* `pathwayPCA` will be covered in depth
* [`rWikiPathways`](https://bioconductor.org/packages/release/bioc/html/rWikiPathways.html)
* [`tidyverse` package suite](https://www.tidyverse.org/) (We will make use of some data constructs and ideas from these packages, but we do not expect users to be intimately familiar with Tidyverse)

## Time outline

| Activity                                                             | Time |
|----------------------------------------------------------------------|------|
| Introduction to `pathwayPCA`                                         | 10m  |
| Case Study 1: Estimating sample-specific pathway activities          | 10m  |
| Case Study 2: Pathway analysis for multi-omics datasets              | 10m  |
| Case Study 3: Analyzing experiments with covariates and interactions | 10m  |
| Summary and Conclusion                                               | 5m   |
| Questions and Comments                                               | 5m   |


# Workshop goals and objectives

## Learning goals

* Describe PCA-based pathway analysis
* Apply the pathwayPCA workflow in typical gene expression and copy number variations experiments
* Perform integrative pathway-based analysis for multiple types of omics data jointly
* Apply pathwayPCA to experiments with complex designs, such as those with covariate and/or interaction effects

## Learning objectives

* Identify pathways significantly associated with survival, binary, or continuous outcome using the pathwayPCA approach
* Prioritize genes most likely to be driver (instead of passenger) genes using integrative analysis 
* Estimate and visualize individual gene effects within a significant pathway
* Estimate and visualize sample-specific pathway activities
* Identify pathways with significant sex-specific effects

