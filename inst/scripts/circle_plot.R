# Create Basic Circos Plot for Protein x CNV
# Lizhong Liu and Gabriel Odom
# 2018-12-17
# Last updated: 2019-05-13



######  1. Input Data  ########################################################
# setwd("~/Dropbox (BBSR)/Ban and Odom - Bioconductor Package/Example Data/Xena Multi-Omics Ovarian")
library(tidyverse)
library(pathwayPCA)

data("ovPheno_df")
data("ovProteome_df")
data("ovCNV_df")

dataDir_path <- system.file(
  "extdata", package = "pathwayPCA", mustWork = TRUE
)
wikipathways_PC <- read_gmt(
  paste0(dataDir_path, "/wikipathways_human_symbol.gmt"),
  description = TRUE
)

# copyNumberClean_df <- read_csv("OV_surv_x_copyNumber.csv")
# proteo_df <- read_csv("OV_surv_x_PNNLproteome_imputed.csv")
#
#
# ##choose c2 pathway
# PWCs_char <- c(c2   = "c2.cp.v6.0.symbols.gmt")
# pw <- "c2"
# PWC <- read_gmt(
#   paste0("../PathwayCollections/", PWCs_char[pw])
# )
#
# ####2. check whether copynumber and Proteomics are match samples
# nrow(proteo_df)  #83 rows
# nrow(copyNumberClean_df) #564 rows
# copyNumberClean_df$Sample[1:10]
# proteo_df$Sample[1:10]
#
# ##subset "-" then convert into "."
# loci<-str_locate_all(copyNumberClean_df$Sample[1],"-")[[1]][3,1]
# name1<-str_sub(copyNumberClean_df$Sample,1,loci-1)
# name2<-str_replace_all(name1,"-",".")
# sum(name2 %in% proteo_df$Sample) #83 = nrow of proteo_df
# copyNumberClean_df$Sample <-name2
# #so they are matched samples!!!
#
# ####3. match samples
# logi<-name2 %in% proteo_df$Sample
# copyNumberSubset<-copyNumberClean_df[logi,]
# copyNumberSubset<-copyNumberSubset[order(copyNumberSubset$Sample),]
# proteo_df<-proteo_df[order(proteo_df$Sample),]



######  2. AESPCA  ############################################################

###  Copy Number  ###
copyNum_Omics <- CreateOmics(
  assayData_df = ovCNV_df,
  pathwayCollection_ls = wikipathways_PC,
  response = ovPheno_df,
  respType = "surv",
  minPathSize = 5
)

a1 <- Sys.time()
c2copyNum_aespca <- AESPCA_pVals(
  copyNum_Omics,
  numReps = 0,
  parallel = TRUE,
  numCores = 20,
  adjustment = "BH"
)
Sys.time() - a1 # 2.8 min


###  Proteome  ###
proteo_Omics <- CreateOmics(
  assayData_df = ovProteome_df,
  pathwayCollection_ls = wikipathways_PC,
  response = ovPheno_df,
  respType = "surv",
  minPathSize = 5
)

a2 <- Sys.time()
c2proteo_aespca <- AESPCA_pVals(
  proteo_Omics,
  numReps = 0,
  parallel = TRUE,
  numCores = 20,
  adjustment = "BH"
)
Sys.time() - a2 # 0.4 min

# saveRDS(c2copyNum_aespca,"lizhong_c2copyNumber_aespca.RDS")
# saveRDS(c2proteo_aespca, "lizhong_c2proteo_aespca.RDS")



######  3. Join PC Values  ####################################################
# The IL-1 signaling pathway is coded as WP195.

cnvIL1PC1_df <- getPathPCLs(c2copyNum_aespca, "WP195")$PCs %>%
  rename(cnvPC1 = V1)
protIL1PC1_df  <- getPathPCLs(c2proteo_aespca, "WP195")$PCs %>%
  rename(protPC1 = V1)

jointPC1_df <- inner_join(
  cnvIL1PC1_df,
  protIL1PC1_df,
  by = "sampleID"
)

jointPC1_df[, -1] <- scale(jointPC1_df[, -1])
jointScaledRanks_df <- jointPC1_df %>%
  mutate(cnvRank = rank(cnvPC1)) %>%
  mutate(protRank = rank(protPC1))


# ####4.select one pathway, extract genes data of this pathway
# ##consider choose "PID_IL2_PI3K_PATHWAY" pathway
# ##34 genes inside this pathway
# PWC$pathways[[which(PWC$TERMS=="PID_IL1_PATHWAY")]]
# PWC[["PID_IL1_PATHWAY"]]
#
# ## for copynumber, find the pc1 scores
# c2copyNumPC1<-getPathPCLs(c2copyNum_aespca, "PID_IL1_PATHWAY")$PCs
#
# ## for protomecs, find the pc1 scores
# c2ProteoPC1<-getPathPCLs(c2proteo_aespca, "PID_IL1_PATHWAY")$PCs
#
# ##merge data and tidy into a data frame
# # NOTE: I (Gabriel) don't think this is correct...
# all.equal(cnvIL1PC1_df$sampleID[1:83], protIL1PC1_df$sampleID)
# # Nope.
# lizhong_totalPC<-cbind(copyNumberSubset$Sample,c2copyNumPC1[,2],c2ProteoPC1[,2])
# colnames(lizhong_totalPC)<-c("sample","c2copyNumPC1","c2ProteoPC1")
#
# ##z-tranform for PC1 scores
# set.seed(5000)
# lizhong_totalPC$TransC2CopyNumPC1 <- scale(lizhong_totalPC$c2copyNumPC1)
# lizhong_totalPC$TransC2ProteoPC1 <- scale(lizhong_totalPC$c2ProteoPC1)
#
##set color scale
#
# lizhong_totalPC <- lizhong_totalPC[order(lizhong_totalPC$TransC2ProteoPC1),]
# lizhong_totalPC$proteinOrder <- order(lizhong_totalPC$TransC2ProteoPC1)
# lizhong_totalPC <- lizhong_totalPC[order(lizhong_totalPC$TransC2CopyNumPC1),]
# lizhong_totalPC$copyNumOrder <- order(lizhong_totalPC$TransC2CopyNumPC1)
#
# sum(lizhong_totalPC$TransC2CopyNumPC1>0)  #38 pos
# sum(lizhong_totalPC$TransC2CopyNumPC1<0)  #45 neg
#
# saveRDS(lizhong_totalPC,"lizhong_totalPC.RDS")



######  4. Colour Set-up  #####################################################

# Define custom colours
sum(jointPC1_df$cnvPC1 > 0)  # 42 pos, 41 negative

library(grDevices)
blues_char <- colorRampPalette(c("blue", "white"))(42)
plot(x = 1:42, y = rep(0, 42), pch = 15, col = blues_char)
reds_char <- colorRampPalette(c("white", "red"))(41)
plot(x = 1:41, y = rep(0, 41), pch = 15, col = reds_char)
palette_char <- c(blues_char, reds_char)
plot(x = 1:83, y = rep(0, 83), pch = 15, col = palette_char)

jointSortedRanks_df <- jointScaledRanks_df %>%
  arrange(cnvRank)

cnvCol <- palette_char[jointSortedRanks_df$cnvRank]
protCol <- palette_char[jointSortedRanks_df$protRank]

showpanel = function(col){
  image(
    z = matrix(seq_len(100), ncol = 1),
    col = col, xaxt = "n", yaxt = "n"
  )
}

showpanel(cnvCol)
showpanel(protCol)

######  5. Circos Plot  #######################################################
library(circlize)

nSamps <- 83
factors <- seq_len(nSamps)

circos.clear()
circos.par(
  "gap.degree" = 0,
  "cell.padding" = c(0, 0, 0, 0),
  start.degree = 360/40,
  track.margin = c(0, 0),
  "clock.wise" = FALSE
)
circos.initialize(factors = factors, xlim = c(0, 3))

# proteomics
circos.trackPlotRegion(
  ylim = c(0, 3),
  factors = factors,
  bg.col = protCol,
  track.height = 0.3
)

# copy number
circos.trackPlotRegion(
  ylim = c(0, 3),
  factors = factors,
  bg.col = cnvCol,
  track.height = 0.3
)

##write characters
suppressMessages(
  circos.trackText(
    rep(-3, nSamps),
    rep(-3.8, nSamps),
    labels = "IL-1 Signaling \n Pathway",
    factors = factors,
    col = "#2d2d2d",
    font = 2,
    adj = par("adj"),
    cex = 1.5,
    facing = "downward",
    niceFacing = TRUE
  )
)

##add legends
par(new = TRUE)
par(mar = c(15, 32, 5, 0))
plot(
  x = c(0, 0.1), y = c(0, 0.1),
  type = "n", xlab = "", ylab = "",
  axes = FALSE,
  main = "pc1 score"
)
text(x = 0.05, y = 0.09, "Red: > 0")
text(x = 0.05, y = 0.07, "Blue: < 0")
# par(new=TRUE)
# par(mar = c(0,0,13,32))
# plot(c(0,0.1),c(0,0.1),type = 'n', xlab = '', ylab = '',axes = F, main = 'outer loop:')
# par(new=TRUE)
# par(mar = c(0,0,15,32))
# plot(c(0,0.1),c(0,0.1),type = 'n', xlab = '', ylab = '',axes = F, main = 'copyNumber',cex=0.5)
# par(new=TRUE)
# par(mar = c(5,0,20,32))
# plot(c(0,0.1),c(0,0.1),type = 'n', xlab = '', ylab = '',axes = F, main = 'inner loop:')
# par(new=TRUE)
# par(mar = c(5,0,22,32))
# plot(c(0,0.1),c(0,0.1),type = 'n', xlab = '', ylab = '',axes = F, main = 'proteOmics',cex=0.5)

circos.clear()

