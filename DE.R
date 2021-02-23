library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)#library(SingleCellExperiment)
library(SummarizedExperiment)
library(HDF5Array)
# cdSk <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")
# #View(cdSk)
# 
# cdScFiltAnnot <-  as(cdSk, "SingleCellExperiment")
# 
# #View(cdScFiltAnnot)
# #install.packages("matrixStats")
# 
# palette(brewer.pal(8, "Dark2"))
# cellLabels <- cdScFiltAnnot$cellType
# cds <- counts(cdScFiltAnnot)
# colnames(cds) <- cellLabels
# counts <- counts(cdScFiltAnnot)
# 
# set.seed(5)
# cdScFiltAnnot <- runPCA(cdScFiltAnnot, ncomponents = 50)
# counts <- as.matrix(counts(cdScFiltAnnot))
# set.seed(5)
# sds <- SlingshotDataSet(cdScFiltAnnot)
# icMat <- evaluateK(counts = counts, sds = cdScFiltAnnot, k = 3:20, nGenes = 200, verbose = FALSE)
# 
# #icMat <- evaluateK(counts = counts, sds = cdScFiltAnnot, k = 3:10, 
#  #                  nGenes = 200, verbose = T)
# 
# 
# plotSmoothers(sce[[sigGeneStart]])
# 
# 
# 
# 
# cdScFiltAnnot <- slingshot(cdScFiltAnnot, clusterLabels = 'cellType',reducedDim = "PCA",
#                            allow.breaks = FALSE)
# summary(cdScFiltAnnot$slingPseudotime_1)
# lnes <- getLineages(reducedDim(cdScFiltAnnot,"PCA"),
#                     cdScFiltAnnot$cellType)
# 
# plot(reducedDims(cdScFiltAnnot)$PCA, col = my_color[as.character(cdScFiltAnnot$cellType)], 
#      pch=16, 
#      asp = 1)
# 
# 
# 
# 
# 
# 
# 
# 
# set.seed(7)
# pseudotime <- slingPseudotime(cdScFiltAnnot, na = FALSE)
# cellWeights <- slingCurveWeights(cdScFiltAnnot)
# sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
#               nknots = 6, verbose = FALSE)
# 
# 
# 
# 
# 
# startRes <- startVsEndTest(sce)
# 
# 
# 
# oStart <- order(startRes$waldStat, decreasing = TRUE)
# sigGeneStart <- names(sce)[oStart[3]]
# plotSmoothers(sce, counts, gene = sigGeneStart)


library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)


cdSk <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")

cdScFiltAnnot <-  as(cdSk, "SingleCellExperiment")

palette(brewer.pal(8, "Dark2"))
cellLabels <- cdScFiltAnnot$cellType
#cds <- counts(cdScFiltAnnot)
colnames(cds) <- cellLabels
counts <- counts(cdScFiltAnnot)
cdScFiltAnnot 

cdScFiltAnnot <- runPCA(cdScFiltAnnot, ncomponents = 50)
counts <- as.matrix(counts(cdScFiltAnnot))
# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
# data(countMatrix, package = "tradeSeq")
# counts <- as.matrix(countMatrix)
# rm(countMatrix)
data(cdScFiltAnnot, package = "tradeSeq")
#data(celltype, package = "tradeSeq")

set.seed(5)
icMat <- evaluateK(counts = counts, sds = cdScFiltAnnot, k = 3:10, 
                   nGenes = 1741 , verbose = T)


set.seed(7)
pseudotime <- slingPseudotime(cdScFiltAnnot, na = FALSE)
cellWeights <- slingCurveWeights(cdScFiltAnnot)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)


table(rowData(sce)$tradeSeq$converged)

assoRes <- associationTest(sce)
head(assoRes)

startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)

plotGeneCount(cdScFiltAnnot, counts, gene = sigGeneStart)

customRes <- startVsEndTest(sce, pseudotimeValues = c(0.1, 0.8))

endRes <- diffEndTest(sce)

o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, counts, sigGene)

plotGeneCount(cdScFiltAnnot, counts, gene = sigGene)

patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])

plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][4])

plotGeneCount(cdScFiltAnnot, counts, gene = rownames(patternRes)[oPat][4])
