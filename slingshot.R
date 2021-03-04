library(SingleCellExperiment)
library(TSCAN)
library(M3Drop)
library(monocle)
library(monocle3)
library(destiny)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
library(Polychrome)
library(slingshot)
library(SLICER)
library("lle")
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")





set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)

assoRes <- associationTest(sce)
head(assoRes)



startRes <- startVsEndTest(sce)
#We can visualize the estimated smoothers for the third most significant gene.

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)

plotGeneCount(crv, counts, gene = sigGeneStart)

