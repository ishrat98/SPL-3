library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
                   nGenes = 200, verbose = T)

set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)





startRes <- startVsEndTest(sce)



oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)
