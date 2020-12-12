library(monocle) # Load Monocle
library(SingleCellExperiment)
library(slingshot)
library(HDF5Array)
library(edgeR)
library(uwot)
library(mclust)
library(RColorBrewer)

cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="F:/SPL-3/updated app/cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")

dim(cdScFiltAnnot)
counts <- assays(cdScFiltAnnot)$counts

dim(counts)

rownames(counts) <- paste0('G',1:12022)
colnames(counts) <- paste0('c',1:1741)

rownames

sim <- SingleCellExperiment(assays = List(counts = counts))
sim
geneFilter <- apply(assays(sim)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sim <- sim[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)


pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


rd2 <- umap(t(log1p(assays(sim)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)


reducedDims(sim) <- SimpleList(PCA = rd1, UMAP = rd2)


cl1 <- Mclust(rd1)$classification
colData(sim)$GMM <- cl1


plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sim)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')


summary(sim$slingPseudotime_1)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]

plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col='black')

plot(reducedDims(sim)$PCA, col = brewer.pal(9,'Set1')[sim$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')





sim


sim5 <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)

sim5

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sim5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim5), lwd=2, col='black')
