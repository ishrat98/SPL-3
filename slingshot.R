library(monocle) # Load Monocle
library(SingleCellExperiment)
library(slingshot)
library(HDF5Array)
library(edgeR)


cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="F:/SPL-3/updated app/cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")

dim(cdScFiltAnnot)
#dimnames(cdScFiltAnnot)
#View(cdScFiltAnnot)
counts <- assays(cdScFiltAnnot)$counts

rownames(counts) <- paste0('G',1:12022)
colnames(counts) <- paste0('c',1:1741)
sim <- SingleCellExperiment(assays = List(counts = counts))
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

rd <- cdScFiltAnnot$rd
cl <- cdScFiltAnnot$cl
dim(rd) # data representing cells in a reduced dimensional space
length(cl) # vector of cluster labels

pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

sim5 <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sim5)$PCA, col = plotcol, pch=16, asp = 1)
lines(cdScFiltAnnot(sim5), lwd=2, col='black')

## next
# reducedDims(sim) <- SimpleList(PCA = rd1, UMAP = rd2)
# 
# library(mclust, quietly = TRUE)
# cl1 <- Mclust(rd1)$classification
# colData(sim)$GMM <- cl1
# 
# library(RColorBrewer)
# plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
# 
# cl2 <- kmeans(rd1, centers = 4)$cluster
# colData(sim)$kmeans <- cl2
# 
# plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)



#  dim(cdScFiltAnnot)
#  dimnames(cdScFiltAnnot)
#  View(cdScFiltAnnot)
#  counts <- assays(cdScFiltAnnot)$counts
#  #View(counts)
# 
# 
# rd <- cdScFiltAnnot$rd
# cl <- cdScFiltAnnot$cl
# 
# View(rd)
# condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
# condition[110:140] <- 'A'
# ls()
# 
# plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
# legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])
# 
# sds <- cdScFiltAnnot(rd, cl)
# sds
# 
# plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
# lines(sds, lwd=3)
# legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])
# 
# n <- nrow(rd); L <- ncol(slingPseudotime(sds))
# noise <- runif(n, -.1,.1)
# plot(as.numeric(slingPseudotime(sds)), rep(1:L, each=n)+noise,pch=16, col = brewer.pal(9,'Set1')[condition], axes=FALSE, xlab='Pseudotime', ylab='Lineage', las=1)
# axis(1); axis(2, at=1:L, las=1)
# 
# 
# t1 <- slingPseudotime(sds, na=FALSE)[,1]
# w1 <- slingCurveWeights(sds)[,1]
# d1 <- weighted.mean(t1[condition=='A'], w1[condition=='A']) - 
#   weighted.mean(t1[condition=='B'], w1[condition=='B'])
# dist1 <- replicate(1e4, {
#   condition.i <- sample(condition)
#   weighted.mean(t1[condition.i=='A'], w1[condition.i=='A']) - 
#     weighted.mean(t1[condition.i=='B'], w1[condition.i=='B'])
# })
# 
# t1 <- slingPseudotime(sds, na=FALSE)[,1]
# w1 <- slingCurveWeights(sds)[,1]
# d1 <- weighted.mean(t1[condition=='A'], w1[condition=='A']) - 
#   weighted.mean(t1[condition=='B'], w1[condition=='B'])
# dist1 <- replicate(1e4, {
#   condition.i <- sample(condition)
#   weighted.mean(t1[condition.i=='A'], w1[condition.i=='A']) - 
#     weighted.mean(t1[condition.i=='B'], w1[condition.i=='B'])
# })
# t2 <- slingPseudotime(sds, na=FALSE)[,2]
# w2 <- slingCurveWeights(sds)[,2]
# d2 <- weighted.mean(t2[condition=='A'], w2[condition=='A']) - 
#   weighted.mean(t2[condition=='B'], w2[condition=='B'])
# dist2 <- replicate(1e4, {
#   condition.i <- sample(condition)
#   weighted.mean(t2[condition.i=='A'], w2[condition.i=='A']) - 
#     weighted.mean(t2[condition.i=='B'], w2[condition.i=='B'])
# })
# 
# layout(matrix(1:2,nrow = 1))
# hist(dist1, breaks=50, col='grey50', xlim = range(c(d1,dist1)), probability = TRUE, xlab = 'Difference of Weighted Means', main = 'Lineage 1 Permutation Test', las=1)
# abline(v=d1, col=2, lwd=2)
# legend('topright','(x,y)',legend = c('Null Distn.','Observed'), fill=c('grey50',NA), border=c(1,NA), lty=c(NA,1), lwd=c(NA,2), col=c(NA,2), merge = TRUE, seg.len = .6)
# 
# hist(dist2, breaks=50, col='grey50', xlim = range(c(d2,dist2)), probability = TRUE, xlab = 'Difference of Weighted Means', main = 'Lineage 2 Permutation Test', las=1)
# abline(v=d2, col=2, lwd=2)
# legend('topright','(x,y)',legend = c('Null Distn.','Observed'), fill=c('grey50',NA), border=c(1,NA), lty=c(NA,1), lwd=c(NA,2), col=c(NA,2), merge = TRUE, seg.len = .6)
# layout(1)

