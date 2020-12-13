library(SingleCellExperiment)
library(HDF5Array)
library(TSCAN)
library(M3Drop)
library(monocle)
library(destiny)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
library(Polychrome)
library(slingshot)
library(SLICER)
library(Seurat)

cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")



cellLabels <- cdScFiltAnnot$cellType
deng <- counts(cdScFiltAnnot)
colnames(deng) <- cellLabels

#cdScFiltAnnot <- scater::runPCA(cdScFiltAnnot,ncomponent = 5)

## change color Palette with library(Polychrome)

cdScFiltAnnot
table(cdScFiltAnnot$cellType)

# Run PCA on Deng data. Use the runPCA function from the SingleCellExperiment package.
cdScFiltAnnot <- runPCA(cdScFiltAnnot, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results. 
pca <- reducedDim(cdScFiltAnnot, "PCA")

# Describe how the PCA is stored in a matrix. Why does it have this structure?
head(pca)
dim(pca)


# Add PCA data to the deng_SCE object.
cdScFiltAnnot$PC1 <- pca[, 1]
cdScFiltAnnot$PC2 <- pca[, 2]


# Plot PC biplot with cells colored by cell_type2. 
# colData(deng_SCE) accesses the cell metadata DataFrame object for deng_SCE.
# Look at Figure 1A of the paper as a comparison to your PC biplot.
ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = PC1, y = PC2, color = cellType)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")












# Read the Slingshot documentation (?slingshot) and then run Slingshot below. 
# Given your understanding of the algorithm and the documentation, what is one 
# major set of parameters we omitted here when running Slingshot?
sce <- slingshot(cdScFiltAnnot, reducedDim = 'PCA')  # no clusters

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

# Plot Slingshot pseudotime vs cell stage. 
ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = sce$slingPseudotime_1, y = cellType, 
                                             colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)


ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = sce$slingPseudotime_1, y = cellType, 
                                             colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

