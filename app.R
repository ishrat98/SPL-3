library(slingshot)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(SingleCellExperiment)

library(HDF5Array)
library(edgeR)


cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")
#KData <- readRDS("D:/SPL3/Single_cell_rnaseq/data/se.rds")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")


dim(cdScFiltAnnot)
counts <- assays(cdScFiltAnnot)$counts
seu <- CreateSeuratObject(counts) %>% 
  SCTransform() # normalize and scale
# Add cell type annotation to metadata
seu <- AddMetaData(seu, setNames(cdScFiltAnnot$labels[ind], cells_use), 
                   col.name = "cell_type")

sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, 
                 start.clus = 4, stretch = 0)


cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(seu$cell_type, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())

plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')

plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
