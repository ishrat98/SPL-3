library(monocle) # Load Monocle
library(SingleCellExperiment)
library(slingshot)
library(HDF5Array)
library(edgeR)


cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")
#KData <- readRDS("D:/SPL3/Single_cell_rnaseq/data/se.rds")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")


dim(cdScFiltAnnot)
dimnames(cdScFiltAnnot)
View(cdScFiltAnnot)
counts <- assays(cdScFiltAnnot)$counts
View(counts)

