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

#data('slingshotExample')
rd <- cdScFiltAnnot$rd
cl <- cdScFiltAnnot$cl
condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
condition[110:140] <- 'A'
ls()

plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])

sds <- cdScFiltAnnot(rd, cl)
sds

plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
lines(sds, lwd=3)
legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])