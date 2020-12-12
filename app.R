library(SingleCellExperiment)
library(HDF5Array)
library(TSCAN)
library(M3Drop)
library(monocle)
#library(destiny)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
library(Polychrome)
library(slingshot)
library(SLICER)
#library(ouija)
set.seed(1)


cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")


cdScFiltAnnot$cellType <- factor(
  cdScFiltAnnot$cellType,
  levels = c("zy", "early2cell", "mid2cell", "late2cell",
             "4cell", "8cell", "16cell", "earlyblast",
             "midblast", "lateblast")
)
cellLabels <- cdScFiltAnnot$cellType
deng <- counts(cdScFiltAnnot)
colnames(deng) <- cellLabels

cdScFiltAnnot <- scater::runPCA(cdScFiltAnnot,ncomponent = 5)

## change color Palette with library(Polychrome)

set.seed(723451) # for reproducibility
my_color <- createPalette(10, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(cdScFiltAnnot$cellType))



cdScFiltAnnot <- slingshot(cdScFiltAnnot, clusterLabels = 'cellType',reducedDim = "UMAP",
                      allow.breaks = FALSE)
summary(cdScFiltAnnot$slingPseudotime_1)
lnes <- getLineages(reducedDim(cdScFiltAnnot,"UMAP"),
                    cdScFiltAnnot$cellType)

plot(reducedDims(cdScFiltAnnot)$UMAP, col = my_color[as.character(cdScFiltAnnot$cellType)], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(cdScFiltAnnot$cellType)]),  
       fill = my_color[levels(cdScFiltAnnot$cellType)])
lines(SlingshotDataSet(cdScFiltAnnot), lwd=2, type = 'lineages', col = c("black"))


## Plotting the pseudotime inferred by slingshot by cell types

slingshot_df <- data.frame(colData(cdScFiltAnnot))

ggplot(slingshot_df, aes(x = sizeFactor, y = cellType, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)


ggplot(slingshot_df, aes(x = sizeFactor, y = cellType, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("Second Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

ggplot(slingshot_df, aes(x = sizeFactor, y = sizeFactor, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("Second Slingshot pseudotime") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

