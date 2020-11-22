
library(scater)
library(plotly)
library(reshape2)
library(circlize)
library(sSeq)
library(shinycssloaders)
library(SingleCellExperiment)
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
library(ouija)
set.seed(1)


deng_SCE <- readRDS("data/deng/deng-reads.rds")

deng_SCE$cell_type2 <- factor(
  deng_SCE$cell_type2,
  levels = c("zy", "early2cell", "mid2cell", "late2cell",
             "4cell", "8cell", "16cell", "earlyblast",
             "midblast", "lateblast")
)
cellLabels <- deng_SCE$cell_type2
deng <- counts(deng_SCE)
colnames(deng) <- cellLabels

deng_SCE <- scater::runPCA(deng_SCE,ncomponent = 5)

## change color Palette with library(Polychrome)

set.seed(723451) # for reproducibility
my_color <- createPalette(10, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(deng_SCE$cell_type2))

pca_df <- data.frame(PC1 = reducedDim(deng_SCE,"PCA")[,1],
                     PC2 = reducedDim(deng_SCE,"PCA")[,2],
                     cell_type2 = deng_SCE$cell_type2)

ggplot(data = pca_df)+geom_point(mapping = aes(x = PC1, y = PC2, colour = cell_type2))+
  scale_colour_manual(values = my_color)+theme_classic()