library(SingleCellExperiment)
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

deng_SCE <- readRDS("data/deng-reads.rds")
deng_SCE
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



deng_SCE <- slingshot(deng_SCE, clusterLabels = 'cell_type2',reducedDim = "PCA",
                      allow.breaks = FALSE)
summary(deng_SCE$slingPseudotime_1)
lnes <- getLineages(reducedDim(deng_SCE,"PCA"),
                    deng_SCE$cell_type2)

plot(reducedDims(deng_SCE)$PCA, col = my_color[as.character(deng_SCE$cell_type2)], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(deng_SCE$cell_type2)]),  
       fill = my_color[levels(deng_SCE$cell_type2)])
lines(SlingshotDataSet(deng_SCE), lwd=2, type = 'lineages', col = c("black"))


## Plotting the pseudotime inferred by slingshot by cell types

slingshot_df <- data.frame(colData(deng_SCE))

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = cell_type2, 
                         colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)


ggplot(slingshot_df, aes(x = slingPseudotime_2, y = cell_type2, 
                         colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("Second Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = slingPseudotime_2, 
                         colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("Second Slingshot pseudotime") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

