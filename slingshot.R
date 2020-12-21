library(SingleCellExperiment)
library(TSCAN)
library(M3Drop)
library(monocle)
#library(monocle3)
library(destiny)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
library(Polychrome)
library(slingshot)
library(SLICER)
library("lle")
#library(ouija)
set.seed(1)


## Dataset

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

##TSCAN

procdeng <- TSCAN::preprocess(counts(deng_SCE))

colnames(procdeng) <- 1:ncol(deng_SCE)

dengclust <- TSCAN::exprmclust(procdeng, clusternum = 10)

TSCAN::plotmclust(dengclust)

dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
deng_SCE$pseudotime_order_tscan <- NA
deng_SCE$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
  dengorderTSCAN$Pseudotime

## TSCAN Psedutime

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_order_tscan, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("TSCAN pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by TSCAN pseudotime")

## slingshot

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


##Monocle2
deng <- counts(deng_SCE)

m3dGenes <- as.character(
  M3DropFeatureSelection(deng)$Gene
)

# component1
d <- deng_SCE[which(rownames(deng_SCE) %in% m3dGenes), ]
d <- d[!duplicated(rownames(d)), ]

colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = geneNames)
fd <- new("AnnotatedDataFrame", data=fd)

dCellData <- newCellDataSet(counts(d), phenoData = pd, featureData = fd)
#
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% m3dGenes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData,reduction_method = "DDRTree", pseudo_expr = 1)
dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
plot_cell_trajectory(dCellDataSet)


# Store the ordering
pseudotime_monocle2 <-
  data.frame(
    Timepoint = phenoData(dCellDataSet)$timepoint,
    pseudotime = phenoData(dCellDataSet)$Pseudotime,
    State = phenoData(dCellDataSet)$State
  )
rownames(pseudotime_monocle2) <- 1:ncol(d)
pseudotime_order_monocle <-
  rownames(pseudotime_monocle2[order(pseudotime_monocle2$pseudotime), ])

deng_SCE$pseudotime_monocle2 <- pseudotime_monocle2$pseudotime

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_monocle2, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle2 pseudotime")



#diffusion map1
deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng)) #library destiny

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() +  scale_color_manual(values = my_color) +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()


#part2

deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color)  + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")


# SLICER

library(lle)
slicer_genes <- select_genes(t(deng))
k <- select_k(t(deng[slicer_genes,]), kmin = 30, kmax=60)
slicer_traj_lle <- lle(t(deng[slicer_genes,]), m = 2, k)$Y
## finding neighbours
## calculating weights
## computing coordinates
reducedDim(deng_SCE, "LLE") <- slicer_traj_lle

plot_df <- data.frame(slicer1 = reducedDim(deng_SCE, "LLE")[,1],
                      slicer2 = reducedDim(deng_SCE, "LLE")[,2],
                      cell_type2 =  deng_SCE$cell_type2)
ggplot(data = plot_df)+geom_point(mapping = aes(x = slicer1, 
                                                y = slicer2, 
                                                color = cell_type2))+
  scale_color_manual(values = my_color)+ xlab("LLE component 1") +
  ylab("LLE component 2") +
  ggtitle("Locally linear embedding of cells from SLICER")+
  theme_classic()

slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")


ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
