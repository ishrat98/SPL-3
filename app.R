library(SingleCellExperiment)
library(SummarizedExperiment)
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
#library(Seurat)
library(monocle3)
library(matrixStats)
library(Matrix)

cdSk <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")
#View(cdSk)

cdScFiltAnnot <-  as(cdSk, "SingleCellExperiment")

#View(cdScFiltAnnot)
#install.packages("matrixStats")


cellLabels <- cdScFiltAnnot$cellType
cds <- counts(cdScFiltAnnot)
colnames(cds) <- cellLabels

#cdScFiltAnnot <- scater::runPCA(cdScFiltAnnot,ncomponent = 5)

## change color Palette with library(Polychrome)

cdScFiltAnnot
table(cdScFiltAnnot$cellType)

# Run PCA on cds data. Use the runPCA function from the SingleCellExperiment package.
cdScFiltAnnot <- runPCA(cdScFiltAnnot, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results. 
pca <- reducedDim(cdScFiltAnnot, "PCA")

# Describe how the PCA is stored in a matrix. Why does it have this structure?
head(pca)
dim(pca)


# Add PCA data to the cdScFiltAnnot object.
cdScFiltAnnot$PC1 <- pca[, 1]
cdScFiltAnnot$PC2 <- pca[, 2]


#View(cdScFiltAnnot)
my_color <- createPalette(14, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(cdScFiltAnnot$cellType))

# Plot PC biplot with cells colored by cellType. 
# colData(cdScFiltAnnot) accesses the cell metadata DataFrame object for cdScFiltAnnot.
# Look at Figure 1A of the paper as a comparison to your PC biplot.
ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = PC1, y = PC2, color = cellType)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")


##Monocle2
cds <- counts(cdScFiltAnnot)

m3dGenes <- as.character(
  M3DropFeatureSelection(cds)$Gene
)

# component1
d <- cdScFiltAnnot[which(rownames(cdScFiltAnnot) %in% m3dGenes), ]
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

cdScFiltAnnot$pseudotime_monocle2 <- pseudotime_monocle2$pseudotime

ggplot(as.data.frame(colData(cdScFiltAnnot)), 
       aes(x = pseudotime_monocle2, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle2 pseudotime")




##monocle3
gene_meta <- rowData(cdScFiltAnnot)
#gene_metadata must contain a column verbatim named 'gene_short_name' for certain functions.
gene_meta$gene_short_name  <- rownames(gene_meta)
cds <- new_cell_data_set(expression_data = counts(cdScFiltAnnot),
                         cell_metadata = colData(cdScFiltAnnot),
                         gene_metadata = gene_meta)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds,num_dim = 5)
plot_pc_variance_explained(cds)


## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
## No preprocess_method specified, using preprocess_method = 'PCA'
## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## change the clusters

## cds@clusters$UMAP$clusters <- cdScFiltAnnot$cellType

## Step 5: Learn a graph
cds <- learn_graph(cds,use_partition = TRUE)

## Step 6: Order cells
cds <- order_cells(cds, root_cells = c("zy","zy.1","zy.2","zy.3") )

plot_cells(cds, color_cells_by="cellType", graph_label_size = 4, cell_size = 2,
           group_label_size = 6)+ scale_color_manual(values = my_color) 



plot_cells(cds,  graph_label_size = 6, cell_size = 1, 
           color_cells_by="pseudotime",
           group_label_size = 6)
pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)

ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle3 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle3 pseudotime")



#diffusion map1
cds <- logcounts(cdScFiltAnnot)
cellLabels <- cdScFiltAnnot$cellType
colnames(cds) <- cellLabels
dm <- DiffusionMap(t(cds)) #library destiny


tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = cdScFiltAnnot$cellType)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() +  scale_color_manual(values = my_color) +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()


#part2

cdScFiltAnnot$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])

ggplot(as.data.frame(colData(cdScFiltAnnot)), 
       aes(x = pseudotime_diffusionmap, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color)  + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")


















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


ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = slingPseudotime_2, y = cellType, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("Second Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)


##tscan

procdeng <- TSCAN::preprocess(counts(cdScFiltAnnot))

colnames(procdeng) <- 1:ncol(cdScFiltAnnot)

dengclust <- TSCAN::exprmclust(procdeng, clusternum = 14)

TSCAN::plotmclust(dengclust)

dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
cdScFiltAnnot$pseudotime_order_tscan <- NA
cdScFiltAnnot$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
  dengorderTSCAN$Pseudotime

ggplot(as.data.frame(colData(cdScFiltAnnot)), 
       aes(x = pseudotime_order_tscan, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("TSCAN pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by TSCAN pseudotime")








ggplot(as.data.frame(colData(cdScFiltAnnot)), 
       aes(x = pseudotime_order_tscan, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("TSCAN pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by TSCAN pseudotime")


proccds <- TSCAN::preprocess(counts(cdScFiltAnnot))

#colnames(proccds) <- 1:ncol(cdScFiltAnnot)

cdsclust <- TSCAN::exprmclust(proccds, clusternum = 14)

TSCAN::plotmclust(cdsclust)

cdsorderTSCAN <- TSCAN::TSCANorder(cdsclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(cdsorderTSCAN$sample_name)
cdScFiltAnnot$pseudotime_order_tscan <- NA
cdScFiltAnnot$pseudotime_order_tscan[as.numeric(cdsorderTSCAN$sample_name)] <- 
  cdsorderTSCAN$Pseudotime

cellLabels[cdsclust$clusterid == 14]

ggplot(as.data.frame(colData(cdScFiltAnnot)), 
       aes(x = pseudotime_order_tscan, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("TSCAN pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by TSCAN pseudotime")



## slingshot

cdScFiltAnnot <- slingshot(cdScFiltAnnot, clusterLabels = 'cellType',reducedDim = "PCA",
                      allow.breaks = FALSE)
summary(cdScFiltAnnot$slingPseudotime_1)
lnes <- getLineages(reducedDim(cdScFiltAnnot,"PCA"),
                    cdScFiltAnnot$cellType)

plot(reducedDims(cdScFiltAnnot)$PCA, col = my_color[as.character(cdScFiltAnnot$cellType)], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(cdScFiltAnnot$cellType)]),  
       fill = my_color[levels(cdScFiltAnnot$cellType)])
  lines(SlingshotDataSet(cdScFiltAnnot), lwd=2, type = 'lineages', col = c("black"))


## Plotting the pseudotime inferred by slingshot by cell types

slingshot_df <- data.frame(colData(cdScFiltAnnot))

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = cellType, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)


ggplot(slingshot_df, aes(x = slingPseudotime_2, y = cellType, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("Second Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = slingPseudotime_2, 
                         colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("Second Slingshot pseudotime") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)


##Heatmap

t <- cdScFiltAnnot$slingPseudotime_1

# for time, only look at the 100 most variable genes 
Y <- log1p(assay(cdScFiltAnnot,"logcounts"))

var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- gam(z ~ lo(t), data=d)
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

## Plot the top 100 genes' expression 

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]

heatdata <- assays(cdScFiltAnnot)$logcounts[topgenes, order(t, na.last = NA)]
heatclus <- cdScFiltAnnot$cellType[order(t, na.last = NA)]

heatmap(heatdata, Colv = NA,
        ColSideColors = my_color[heatclus],cexRow = 1,cexCol = 1)

########## Monocle2 ########
## Part - 1
library(monocle)
#d <- cdScFiltAnnot[m3dGenes,]
## feature selection
cds <- counts(cdScFiltAnnot)

m3dGenes <- as.character(
  M3DropFeatureSelection(cds)$Gene
)

##part -2
d <- cdScFiltAnnot[which(rownames(cdScFiltAnnot) %in% m3dGenes), ]
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

##part-3

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

cdScFiltAnnot$pseudotime_monocle2 <- pseudotime_monocle2$pseudotime

ggplot(as.data.frame(colData(cdScFiltAnnot)), 
       aes(x = pseudotime_monocle2, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle2 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle2 pseudotime")



cds <- logcounts(cdScFiltAnnot)
colnames(cds) <- cellLabels
View(cds)
dm <- DiffusionMap(t(cds))
cdScFiltAnnot

slicer_genes <- select_genes(t(cds))
k <- select_k(t(cds[slicer_genes,]), kmin = 30, kmax=60)
slicer_traj_lle <- lle(t(cds[slicer_genes,]), m = 2, k)$Y


reduceddim(cdScFiltAnnot, "lle") <- slicer_traj_lle

plot_df <- data.frame(slicer1 = reduceddim(cdScFiltAnnot, "lle")[,1],
                      slicer2 = reduceddim(cdScFiltAnnot, "lle")[,2],
                      cellType =  cdScFiltAnnot$cellType)
ggplot(data = plot_df)+geom_point(mapping = aes(x = slicer1,
                                                y = slicer2,
                                                color = cellType))+
  scale_color_manual(values = my_color)+ xlab("lle component 1") +
  ylab("lle component 2") +
  ggtitle("locally linear embedding of cells from slicer")+
  theme_classic()

