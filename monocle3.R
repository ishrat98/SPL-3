library(monocle3)
library(SingleCellExperiment)
library(HDF5Array)
library(TSCAN)
library(M3Drop)
#library(monocle)
library(destiny)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
library(Polychrome)
library(slingshot)
library(SLICER)
##library(Seurat)
##library(gam)

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


set.seed(723451) # for reproducibility
my_color <- createPalette(10, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(deng_SCE$cell_type2))

deng_SCE <- scater::runPCA(deng_SCE,ncomponent = 5)

library(monocle)
#d <- deng_SCE[m3dGenes,]
## feature selection 
deng <- counts(deng_SCE)

m3dGenes <- as.character(
  M3DropFeatureSelection(deng)$Gene
)

gene_meta <- rowData(deng_SCE)
#gene_metadata must contain a column verbatim named 'gene_short_name' for certain functions.
gene_meta$gene_short_name  <- rownames(gene_meta)
cds <- new_cell_data_set(expression_data = counts(deng_SCE),
                         cell_metadata = colData(deng_SCE),
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

## cds@clusters$UMAP$clusters <- deng_SCE$cell_type2

## Step 5: Learn a graph
cds <- learn_graph(cds,use_partition = TRUE)

## Step 6: Order cells
cds <- order_cells(cds, root_cells = c("zy","zy.1","zy.2","zy.3") )

plot_cells(cds, color_cells_by="cell_type2", graph_label_size = 4, cell_size = 2,
           group_label_size = 6)+ scale_color_manual(values = my_color) 



plot_cells(cds,  graph_label_size = 6, cell_size = 1, 
           color_cells_by="pseudotime",
           group_label_size = 6)
pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)

ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle3 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle3 pseudotime")

##diffusion map

deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() +  scale_color_manual(values = my_color) +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()
##second part

deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color)  + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")


##Slicer

library("lle")
slicer_genes <- select_genes(t(deng))
k <- select_k(t(deng[slicer_genes,]), kmin = 30, kmax=60)
