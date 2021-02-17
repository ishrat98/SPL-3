library(monocle3)
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

