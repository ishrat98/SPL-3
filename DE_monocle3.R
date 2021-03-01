library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle)
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

# For reproducibility
palette(brewer.pal(8, "Dark2"))
cdSk <- loadHDF5SummarizedExperiment(dir="updated app/cdScFiltAnnotHDF5", prefix="")
#cdSk <- loadHDF5SummarizedExperiment(dir="F:/SPL-3/updated app/updated app/cdScFiltAnnotHDF5", prefix="")
#View(cdSk)

cdScFiltAnnot <-  as(cdSk, "SingleCellExperiment")
cellLabels <- cdScFiltAnnot$cellType
cds <- counts(cdScFiltAnnot)
counts <- as.matrix(cds)

cdScFiltAnnot <- runPCA(cdScFiltAnnot, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results. 
pca <- reducedDim(cdScFiltAnnot, "PCA")
#Fit trajectories using

set.seed(200)
pd <- data.frame(cells = colnames(counts), cellType = cellLabels)
rownames(pd) <- pd$cells
fd <- data.frame(gene_short_name = rownames(counts))
rownames(fd) <- fd$gene_short_name
cds <- newCellDataSet(counts, phenoData = new("AnnotatedDataFrame", data = pd),
                      featureData = new("AnnotatedDataFrame", data = fd))
cds <- estimateSizeFactors(cds)
#cds <- reduceDimension(cds, max_components = 3)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5)
plot_cell_trajectory(cds)


# info <- extract_monocle_info(cds)
# sce <- fitGAM(counts = Biobase::exprs(cds),
#               cellWeights = info$cellWeights,
#               pseudotime = info$pseudotime)


sce <- fitGAM(cds, verbose = TRUE)


set.seed(22)
library(monocle3)
# Create a cell_data_set object
cds <- new_cell_data_set(counts, cell_metadata = pd,
                         gene_metadata = data.frame(gene_short_name = rownames(counts),
                                                    row.names = rownames(counts)))
# Run PCA then UMAP on the data
cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, preprocess_method = "PCA",
                        reduction_method = "UMAP")

# First display, coloring by the cell types from Paul et al
plot_cells(cds, label_groups_by_cluster = FALSE, cell_size = 1,
           color_cells_by = "cellType")

# Running the clustering method. This is necessary to the construct the graph
cds <- cluster_cells(cds, reduction_method = "UMAP")
# Visualize the clusters
plot_cells(cds, color_cells_by = "cluster", cell_size = 1)

# Construct the graph
# Note that, for the rest of the code to run, the graph should be fully connected
cds <- learn_graph(cds, use_partition = FALSE)

# We find all the cells that are close to the starting point
cell_ids <- colnames(cds)[pd$cellType ==  "Multipotent progenitors"]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

# We compute the trajectory
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
plot_cells(cds, color_cells_by = "pseudotime")



library(magrittr)
# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

# Get the root vertices
# It is the same node as above
root <- cds@principal_graph_aux$UMAP$root_pr_nodes

# Get the other endpoints
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

# For each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  # We find the cells that map along that path
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()
rownames(cellWeights) <- colnames(cds)
pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
                     nrow = ncol(cds), byrow = FALSE)
sce <- fitGAM(counts = counts,
              pseudotime = pseudotime,
              cellWeights = cellWeights)


gene_meta <- rowData(cds)
#gene_metadata must contain a column verbatim named 'gene_short_name' for certain functions.
gene_meta$gene_short_name  <- rownames(gene_meta)
cds <- new_cell_data_set(expression_data = counts(cds),
                         cell_metadata = colData(cds),
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

## cds@clusters$UMAP$clusters <- cds$cellType

## Step 5: Learn a graph
cds <- learn_graph(cds,use_partition = TRUE)

## Step 6: Order cells
#cds <- order_cells(cds, root_cells = cds )
cds <- order_cells(cds, root_cells = c("NPCs","NPCs.1","NPCs.2","NPCs.3") )
plot_cells(cds, color_cells_by="cellType", graph_label_size = 4, cell_size = 2,
           group_label_size = 6)+ scale_color_manual(values = my_color)

plot_cells(cds,  graph_label_size = 6, cell_size = 1, 
           color_cells_by="pseudotime",
           group_label_size = 6)


cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cellType",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)


cds <- order_cells(cds)
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))







pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)

ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle3 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle3 pseudotime")


##,monocle3
cds <- preprocess_cds(cds, method = "PCA", num_dim = 50)
plot_pc_variance_explained(cds)

plot_cells(cds, reduction_method = "PCA",
           color_cells_by = "cellType", group_label_size = 3.5,
           label_groups_by_cluster = FALSE) +
  scale_color_d3(palette = "category20b")


# Seed for random initiation of UMAP
set.seed(4837)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", init = "random")
plot_cells(cds, color_cells_by = "cellType", group_label_size = 3.5,
           label_groups_by_cluster = FALSE) +
  scale_color_d3(palette = "category20b")

set.seed(723451) # for reproducibility
my_color <- createPalette(14, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(cds$cellType))

pdata_cds <- pData(cds)
pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
cds$pseudotime_monocle3 <- pdata_cds$pseudotime_monocle3
ggplot(as.data.frame(pdata_cds), 
       aes(x = pseudotime_monocle3, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("monocle3 pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle3 pseudotime")




##diffusion map
#dengg<- logcounts(cdScFiltAnnot)

dengg <- as.matrix(logcounts(cdScFiltAnnot))
colnames(dengg) <- cellLabels                             
dm <- DiffusionMap(t(dengg))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = cds$cellType)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() +  scale_color_manual(values = my_color) +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()


cds$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])

ggplot(as.data.frame(colData(cds)), 
       aes(x = pseudotime_diffusionmap, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color)  + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")





##slicer

library("lle")
slicer_genes <- select_genes(t(dengg))
k <- select_k(t(dengg[slicer_genes,]), kmin = 30, kmax=60)

slicer_traj_lle <- lle(t(dengg[slicer_genes,]), m = 2, k)$Y

reducedDim(cds, "LLE") <- slicer_traj_lle

plot_df <- data.frame(slicer1 = reducedDim(cds, "LLE")[,1],
                      slicer2 = reducedDim(cds, "LLE")[,2],
                      cellType =  cds$cellType)
ggplot(data = plot_df)+geom_point(mapping = aes(x = slicer1, 
                                                y = slicer2, 
                                                color = cellType))+
  scale_color_manual(values = my_color)+ xlab("LLE component 1") +
  ylab("LLE component 2") +
  ggtitle("Locally linear embedding of cells from SLICER")+
  theme_classic()


slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")

ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)

start <- ends[1]
# Having defined a start cell we can order the cells in the estimated pseudotime.

pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
branches <- assign_branches(slicer_traj_graph, start)

pseudotime_slicer <-
  data.frame(
    Timepoint = cellLabels,
    pseudotime = NA,
    State = branches
  )
pseudotime_slicer$pseudotime[pseudotime_order_slicer] <-
  1:length(pseudotime_order_slicer)
cds$pseudotime_slicer <- pseudotime_slicer$pseudotime
# We can again compare the inferred pseudotime to the known sampling timepoints. SLICER does not provide a pseudotime value per se, just an ordering of cells.

ggplot(as.data.frame(colData(cds)), 
       aes(x = pseudotime_slicer, 
           y = cellType, colour = cellType)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = my_color) + theme_classic() +
  xlab("SLICER pseudotime (cell ordering)") +
  ylab("Timepoint") +
  theme_classic()