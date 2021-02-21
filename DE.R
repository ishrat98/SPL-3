library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle)

# For reproducibility
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(celltype, package = "tradeSeq")

set.seed(200)
pd <- data.frame(cells = colnames(counts), cellType = celltype)
rownames(pd) <- pd$cells
fd <- data.frame(gene_short_name = rownames(counts))
rownames(fd) <- fd$gene_short_name
cds <- newCellDataSet(counts, phenoData = new("AnnotatedDataFrame", data = pd),
                      featureData = new("AnnotatedDataFrame", data = fd))
cds <- estimateSizeFactors(cds)
cds <- reduceDimension(cds, max_components = 2)
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5)
plot_cell_trajectory(cds)


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
