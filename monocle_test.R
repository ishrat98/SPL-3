library(monocle)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd)

cth <- newCellTypeHierarchy()

MYF5_id <- row.names(subset(fData(cds), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(cds), gene_short_name == "ANPEP"))

cth <- addCellType(cth, "Myoblast", classify_func =
                     function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func =
                     function(x) { x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 } )

cds <- classifyCells(cds, cth, 0.1)

cds <- clusterCells(cds)

disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Media")
sig_genes <- subset(diff_test_res, qval < 0.1)