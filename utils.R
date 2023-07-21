getCiliatedCds <- function(sub_cds){
  sub_ciliated_cds_pr_test_res <- graph_test(sub_cds, neighbor_graph="principal_graph", cores=1)
  return(sub_ciliated_cds_pr_test_res)
}

getAggregatedMatrix <- function(sub_cds, sub_ciliated_cds_pr_test_res, clusterIdOfCellsInChosenCluster, qValueChosen){
  sub_pr_deg_ids <- row.names(subset(sub_ciliated_cds_pr_test_res, q_value < qValueChosen))
  print('Number of genes under this p value: ')
  print(length(sub_pr_deg_ids))
  sub_gene_module_df <- find_gene_modules(sub_cds[sub_pr_deg_ids,], resolution=0.01)
  sub_cell_group_df <- tibble::tibble(cell=row.names(colData(sub_cds)),
                                      cell_group=clusterIdOfCellsInChosenCluster)
  sub_agg_mat <- aggregate_gene_expression(sub_cds, sub_gene_module_df, sub_cell_group_df)
  row.names(sub_agg_mat) <- stringr::str_c("Module ", row.names(sub_agg_mat))
  return(list(sub_gene_module_df, sub_agg_mat))
}

showAggregatedMat <- function(agg_mat){
  pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
}

plotGeneEvolution <- function(cds,gene_module_df, wanted_Modules, label = TRUE){
  if (label == FALSE){
    plot_cells(cds,
               genes=gene_module_df%>% filter(module %in% wanted_Modules),
               show_trajectory_graph=FALSE,
               label_cell_groups = FALSE,
               reduction_method = "UMAP"
    )
  }
  else {
    plot_cells(cds,
               genes=gene_module_df%>% filter(module %in% wanted_Modules),
               show_trajectory_graph=FALSE,
               label_cell_groups = TRUE,
               reduction_method = "UMAP"
    )

  }
}


organize_files_in_folders <- function(){
  dir.create(file.path(".", "WT1"), showWarnings = FALSE)
  dir.create(file.path(".", "WT2"), showWarnings = FALSE)
  dir.create(file.path(".", "WT3"), showWarnings = FALSE)
  file.move("WT1_barcodes.tsv.gz", "WT1/barcodes.tsv.gz")
  file.move("WT2_barcodes.tsv.gz", "WT2/barcodes.tsv.gz")
  file.move("WT3_barcodes.tsv.gz", "WT3/barcodes.tsv.gz")

  file.copy("features.tsv.gz", "WT1/")
  file.copy("features.tsv.gz", "WT2/")
  file.move("features.tsv.gz", "WT3/features.tsv.gz")

  file.move("WT1_matrix.mtx.gz", "WT1/matrix.mtx.gz")
  file.move("WT2_matrix.mtx.gz", "WT2/matrix.mtx.gz")
  file.move("WT3_matrix.mtx.gz", "WT3/matrix.mtx.gz")
}
