library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(monocle3)
library(SeuratWrappers)
library(ff)
source(file="utils.R")

## list of cell type markers
# Guard_cell = c('GC-AT5G25980','GC-AT5G48485','GC-AT1G62480','GC-AT3G16400','GC-AT2G15830','GC-AT1G71050','GC-AT2G19810','GC-AT4G37870','GC-AT5G66400','GC-AT3G23730','GC-AT3G24140','GC-AT5G66440','GC-AT3G56620','GC-AT4G37430','GC-AT2G34655','GC-AT2G47260','GC-AT5G42970','GC-AT3G58640','GC-AT1G23170','GC-AT1G29050')
# Companion_cell = c('CC-AT1G23130','CC-AT1G67865','CC-AT1G64370','CC-AT4G19840','CC-AT2G18328','CC-AT5G18600','CC-AT1G67860','CC-AT1G67870','CC-AT5G45350','CC-AT2G32870','CC-AT5G04080','CC-AT5G22090','CC-AT4G00780','CC-AT1G07610','CC-AT4G16008','CC-AT1G06830','CC-AT2G16740','CC-AT2G30540','CC-AT4G16000','CC-AT4G15690')
# Epidermia_cell = c('EC-AT2G38540','EC-AT1G66100','EC-AT1G09310','EC-AT3G51600','EC-AT5G25610','EC-AT5G44020','EC-AT3G16370','EC-AT2G27385','EC-AT3G26450','EC-AT1G68530','EC-AT2G32690','EC-AT4G04840','EC-AT4G23670','EC-AT1G29660','EC-AT5G13930','EC-AT5G64770','EC-AT4G39330','EC-AT1G29670','EC-AT1G55260','EC-AT2G26250')
# Mesophyl_cell = c('MC-AT2G10940','MC-AT5G38430','MC-AT3G08940','MC-AT1G72610','MC-AT3G27690','MC-AT2G05070','MC-AT1G12090','MC-AT1G29910','MC-AT2G34420','MC-AT2G34430','MC-AT1G15820','MC-AT2G21330','MC-AT1G06680','MC-AT3G59400','MC-AT2G05100','MC-AT1G67090','MC-AT3G54890','MC-AT5G66570','MC-AT4G38970','MC-AT1G44575')
# Mesophyl_cell_2 = c('MC2-AT1G18740','MC2-AT1G74930','MC2-AT1G27730','MC2-AT2G44840','MC2-AT1G80840','MC2-AT3G44260','MC2-AT5G12030','MC2-AT5G12020','MC2-AT1G74450','MC2-AT4G24570','MC2-AT3G56880','MC2-AT1G71000','MC2-AT5G66650','MC2-AT4G27652','MC2-AT3G46230','MC2-AT3G12580','MC2-AT3G55980','MC2-AT4G34410','MC2-AT5G52050','MC2-AT1G07400')
# Hydathode_cell = c('HC-AT3G16670','HC-AT3G05730','HC-AT3G16660','HC-AT1G56710','HC-AT3G09330','HC-AT1G22900','HC-AT1G08090','HC-AT4G36260','HC-AT4G32950','HC-AT2G43610','HC-AT4G23550','HC-AT2G19990','HC-AT1G62510','HC-AT2G33175','HC-AT2G38940','HC-AT3G14060','HC-AT3G60700','HC-AT1G19610','HC-AT5G60910')#,'HC-AT1G08757')
# S_cell = c('SC-AT1G78370','SC-AT3G19710','SC-AT2G30860','SC-AT1G80520','SC-AT2G43100','SC-AT5G23010','SC-AT5G02380','SC-AT2G22330','SC-AT3G14990','SC-AT2G46650','SC-AT2G26690','SC-AT5G14200','SC-AT2G22860','SC-AT5G01600','SC-AT4G14040','SC-AT3G11930','SC-AT2G37170','SC-AT3G15450','SC-AT5G03610','SC-AT1G11580')

#name of the output pdf file and output rds file
name_rep = "WT_rep123"
pdf(paste('plots_', name_rep,'.pdf',sep =""), width = 10, height = 10)
list_reps <- list("rep_1", "rep_2", "rep_3")


organize_files_in_folders()
# import the data
sc.data_WT_rep1 <- Read10X(data.dir = "WT1")
sc.data_WT_rep2 <- Read10X(data.dir = "WT2")
sc.data_WT_rep3 <- Read10X(data.dir = "WT3")

# remember the columns to label each replicate
col_1 <- colnames(sc.data_WT_rep1)
col_2 <- colnames(sc.data_WT_rep2)
col_3 <- colnames(sc.data_WT_rep3)

# Build sc_data as the concatenation of each dataset
sc.data = Reduce("cbind",list(sc.data_WT_rep1, sc.data_WT_rep2,sc.data_WT_rep3))

# Create the Seurat Object
sc_full <- CreateSeuratObject(counts = sc.data, prejct = "scRNA", min.cells = 3, min.features = 200)
sc_full[['rep_label']] = 0
# assigning to the Seurat object the replicate label for each cell
sc_full[['rep_label']][col_1,] = 'rep_1'
sc_full[['rep_label']][col_2,] = 'rep_2'
sc_full[['rep_label']][col_3,] = 'rep_3'

# Extract the pt and mt

sc_full[["percent.pt"]]<-PercentageFeatureSet(sc_full, pattern="^ATCG")
sc_full[["percent.mt"]]<-PercentageFeatureSet(sc_full, pattern="^ATMG")


# Removing the useless cells.
sc_full <- subset(sc_full, percent.pt < 5)


sc_full$rep_label <- c(sc_full[['rep_label']])



# Doing a SCTransform
sc_full <- SCTransform(sc_full, vars.to.regress="percent.pt")

# PCA only on the interested features
sc_full <- RunPCA(sc_full, features = VariableFeatures(object = sc_full))

# ElbowPlot to determine the dimensionality of the dataset
ElbowPlot(sc_full)

sc_full <- RunUMAP(sc_full, dims = 1:30)


# Check which features are the most important.
sc_full <- FindVariableFeatures(sc_full, selection.method = "vst", nfeatures = 2000)



# Get the neighbors

sc_full <- FindNeighbors(sc_full,  reduction = 'pca', dims = 1:30)
# Get the clusters via SNN

sc_full <- FindClusters(sc_full, resolution =1)


saveRDS(sc_full, file = paste("sc_full", name_rep, ".rds", sep = ""))

clusters <- Idents(sc_full)

# Visiualization with UMAP, and ACP
# Visualization of cell specific expression of selected genes in selected clusters

DimPlot(sc_full, reduction = "umap", pt.size = 0.2)

DimPlot(sc_full, reduction = 'umap', group.by = 'rep_label', pt.size = 0.2)


sc_full_cluster <- subset(x = sc_full, idents = c(1,2,5,8,9,10,11,13,15,16))
DimPlot(sc_full_cluster, reduction = "umap", pt.size = 0.2)

# Next we build a dataframe that will contain the number of cells in each cluster for each rep
rep_label_frame  <- data.frame(sc_full[['rep_label']])
cluster_frame <- data.frame(clusters)
df = data.frame(list(rep_label_frame, cluster_frame))

nb_clusters <- max(as.numeric(clusters))
# aggregate the dataframe so that we have the number of cells for each cluster for each rep
aggreg_df <- aggregate(df$rep_label, list(df$clusters, df$rep_label), FUN=length, drop = FALSE)
# Replace the NA with zeros
aggreg_df[is.na(aggreg_df)] <- 0
abs_repart_cluster = data.frame( 0:(nb_clusters -1), list_reps, check.names = FALSE)
rownames(abs_repart_cluster) <- 0:(nb_clusters -1)
# set all the entries to zero for clarity
abs_repart_cluster[-1] <- 0
# Remove the unwanted col
abs_repart_cluster <- abs_repart_cluster[-1]
# set the right names of the columns
colnames(abs_repart_cluster) <- list_reps


for (col in colnames(abs_repart_cluster)){
  ind_rep_col <- (aggreg_df$Group.2 == col)
  abs_repart_cluster[col] <- aggreg_df$x[ind_rep_col]
}

rel_repart_cluster <- mapply("/", abs_repart_cluster, colSums(abs_repart_cluster))

plot_list_hist <- vector('list', nrow(rel_repart_cluster))
group <- unlist(list_reps)
for (i in 1:nb_clusters){
  rel_repart <- rel_repart_cluster[i,]
  abs_repart <- abs_repart_cluster[i,]
  rel_value <- as.numeric(rel_repart)
  abs_value <- as.numeric(abs_repart)
  rel_df <- data.frame(group,rel_value)
  abs_df <- data.frame(group, abs_value)

  bh <-ggplot(abs_df, aes(x=group, y=abs_value))+
    geom_bar(width = 1, stat = "identity" ) +
    geom_text(aes(label=abs_value), vjust = 1,size = 3) + xlab(paste("cluster", as.character(i-1)))

  plot_list_hist[[i]] <- bh
}



ggarrange(plotlist = plot_list_hist)

# Cluster repartition
ggplot(data.frame(clusters), aes(x=as.numeric(clusters)))  +geom_bar(aes(x = clusters), position = "dodge", stat = "count") +
  stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, binwidth =1)

# Get the markers for each cluster
sc_full.markers <- FindAllMarkers(sc_full, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(sc_full.markers, file=paste("violin_plot", name_rep ,"genes_markers.txt", sep = '/'),sep="\t",quote=FALSE)

dev.off()


###### MONOCLE

setwd("~/")
rdsFile <- "sc_full WT_rep123.rds"
Seurat_obj <- readRDS(file = rdsFile)


nbClusters <- length(summary(Idents(Seurat_obj)))
print('The clusters are: ')
print(summary(Idents(Seurat_obj)))

clusterChosen <- c(1,2,5,8,9,10,11,13,15,16)


pdf(file = paste('first_plots_', rdsFile, paste(as.character(clusterChosen), collapse = "-"), '.pdf', sep = ""))
Seurat_obj[['clusters']] = Idents(Seurat_obj)
Seurat_obj <- subset(Seurat_obj, clusters %in% clusterChosen)

# From Seurat object to Monocle object
Seurat_obj@active.assay = 'RNA'
cds <- SeuratWrappers::as.cell_data_set(Seurat_obj)
cds <- estimate_size_factors(cds)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

cds <- preprocess_cds(cds, num_dim = 100)
## Cluster the cells and compute the partitions. We will change the cluster with Seurat clusters.

cds <- cluster_cells(cds)

## Assign the clusters of the cds to the Seurat clusters
cds@clusters[['UMAP']]$clusters <- Idents(Seurat_obj)

# Learn the graph, takes some time
cds <- learn_graph(cds)

## Plot the cells, according to their  clusters.
plot_cells(cds, color_cells_by = 'cluster', label_groups_by_cluster = TRUE, labels_per_group = 0)


get_earliest_principal_node <- function(cds, time_bin="130-170"){
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex))))]
  root_pr_nodes
}

cds  <- order_cells(cds)

# plot the cells by pseudotime
plot_cells(cds, color_cells_by = 'pseudotime')


nbClusters <- length(summary(Idents(Seurat_obj)))
print(paste('There are ', nbClusters, ' clusters in total.'))
print('The clusters are:')
print(clusterChosen)
nbClustersChosen <- length(summary(Idents(Seurat_obj)))

#### IF YOU WANT TO GET ALL THE CLUSTERS ALREADY CHOSEN BEFORE RUN THIS
subClusterChosen <- clusterChosen
# get only the wanted clusters the cell_data_set. It will be stored in sub_cds
BoolInChosenCluster <- as.vector(cds@clusters[['UMAP']][3][1]$clusters) %in% subClusterChosen
clusters = as.vector(cds@clusters[['UMAP']][3][1]$clusters)
clusterIdOfCellsInChosenCluster <- clusters[clusters %in% subClusterChosen]
cellsInChosenCluster <- row.names(subset(pData(cds),BoolInChosenCluster))
sub_cds <- cds[,cellsInChosenCluster]



## Takes a lot of time
sub_ciliated_cds_pr_test_res <- getCiliatedCds(sub_cds)


qValueChosen = 0.0000000005
geneModuleAndAggMat <- getAggregatedMatrix(sub_cds, sub_ciliated_cds_pr_test_res, clusterIdOfCellsInChosenCluster, qValueChosen)
gene_module_df <- geneModuleAndAggMat[[1]]
agg_mat <- geneModuleAndAggMat[[2]]
showAggregatedMat(agg_mat)
agg_mat
gene_module_df

dev.off()

##If you want to select genes and plot their expression levels in the different cells
pdf(file = paste('first_plots_', rdsFile, paste(as.character(clusterChosen), collapse = "-"), '.pdf', sep = ""))
plot_cells(cds, genes=c("AT2G23770","SC-AT1G78370","AT2G39518","MC2-AT3G46230","MC2-AT3G56880", "AT3G16570", "AT5G42980", "AT2G33380"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)


## If you want to plot mean expression of genes (here Marker_cell) in  different clusters
Marker_cell = c('AT5G25980',
                'AT5G48485',
                'AT1G62480',
                'AT3G16400',
                'AT2G15830',
                'AT1G71050',
                'AT2G19810',
                'AT4G37870',
                'AT5G66400',
                'AT3G23730',
                'AT1G23130',
                'AT1G67865',
                'AT1G64370',
                'AT4G19840',
                'AT2G18328',
                'AT5G18600',
                'AT1G67860',
                'AT1G67870',
                'AT5G45350',
                'AT2G32870',
                'AT2G38540',
                'AT1G66100',
                'AT1G09310',
                'AT3G51600',
                'AT5G25610',
                'AT5G44020',
                'AT3G16370',
                'AT2G27385',
                'AT3G26450',
                'AT1G68530',
                'AT2G10940',
                'AT5G38430',
                'AT3G08940',
                'AT1G72610',
                'AT3G27690',
                'AT2G05070',
                'AT1G12090',
                'AT1G29910',
                'AT2G34420',
                'AT2G34430',
                'AT1G15820',
                'AT2G21330',
                'AT1G06680',
                'AT3G59400',
                'AT2G05100',
                'AT1G67090',
                'AT3G54890',
                'AT5G66570',
                'AT4G38970',
                'AT1G44575',
                'AT1G18740',
                'AT1G74930',
                'AT1G27730',
                'AT2G44840',
                'AT1G80840',
                'AT3G44260',
                'AT5G12030',
                'AT5G12020',
                'AT1G74450',
                'AT4G24570',
                'AT3G16670',
                'AT3G05730',
                'AT3G16660',
                'AT1G56710',
                'AT3G09330',
                'AT1G22900',
                'AT1G08090',
                'AT4G36260',
                'AT4G32950',
                'AT2G43610',
                'AT1G78370',
                'AT3G19710',
                'AT2G30860',
                'AT2G43100',
                'AT5G23010',
                'AT2G46650',
                'AT5G14200',
                'AT4G14040',
                'AT1G16410',
                'AT3G02020'
)

plot_genes_by_group(cds,
                    Marker_cell,
                    group_cells_by="clusters",
                    ordering_type="cluster_row_col",
                    max.size=3)


dev.off()
get_citation(cds)

