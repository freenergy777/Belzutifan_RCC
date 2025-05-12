suppressPackageStartupMessages({
  library(R.utils)
  library(openxlsx)
  library(Seurat)
  library(ggtext)
  library(SeuratDisk)
  library(reticulate)
  library(ggpubr)
  library(SeuratWrappers)
  library(monocle3)
  library(dplyr)
  library(pheatmap)
})

GSE269823_total <- readRDS('GSE269826_Cultured.rds')
data_seurat_umap_coords <- as.data.frame(GSE269823_total@reductions[["umap"]]@cell.embeddings) 
keep_cells <- rownames(data_seurat_umap_coords)[
  data_seurat_umap_coords[,1] <= 10 &
    data_seurat_umap_coords[,2] <= 0]

GSE269823_total <- subset(GSE269823_total, cells = keep_cells)
GSE269823_total_meta <- GSE269823_total@meta.data
write.xlsx(GSE269823_total_meta, 'D:/2025.04/belzutifan/GSE269823_cultured_metadata.xlsx')
write.xlsx(GSE269823_17_meta, 'D:/2025.04/belzutifan/GSE269823_T17_metadata.xlsx')


DimPlot(GSE269823_total, reduction = "umap", label =T, label.size = 5,  repel = T, pt.size = 2,
        group.by = "cell_type11") + ggtitle('Clusters')
DimPlot(GSE269823_total, reduction = "tsne", label =T, label.size = 5,  repel = T, pt.size = 2,
        group.by = "cell_type11") + ggtitle('Clusters')


seurat_subset <- subset(GSE269823_total, subset = cell_type11 %in% c('ccRCC_DMSO', 'ccRCC_PT2385'))

# edit gene expression data
seurat_subset = seurat_subset[!grepl("^RP", seurat_subset@assays[["RNA"]]@data@Dimnames[[1]]),]
cds <- as.cell_data_set(seurat_subset, assay = "RNA") # using SeuratWrappers
cds_RCC_epi <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3) # resolution = 1e-3 is reasonable
cds_RCC_epi <- learn_graph(cds_RCC_epi)
plot_cells(cds_RCC_epi,
           color_cells_by = "cell_type11", trajectory_graph_segment_size = 2,
           label_cell_groups=FALSE,
           label_leaves=FALSE, cell_size = 0.5,
           label_branch_points=FALSE,
           graph_label_size=3)

root_candidates <- which(
  cds_RCC_epi@colData@listData[["cell_type11"]] == "ccRCC_DMSO" 
)

colData(cds_RCC_epi)[['root_or_not']] <- FALSE  
colData(cds_RCC_epi)[['root_or_not']][root_candidates] <- TRUE 

closest_vertex <- cds_RCC_epi@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex 
closest_vertex <- as.matrix(closest_vertex[colnames(cds_RCC_epi),])
root_pr_nodes <- igraph::V(principal_graph(cds_RCC_epi)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[root_candidates,]))))]
cds_RCC_epi <- order_cells(cds_RCC_epi, root_pr_nodes=root_pr_nodes) # designate starting principal point where pseudo-time begin

cds_RCC_epi@principal_graph_aux[["UMAP"]]$pseudotime <- pseudotime(cds_RCC_epi)
colData(cds_RCC_epi)$pseudotime <- cds_RCC_epi@principal_graph_aux[["UMAP"]]$pseudotime

plot_cells(cds_RCC_epi,
           color_cells_by = "pseudotime", trajectory_graph_segment_size = .7,
           label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
           cell_size = 0.8, alpha = 0.5, scale_to_range = T, trajectory_graph_color = "black", label_principal_points = F)

# choose sub trajectory manually
root_nodes <- function(cds_RCC_epi, reduction_method="UMAP"){
  g = principal_graph(cds_RCC_epi)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                           cds_RCC_epi@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds_RCC_epi@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}
branch_nodes <- function(cds_RCC_epi,reduction_method="UMAP"){
  g = principal_graph(cds_RCC_epi)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds_RCC_epi, reduction_method) == FALSE]
  return(branch_points)
}
leaf_nodes <- function(cds_RCC_epi,reduction_method="UMAP"){
  g = principal_graph(cds_RCC_epi)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds_RCC_epi, reduction_method) == FALSE]
  return(leaves)
}
dp_mst <- principal_graph(cds_RCC_epi)[["UMAP"]]
mst_root_nodes <- root_nodes(cds_RCC_epi, "UMAP")
mst_branch_nodes <- branch_nodes(cds_RCC_epi, "UMAP")
mst_leaf_nodes <- leaf_nodes(cds_RCC_epi, 'UMAP')

# designate Y_59 as leaf node
ENDnodes <- c("Y_22")
path.nodes <- NULL
for (end in ENDnodes){
  path <- igraph::shortest_paths(
    dp_mst, from <- root_pr_nodes, to = end, mode = "all",
    algorithm = "unweighted")
  nodes <- path$vpath[[1]]
  path.nodes <- c(path.nodes, nodes)
}
path.nodes <- unique(path.nodes)
cells.branch <- closest_vertex[closest_vertex %in% path.nodes,] # select cell that included in path.node
cds_RCC_epi@colData$sub_trajectory <- names(cds_RCC_epi@clusters$UMAP$clusters) %in% names(cells.branch) # whether the cell 
cds_RCC_epi@colData$optimal_path <- rownames(cds_RCC_epi@colData) %in% names(cells.branch)
plot_cells(cds_RCC_epi, color_cells_by = "optimal_path")
plot_cells(cds_RCC_epi, color_cells_by = "sub_trajectory")
table(cds@colData$sub_trajectory, cds@clusters$UMAP$partitions) 
table(cds@clusters$UMAP$partitions) 
table(cds@colData$sub_trajectory) 

# differential expression gene across a single-cell trajectory 
cds_sub <- cds_RCC_epi[, cds_RCC_epi@colData$sub_trajectory]
table(cds_sub@colData@listData[["cell_type11"]])
cds_sub_principal_graph_segment <- as.data.frame(cds_RCC_epi@principal_graph_aux@listData[["UMAP"]][["R"]])
write.xlsx(cds_sub_principal_graph_segment, 'D:/2025.04/monocle3_practice/20250415_3k/20250417_cds_sub_principal_graph_segment.xlsx',
           rowNames=T)

deg_res <- graph_test(cds_sub, neighbor_graph="principal_graph", cores = 4) # identify genes with interesting patterns of expression that fall only within the region of the trajectory
write.xlsx(deg_res, 'D:/2025.05/belzutifan/20250505_Y_1_trajectory.xlsx', rowNames=T)

pr_deg_ids <- rownames(subset(deg_res, morans_I > 0.1))
pr_deg_ids_1 <- rownames(subset(deg_res, morans_I > 0))
cds_sub <- preprocess_cds(cds_sub, num_dim = 10) # PCA re 
cds_sub <- reduce_dimension(cds_sub) # UMAP re 

genes_to_plot <- rownames(deg_res)[order(deg_res$morans_I, decreasing = TRUE)][1:70]
exprs_mat <- logcounts(cds_sub)[genes_to_plot, ]
summary(exprs_mat@x)
clustering <- kmeans(exprs_mat, centers = 7)
exprs_ordered <- exprs_mat[, order(cds_sub@colData$pseudotime, na.last = NA)] # make pseudotime order
gene_order <- names(sort(clustering$cluster))  # gene names
exprs_final <- exprs_ordered[gene_order, ]

label_col <- rep("", ncol(exprs_final))
label_col[which(colnames(exprs_final) == "5739STDY8351218_GCCAAATGTCTTGTCC-1")] <- "RCC start"

pheatmap(exprs_final,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         labels_col = label_col,
         annotation_row = data.frame(Cluster = factor(clustering$cluster[gene_order])))

pheatmap(exprs_final,
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         labels_col = label_col,
         annotation_row = data.frame(Cluster = factor(clustering$cluster[gene_order])))

cds_sub_df <- cds_sub_df[cds_sub_df$broad_type %in% c('Epi_PT', 'RCC'),]

table(cds_sub_df$broad_type)
cds_sub_metadata <- as.data.frame(cds_sub@colData)
write.xlsx(cds_sub_metadata, 'D:/2025.04/monocle3_practice/20250421_3k_1e-3/20250421_cds_sub_metadata_3k.xlsx')

# make cell cluster equal number for pySCENIC
data <- data[, data$broad_type %in% c('RCC', 'Epi_PT')]
table(data$broad_type)
data_equal = do.call(cbind, lapply(unique(data$broad_type), function(x) {
  subset = data[, data$broad_type == x]
  if (ncol(subset) > 521) {
    set.seed(1)
    subset = subset[, sample(1:ncol(subset), 521)]
  }
  return(subset)
}))
table(data_equal$broad_type)

#
data_equal_count <- cds_sub@assays@data@listData[["counts"]]
data_equal_count <- as.data.frame(data_equal_count)
data_equal_count_edited <- data_equal_count[rownames(data_equal_count) %in% pr_deg_ids_1,]

write.csv(data_equal_count_edited, file = "D:/2025.05/belzutifan/20250507_data_count.csv", quote = FALSE)
write.csv(trajectory_count, file = "D:/2025.04/monocle3_practice/20250421_3k_1e-3/20250422_trajectory_count.csv", quote = FALSE)

#
cds_sub_df_1 <- data.frame(
  pseudotime = cds_sub@colData@listData[["pseudotime"]],
  cell_type = cds_sub@colData@listData[["cell_type11"]],
  row.names = colnames(cds_sub)
)
write.csv(cds_sub_df, 'D:/2025.05/monocle3_practice/resolution_1e-3/20250507_order.csv')

# barplot
cds_sub_df <- data.frame(
  pseudotime = cds_sub@colData@listData[["pseudotime"]],
  cell_type = cds_sub@colData@listData[["cell_type11"]],
  cycle = cds_sub@colData@listData[["Phase"]],
  HIF2a = cds_sub@colData@listData[["HIF2_specific_1"]],
  HIF1a = cds_sub@colData@listData[["HIF1_specific_1"]],
  row.names = colnames(cds_sub)
)
seurat_subset_df$HIF2a_bin <- cut(seurat_subset_df$HIF2a, breaks = 2, labels = c('Low','High'))
seurat_subset_df$HIF2a_bin <- cut(seurat_subset_df$HIF2a, breaks = 3, labels = c('Low','Medium','High'))
seurat_subset_df$HIF1a_bin <- cut(seurat_subset_df$HIF1a, breaks = 3, labels = c('Low','Medium','High'))
cds_sub_df$pseudotime_bin <- cut(cds_sub_df$pseudotime, breaks = 5, labels = c('1','2','3','4','5'))

plot_df <- cds_sub_df %>%
  group_by(pseudotime_bin, cycle) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(pseudotime_bin) %>%
  mutate(freq = n / sum(n))

plot_df_1 <- cds_sub_df %>%
  group_by(pseudotime_bin, cell_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(pseudotime_bin) %>%
  mutate(freq = n / sum(n))

plot_df_2 <- cds_sub_df %>%
  group_by(pseudotime_bin, HIF1a) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(pseudotime_bin) %>%
  mutate(freq = n / sum(n))

ggplot(plot_df, aes(x = pseudotime_bin, y = freq, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_brewer(palette = "PuBu") +
  labs(x = "pseudo-time", y = "Proportion", fill = "Cell Cycle") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = 'bold'))
