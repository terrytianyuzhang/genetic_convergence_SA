# This script was split from D1_establish_clustering.R.
# It loads gene-selection metadata and the CSCORE result created by
# D1_establish_clustering.R.

gene_selection <- readRDS(
  file = file.path(
    work_directory_yao_2023,
    "data", "intermediate_data",
    "D1_gene_selection_2000_gene.rds"
  )
)

gene_count <- gene_selection$gene_count
genes_selected <- gene_selection$genes_selected

CSCORE_result <- readRDS(
  file = file.path(
    work_directory_yao_2023,
    "data", "intermediate_data",
    paste0("D1_CSCORE_result_", gene_count, "_gene.rds")
  )
)

CSCORE_coexp <- CSCORE_result$est

CSCORE_p <- CSCORE_result$p_value
p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(
  CSCORE_p[upper.tri(CSCORE_p)],
  method = "BH"
)
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

# Set co-expression entries with BH-adjusted p-values greater than 0.05 to 0.
CSCORE_coexp[p_matrix_BH > 0.05] <- 0

# Compute the adjacency matrix based on the co-expression matrix.
adj <- WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 10)

# Compute the topological overlap matrix.
TOM <- WGCNA::TOMsimilarity(adj)
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected

# Run hierarchical clustering as in the WGCNA workflow.
hclust_dist <- hclust(as.dist(dissTOM), method = "average")
memb <- dynamicTreeCut::cutreeDynamic(
  dendro = hclust_dist,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 10
)

# For more instructions on how to tune the parameters in the WGCNA workflow,
# please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

names(memb) <- genes_selected
memb_tab <- table(memb)
module_list <- lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))
module_list

module_list_df <- tibble()
# I skipped the module 1, let's see if the result is consistent still.
for(module_index in 2:length(module_list)){
  module_list_df <- rbind(module_list_df,
                          tibble(gene_name = module_list[[module_index]],
                                     cluster_index = module_index))
}

gene_count <- nrow(module_list_df)

saveRDS(
  module_list_df,
  file = file.path(
    work_directory_yao_2023,
    "data", "intermediate_data",
    paste0("D2_module_list_", gene_count, "_genes.csv")
  )
)
