work_directory_yao_2023 <- paste0(work_directory, '/yao_2023')

file_name <- file.path(work_directory_yao_2023, 'data/GSM6858447_KO_conventional.rds')

###read in the original data
Cleary_raw <- readRDS(file_name)
Cleary_raw@meta.data
Cleary_subset <- Cleary_raw[, Cleary_raw$Guides_collapsed_by_gene %in% c('non-targeting')]
remove(Cleary_raw)

###keep high expression level genes
mean_exp <- rowMeans(Cleary_subset@assays$RNA@counts/Cleary_subset$Total_RNA_count)
gene_count <- 2000
genes_selected <- names(sort.int(mean_exp, decreasing = T))[1:gene_count]

library(CSCORE)
CSCORE_result <- CSCORE(Cleary_subset, genes = genes_selected)

saveRDS(CSCORE_result, file = file.path(work_directory_yao_2023, 'data/CSCORE_result_', gene_count, '_gene.rds'))
CSCORE_result <- readRDS(file = file.path(work_directory_yao_2023, 'data/CSCORE_result_', gene_count, '_gene.rds'))
CSCORE_coexp <- CSCORE_result$est

CSCORE_p <- CSCORE_result$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp[p_matrix_BH > 0.05] <- 0

# Compute the adjacency matrix based on the co-expression matrix
adj = WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 10)
# Compute the topological overlap matrix
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
# Run hierarchical clustering as in the WGCNA workflow
hclust_dist = hclust(as.dist(dissTOM), method = "average") 
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                     distM = dissTOM, 
                                     deepSplit = 2,
                                     pamRespectsDendro = FALSE,
                                     minClusterSize = 10)
# For more instructions on how to tune the parameters in the WGCNA workflow,
# please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

names(memb) = genes_selected
memb_tab <- table(memb)
module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))
module_list
saveRDS(module_list, file = file.path(work_directory_yao_2023, 'data/module_list_', gene_count,'_genes.rds'))

library(clusterProfiler)
library(org.Hs.eg.db)
library(parallel)
universe <- genes_selected
ego_result <- mclapply(1:length(module_list), function(i){
  enrichGO(gene = module_list[[i]],
           OrgDb = 'org.Hs.eg.db', # human
           keyType = "SYMBOL",
           ont = "ALL",
           pAdjustMethod = "BH",
           universe = universe,
           pvalueCutoff = 0.5)
}, mc.cores = 3, mc.cleanup = TRUE)
saveRDS(ego_result, file = file.path(work_directory_yao_2023, 'data/ego_result_', gene_count,'_genes.rds'))

###based on the GO term concentration, we further split the cluster 1 in module_list

gray_gene_num <- length(module_list[[1]])
GO_correlation <- matrix(0, gray_gene_num, gray_gene_num)
rownames(GO_correlation) <- module_list[[1]]
colnames(GO_correlation) <- module_list[[1]]

gray_cluster <- data.table(all_ego[all_ego$cluster_index == 1,])
for(i in 1:nrow(gray_cluster)){
  # gray_cluster[i, 'gene_number'] <- length(strsplit(as.character(gray_cluster[i, ]$geneID), "/")[[1]])
  name_overlap <- strsplit(as.character(gray_cluster[i, ]$geneID), "/")[[1]]
  GO_correlation[name_overlap, name_overlap] <- GO_correlation[name_overlap, name_overlap] + 1
  
}
GO_correlation <- 2*exp(GO_correlation)/(1+exp(GO_correlation)) - 1
diag(GO_correlation) <- 1

# Compute the adjacency matrix based on the co-expression matrix
adj = WGCNA::adjacency.fromSimilarity(abs(GO_correlation), power = 10)
# Compute the topological overlap matrix
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- module_list[[1]]
# Run hierarchical clustering as in the WGCNA workflow
hclust_dist = hclust(as.dist(dissTOM), method = "average") 
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                     distM = dissTOM, 
                                     deepSplit = 2,
                                     pamRespectsDendro = FALSE,
                                     minClusterSize = 10)
# For more instructions on how to tune the parameters in the WGCNA workflow,
# please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
names(memb) = module_list[[1]]
memb_tab <- table(memb)
GO_module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))
GO_module_list

finer_module_list <- append(GO_module_list, module_list[-1])
saveRDS(finer_module_list, file = file.path(work_directory_yao_2023, 'data/finer_module_list_', gene_count,'_genes.rds'))

finer_module_list <- readRDS(file = file.path(work_directory_yao_2023, 'data/finer_module_list_', gene_count,'_genes.rds'))
module_list_df <- data.frame()
for(module_index in 1:length(finer_module_list)){
  module_list_df <- rbind(module_list_df, 
                          data.frame(gene_name = finer_module_list[[module_index]],
                                     cluster_index = module_index))
}
# cluster1 <- sum(module_list_df$cluster_index == 1)
# module_list_df[module_list_df$cluster_index == 1, ]$cluster_index <- (length(module_list) + 1): (length(module_list) + cluster1)
# module_list_df$cluster_index <- module_list_df$cluster_index - 1
# module_list_df <- module_list_df[order(module_list_df$cluster_index), ] 
write.csv(module_list_df, file.path(work_directory_yao_2023, 'data/module_list_df_', gene_count,'_genes.csv'), row.names = FALSE)






universe <- module_list[[1]]
ego_result <- mclapply(1:length(GO_module_list), function(i){
  enrichGO(gene = GO_module_list[[i]],
           OrgDb = 'org.Hs.eg.db', # human
           keyType = "SYMBOL",
           ont = "ALL",
           pAdjustMethod = "BH",
           universe = universe,
           pvalueCutoff = 0.5)
}, mc.cores = 3, mc.cleanup = TRUE)
saveRDS(ego_result, file = file.path(work_directory_yao_2023, 'data/ego_result_', gene_count,'_genes_second.rds'))

all_ego <- data.frame()
for(cluster_index in 1:length(ego_result)){
  one_ego <- ego_result[[cluster_index]]
  if(is.null(one_ego)) next
  one_ego <- one_ego[one_ego$p.adjust < 0.05, c('Description', 'GeneRatio', 'p.adjust', 'geneID')]
  
  if(nrow(one_ego) == 0) next
  
  one_ego$cluster_index <- cluster_index
  all_ego <- rbind(all_ego, one_ego)
}

for(i in 15:20){
  gene_names1 <- finer_module_list[[i]]
  gene_names2 <- finer_module_list[[25]]
  gene_names1 <- sample(gene_names1, size = min(50, length(gene_names1)))
  gene_names2 <- sample(gene_names2, size = min(50, length(gene_names2)))
  gene_names <- union(gene_names1, gene_names2)
  
  mat <- CSCORE_coexp[gene_names,gene_names]
  corrplot(mat) 
}

gene_string <- as.character(gray_cluster[2,'geneID'])
gene_count <- length(strsplit(gene_string, "/")[[1]])
print(gene_count)

sink(file.path(work_directory_yao_2023, 'data/top_ego_result_', gene_count,'_genes.txt'))
for(cluster_index in 1:length(ego_result)){
  print(cluster_index)
  cat('\n')
  
  print(all_ego[all_ego$cluster_index == cluster_index, c('Description', 'GeneRatio', 'p.adjust', 'geneID')])
  cat('\n')
}

sink()

# 
# error_index <- which(sapply(ego_result, is.null))
# ego_result_clean <- ego_result[-error_index]
# top_enrich_clusters <- which(sapply(ego_result_clean, function(x) 
#   (x@result$p.adjust[1] < 0.2) ))
# top_enrich_go <- lapply(top_enrich_clusters, function(i) ego_result_clean[[i]]@result[1:10,])

sink(file.path(work_directory_yao_2023, 'data/top_ego_result_', gene_count,'_genes.txt'))
index_to_fetch <- 1
for(i in 1:length(ego_result)){
  print(i)
  cat('\n')
  
  if(i == error_index){
    print('this cluster occured error, likely all genes are missing from the data base')
    cat('\n')
  }
  else{
    print(top_enrich_go[[index_to_fetch]][, c('Description', 'GeneRatio', 'p.adjust', 'geneID')])
    cat('\n')
    index_to_fetch <- index_to_fetch + 1
  }
}
sink()


library(corrplot)
for(i in 1:10){
  gene_names1 <- module_list[[i]]
  gene_names2 <- module_list[[i+1]]
  gene_names1 <- sample(gene_names1, size = min(50, length(gene_names1)))
  gene_names2 <- sample(gene_names2, size = min(50, length(gene_names2)))
  gene_names <- union(gene_names1, gene_names2)
  
  mat <- CSCORE_coexp[gene_names,gene_names]
  corrplot(mat) 
}
print(length(module_list))


ego_result <- readRDS(file = file.path(work_directory_yao_2023, 
                                       'data/ego_result_', gene_count,'_genes.rds'))
all_ego <- data.frame()
for(cluster_index in 1:length(ego_result)){
  one_ego <- ego_result[[cluster_index]]
  if(is.null(one_ego)) next
  one_ego <- one_ego[one_ego$p.adjust < 0.05, c('Description', 'GeneRatio', 'p.adjust', 'geneID')]
  
  if(nrow(one_ego) == 0) next
  
  one_ego$cluster_index <- cluster_index
  all_ego <- rbind(all_ego, one_ego)
}
