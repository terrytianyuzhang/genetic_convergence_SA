library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(glue)

work_directory_yao_2023 <- paste0(work_directory, '/yao_2023')

clustering_file <- file.path(work_directory_yao_2023, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))

cluster_ids <- unique(clustering$cluster_index)
universe <- clustering$gene_name
# Set up parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

dir.create(file.path(work_directory_yao_2023, "/data/final_data/module_GO"), showWarnings = FALSE)

foreach(cluster_id = cluster_ids, .packages = c("data.table", "clusterProfiler", "org.Hs.eg.db", "glue")) %dopar% {
  
  print(cluster_id)
  # Subset gene list for this cluster
  gene_list <- clustering[cluster_index == cluster_id, gene_name]
  
  ego <- enrichGO(gene = gene_list,
           OrgDb = 'org.Hs.eg.db', # human
           keyType = "SYMBOL",
           ont = "ALL",
           pAdjustMethod = "BH",
           universe = universe,
           pvalueCutoff = 0.2)
  
  # Save result
  if (!is.null(ego) && nrow(ego) > 0) {
    saveRDS(ego, file = file.path(work_directory_yao_2023, "data/final_data/module_GO/module_", cluster_id,".rds"))
    write.csv(as.data.frame(ego), 
              file = file.path(work_directory_yao_2023, "data/final_data/module_GO/module_", cluster_id,".csv"), row.names = FALSE)
  }
}

# Stop the cluster
stopCluster(cl)