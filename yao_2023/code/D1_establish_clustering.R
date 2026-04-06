# This script was created from 31_establish_clustering.R.
work_directory_yao_2023 <- file.path(work_directory, "yao_2023")

file_name <- file.path(
  work_directory_yao_2023,
  "data", "raw_data", "GSM6858447_KO_conventional.rds"
)

# Read the original data and keep non-targeting control cells.
Cleary_raw <- readRDS(file_name)
Cleary_raw@meta.data
Cleary_subset <- Cleary_raw[, Cleary_raw$Guides_collapsed_by_gene %in% "non-targeting"]
remove(Cleary_raw)

###keep high expression level genes
mean_exp <- rowMeans(Cleary_subset@assays$RNA@counts/Cleary_subset$Total_RNA_count)
gene_count <- 2000
genes_selected <- names(sort.int(mean_exp, decreasing = T))[1:gene_count]

CSCORE_result <- CSCORE(Cleary_subset, genes = genes_selected)
saveRDS(CSCORE_result, 
        file = file.path(
          work_directory_yao_2023,
          "data", "intermediate_data", 
          paste0("D1_CSCORE_result_", gene_count, "_gene.rds")
        ))
