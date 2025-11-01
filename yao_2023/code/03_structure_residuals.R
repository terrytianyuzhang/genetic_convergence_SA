   
file_name <- file.path(work_directory, 
                       'yao_2023/data/intermediate_data/after_subsetting.rds')
Cleary <- readRDS(file_name)

residual_list <- readRDS(file = file.path(work_directory, 
                                          "yao_2023/data/intermediate_data/poisson_residual.rds"))

residual_matrix <- data.frame()

for(chunk_index in 1:length(residual_list)){
  if(chunk_index == 1){
    residual_matrix <- residual_list[[chunk_index]]
  }else{
    one_chunk_matrix <- residual_list[[chunk_index]]
    if(identical(residual_matrix$ID, one_chunk_matrix$ID)){
      residual_matrix <- cbind(residual_matrix, one_chunk_matrix[, !colnames(one_chunk_matrix) %in% c("ID")])
    }else{
      warning("the IDs do not match")
    }
  }
}

meta_data <- as.data.table(Cleary@meta.data[,c("Guides_collapsed_by_gene", "Cell_cycle_phase")])
meta_data$ID <- rownames(Cleary@meta.data)

#confirm the batch effects are removed
if(0){
  seurat_obj <- CreateSeuratObject(counts = t(residual_matrix[, !colnames(residual_matrix) == "ID"]), 
                                   project = "MyProject")
  seurat_obj$group <- meta_data$Cell_cycle_phase
  subset_cells <- sample(colnames(seurat_obj), size = 2000)
  subset_obj <- subset(seurat_obj, cells = subset_cells)
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 2)  # Adjust number of PCs as needed
  
  # Plot PCA colored by cell groups
  DimPlot(seurat_obj, reduction = "pca", group.by = "group") + 
    ggtitle("PCA Plot Colored by Cell Groups")
  ggsave(filename = file.path(work_directory, 
                              'yao_2023/report/PCA_cell_cycle_after.pdf'))
}
#####

residual_matrix <- merge(meta_data, residual_matrix, by = 'ID')
saveRDS(residual_matrix, file = file.path(work_directory, 
                                          'yao_2023/data/intermediate_data/residual_matrix.rds'))
