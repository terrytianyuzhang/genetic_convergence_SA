
file_name <- file.path(work_directory, 
                       'yao_2023/data/raw_data/after_subsetting.rds')

###remove low quanlity cells
Cleary <- readRDS(file_name)

subset_cells <- sample(colnames(Cleary), size = 500)  # Adjust the sample size as needed
Cleary_sub <- subset(Cleary, cells = subset_cells)
VlnPlot(Cleary_sub, features = c("nFeature_RNA", "nCount_RNA", "Percent_mitochondrial_reads"), ncol = 3)
ggsave(file = file.path(work_directory, 
                        'yao_2023/report/violin_plot_QC.pdf'))

Cleary <- subset(Cleary, 
                 subset = nFeature_RNA > 500 & nCount_RNA < 2e4 & Percent_mitochondrial_reads < 10)
dim(Cleary@assays$RNA@counts)

####show batch effect
if(0){
  subset_cells <- sample(colnames(Cleary), size = 500)  # Adjust the sample size as needed
  Cleary_sub <- subset(Cleary, cells = subset_cells)
  Cleary_sub <- NormalizeData(Cleary_sub)
  Cleary_sub <- FindVariableFeatures(Cleary_sub, selection.method = "vst", nfeatures = 2000)
  Cleary_sub <- ScaleData(Cleary_sub)
  
  # Perform PCA
  Cleary_sub <- RunPCA(Cleary_sub)
  pca_plot <- DimPlot(Cleary_sub, reduction = "pca", group.by = "Cell_cycle_phase") + ggtitle("PCA colored by cell cycle")
  print(pca_plot)
  ggsave(file = file.path(work_directory, 
                          'yao_2023/report/PCA_cell_cycle_before.pdf'))
  pca_plot <- DimPlot(Cleary_sub, reduction = "pca", group.by = "10X_channel") + ggtitle("PCA colored by cell cycle")
  print(pca_plot)
  ggsave(file = file.path(work_directory, 
                          'yao_2023/report/PCA_chennel_before.pdf'))
}
#
count_matrix <- Cleary@assays$RNA@counts
# nCount_RNA_subset <- colSums(count_matrix)
meta_data <- as.data.table(Cleary@meta.data[,c("Cell_cycle_phase", "nCount_RNA")])
meta_data$ID <- colnames(count_matrix)
# meta_data$nCount_RNA_subset <- nCount_RNA_subset


chunk_size <- 100
residual_list <- list()
for(chunk_index in 1:ceiling(nrow(count_matrix)/chunk_size)){
  print(glue("chunk_index is ", chunk_index))
  
  working_df <- count_matrix[((chunk_index-1)*chunk_size + 1) : (chunk_index*chunk_size),]
  working_df <- as.data.table(t(working_df))
  working_df[, ID:= colnames(count_matrix)]
  
  working_df <- merge(working_df, meta_data, by = "ID")
  working_df[, nCount_RNA:= as.numeric(nCount_RNA)]
  working_df[, Cell_cycle_phase:= as.factor(Cell_cycle_phase)]
  working_df <- data.frame(working_df)
  
  adjusted_df <- working_df
  # ggplot()+
  #   geom_histogram(aes(CD63, y = ..density.., group = Cell_cycle_phase, fill = Cell_cycle_phase),
  #                  position = "identity", alpha = 0.5, data = working_df)
  
  for(gene_index in 1:ncol(working_df[, !colnames(working_df) %in% c("ID", "nCount_RNA", "Cell_cycle_phase")])){
    print(glue("gene_index is ", gene_index))
    model<-glm(working_df[, gene_index + 1]~Cell_cycle_phase+offset(log(nCount_RNA)), 
               family=poisson, data=working_df)
    coef_G2M<-coef(model)["Cell_cycle_phaseG2M"]
    coef_S<-coef(model)["Cell_cycle_phaseS"]    
    
    adjusting_coef <- coef_G2M*(working_df$Cell_cycle_phase == "G2M") + coef_S*(working_df$Cell_cycle_phase == "S")
    residual <- working_df[, gene_index + 1]*exp(-adjusting_coef)/working_df$nCount_RNA
    adjusted_df[, gene_index + 1] <- residual
  }
  residual_list[[chunk_index]] <- adjusted_df[, !colnames(adjusted_df) %in% c("nCount_RNA", "Cell_cycle_phase")]
}
saveRDS(residual_list, file = file.path(work_directory, 
                                   "yao_2023/data/intermediate_data/poisson_residual.rds"))

