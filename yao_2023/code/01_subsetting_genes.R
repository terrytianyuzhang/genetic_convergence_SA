file_name <- file.path(work_directory, 
                       'yao_2023/data/raw_data/GSM6858447_KO_conventional.rds')

###read in the original data
Cleary_raw <- readRDS(file_name)
Cleary_raw@meta.data
Cleary_subset <- Cleary_raw[, Cleary_raw$Guides_collapsed_by_gene %in% c('non-targeting')] #not subsetting non-targeting cell will give similar results

###keep high expression level genes
mean_exp <- rowMeans(Cleary_subset@assays$RNA@counts/Cleary_subset$Total_RNA_count)
gene_count <- 2000
genes_selected <- names(sort.int(mean_exp, decreasing = T))[1:gene_count]

Cleary_subset_gene <- subset(Cleary_raw, features = genes_selected)
saveRDS(Cleary_subset_gene, 
        file = file.path(work_directory, 
                         'yao_2023/data/raw_data/after_subsetting.rds'))
