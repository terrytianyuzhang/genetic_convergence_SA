file_name <- file.path(
  work_directory_yao_2023,
  "data", "module_list_df_2000_genes.csv"
)

original_gene_module <- read_csv(file_name)

#only keep the immunity relevant modules
sensitivity_gene_module <- original_gene_module |> 
  filter(cluster_index %in% c(6, 12, 13, 20, 31))

gene_count <- sensitivity_gene_module |> nrow()

write_csv(
  sensitivity_gene_module,
  file = file.path(
    work_directory_yao_2023,
    "data", "intermediate_data",
    paste0("F1_module_list_", gene_count, "_genes.csv")
  )
)

