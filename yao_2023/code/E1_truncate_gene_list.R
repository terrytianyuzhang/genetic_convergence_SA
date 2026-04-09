file_name <- file.path(
  work_directory_yao_2023,
  "data", "module_list_df_2000_genes.csv"
)

original_gene_module <- read_csv(file_name)

#modules 24 - 42 are based on correlation, the first 23 are based on known GO function.
sensitivity_gene_module <- original_gene_module |> filter(cluster_index %in% 1:23)

gene_count <- sensitivity_gene_module |> nrow()

write_csv(
  sensitivity_gene_module,
  file = file.path(
    work_directory_yao_2023,
    "data", "intermediate_data",
    paste0("E1_module_list_", gene_count, "_genes.csv")
  )
)

