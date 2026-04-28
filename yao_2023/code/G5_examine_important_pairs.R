batch_name <- "G1_pairwise"
processed_results <- fread(file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", batch_name, "G2_processed_results.csv"
))

processed_results <- processed_results[, .(comparison, p_value, active_genes)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(
  processed_results,
  processed_results,
  by.x = c("gene1", "gene2"),
  by.y = c("gene2", "gene1")
)
merged_results[, p_max := pmax(p_value.x, p_value.y)]

threshold <- 0.025
merged_results[, is_significant := 0]
merged_results[, is_significant := is_significant + (p_max < threshold)]
merged_results[, is_significant := is_significant + (p_max < threshold / 2450)]
merged_results[, is_significant := factor(is_significant, levels = c(0, 1, 2))]

merged_results |> filter(gene1 %in% c("TRIB1", "TRAF6", "HSP90B1") & gene2 %in% c("TRIB1", "TRAF6", "HSP90B1"))
merged_results |> filter(gene1 %in% c("STAT1", "STAT2", "TYK2") & gene2 %in% c("STAT1", "STAT2", "TYK2"))

merged_results |> filter(gene1 %in% c("MYD88", "IRAK1", "IRAK4") & 
                         gene2 %in% c("MYD88", "IRAK1", "IRAK4"))

merged_results |> filter(gene1 %in% c("IRAK1", "IRAK4", "IKBKB", "IKBKG", "RELA") & 
                           gene2 %in% c("IRAK1", "IRAK4", "IKBKB", "IKBKG", "RELA")) |>
                  filter(is_significant == 2)

