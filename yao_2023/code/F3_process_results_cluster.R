library(glue)
library(data.table)
library(here)

work_directory <- here::here()
source(file.path(work_directory, "R", "collect_and_structure_results.R"))
source(file.path(work_directory, "R", "convergence.R"))

batch_name <- "F2_pairwise"
summary_file <- file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", batch_name, "all_results_summary.csv"
)

clustering_file <- file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", "F1_module_list_201_genes.csv"
)
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

if (file.exists(summary_file)) {
  results_summary <- fread(summary_file)
} else {
  stop("Summary file not found! Ensure that the results have been generated.")
}

all_results <- lapply(results_summary$file_path, load_results)
names(all_results) <- paste(results_summary$treatment1, results_summary$treatment2, sep = "_vs_")

processed_results <- lapply(all_results, process_result, group = clustering$cluster_index)
processed_results_dt <- rbindlist(processed_results, fill = TRUE, idcol = "comparison")

fwrite(
  processed_results_dt,
  file.path(
    work_directory,
    "yao_2023", "data", "intermediate_data", batch_name, "F3_processed_results.csv"
  )
)

message("Processed results saved.")
