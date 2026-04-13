library(glue)
library(data.table)
library(here)

work_directory <- '/raid6/Tianyu/genetic_convergence_SA/'
source(file.path(work_directory, "R", "collect_and_structure_results.R"))
source(file.path(work_directory, "R", "convergence.R"))

batch_name <- "G1_pairwise"
summary_file <- file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", batch_name, "all_results_summary.csv"
)

if (file.exists(summary_file)) {
  results_summary <- fread(summary_file)
} else {
  stop("Summary file not found! Ensure that the results have been generated.")
}

process_lasso_result <- function(result) {
  if (is.null(result)) return(NULL)

  active_genes <- collect_active_features(result)
  active_genes <- paste(active_genes, collapse = "/")

  list(
    p_value = result$p_value,
    active_genes = active_genes
  )
}

all_results <- lapply(results_summary$file_path, load_results)
names(all_results) <- paste(results_summary$treatment1, results_summary$treatment2, sep = "_vs_")

processed_results <- lapply(all_results, process_lasso_result)
processed_results_dt <- rbindlist(processed_results, fill = TRUE, idcol = "comparison")

fwrite(
  processed_results_dt,
  file.path(
    work_directory,
    "yao_2023", "data", "intermediate_data", batch_name, "G2_processed_results.csv"
  )
)

message("Processed results saved.")
