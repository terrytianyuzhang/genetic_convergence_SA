library(glue)
library(data.table)
library(grpreg)
library(HMC)
library(doParallel)
library(foreach)
library(here)

work_directory <- here::here()
source(file.path(work_directory, "R", "convergence.R"))

residual_subset <- readRDS(file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", "residual_matrix_all_in_paper.rds"
))

residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

clustering_file <- file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", "F1_module_list_201_genes.csv"
)
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

set.seed(1)
control_subsetting_sample_size <- 500
control <- residual_subset[
  Guides_collapsed_by_gene == "non-targeting",
  -"Guides_collapsed_by_gene"
]
control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]

output_dir <- file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", "F2_pairwise"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(work_directory, "yao_2023", "log"),
  recursive = TRUE,
  showWarnings = FALSE
)

preprocess_one_setting <- function(
  residual_subset,
  treatment1_name,
  treatment2_name,
  clustering
) {
  treatment1 <- residual_subset[
    Guides_collapsed_by_gene == treatment1_name,
    -"Guides_collapsed_by_gene"
  ]
  treatment2 <- residual_subset[
    Guides_collapsed_by_gene == treatment2_name,
    -"Guides_collapsed_by_gene"
  ]

  treatment1 <- treatment1[, match(clustering$gene_name, colnames(treatment1)), with = FALSE]
  treatment2 <- treatment2[, match(clustering$gene_name, colnames(treatment2)), with = FALSE]

  list(treatment1 = treatment1, treatment2 = treatment2)
}

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != "non-targeting"]

log_file <- file.path(work_directory, "yao_2023", "log", "F2_logs.txt")

treatment_pairs <- expand.grid(
  treatment1 = all_treatment_names,
  treatment2 = all_treatment_names
)
treatment_pairs <- treatment_pairs[treatment_pairs$treatment1 != treatment_pairs$treatment2, ]

summary_results <- foreach(
  i = 1:nrow(treatment_pairs),
  .combine = rbind,
  .packages = c("data.table", "glue", "grpreg", "HMC")
) %dopar% {
  treatment1_name <- treatment_pairs$treatment1[i]
  treatment2_name <- treatment_pairs$treatment2[i]

  log_message <- glue(
    "{Sys.time()} - Processing {treatment1_name} vs {treatment2_name} in parallel\n"
  )
  write(log_message, file = log_file, append = TRUE)

  process_data <- preprocess_one_setting(
    residual_subset,
    treatment1_name,
    treatment2_name,
    clustering
  )
  test_result <- convergence_testing(
    control,
    process_data$treatment1,
    process_data$treatment2,
    pca_method = "dense_pca",
    classifier_method = "group_lasso",
    lambda_type = "lambda.min",
    n_folds = 5,
    group = clustering$cluster_index,
    standardize_feature = FALSE,
    verbose = TRUE
  )

  output_filename <- file.path(
    output_dir,
    paste0(treatment1_name, "_", treatment2_name, ".rds")
  )
  saveRDS(test_result, file = output_filename)

  data.table(
    treatment1 = treatment1_name,
    treatment2 = treatment2_name,
    file_path = output_filename
  )
}

fwrite(summary_results, file.path(output_dir, "all_results_summary.csv"))
stopCluster(cl)

message("Processing completed. Summary saved.")
