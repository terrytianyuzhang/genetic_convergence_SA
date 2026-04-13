library(glue)
library(data.table)
library(here)

work_directory <- here::here()
source(file.path(work_directory, "R", "collect_and_structure_results.R"))
source(file.path(work_directory, "R", "convergence.R"))

batch_name <- "G1_pairwise"
one_pair_file <- file.path(work_directory, "yao_2023", "data", "intermediate_data", batch_name, "STAT1_STAT2.rds")

one_pair <- readRDS(one_pair_file)


