library(here)
library(splatter)
library(SingleCellExperiment)
library(ggplot2)
library(HMC)

work_directory <- file.path(here(), "converge_by_module_real_data", "functions")

# Generate Splatter counts and apply the A/B perturbations.
file.path(work_directory, "10_prepare_splatter_data.R")

# Draw the module heatmap.
file.path(work_directory, "12_draw_heatmap.R") 

# Run XConTest on A versus B in both directions.
file.path(work_directory, "11_run_XConTest.R")

# Draw the heatmap and XConTest convergence plot together.
file.path(work_directory, "13_draw_heatmap_convergence.R") 
