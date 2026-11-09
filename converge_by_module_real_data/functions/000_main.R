library(here)
library(splatter)
library(SingleCellExperiment)
library(ggplot2)
library(HMC)

work_directory <- file.path(here(), "converge_by_module_real_data", "functions")

# Generate Splatter counts and apply the A/B perturbations.
file.path(work_directory, "10_prepare_splatter_pilot.R") |>
  source(chdir = TRUE, local = new.env(parent = globalenv()))

# Draw the module heatmap.
file.path(work_directory, "12_draw_splatter_pilot.R") |>
  source(chdir = TRUE, local = new.env(parent = globalenv()))

# Run XConTest on A versus B in both directions.
file.path(work_directory, "11_test_splatter_pilot.R") |>
  source(chdir = TRUE, local = new.env(parent = globalenv()))

# Draw the heatmap and XConTest convergence plot together.
file.path(work_directory, "13_draw_splatter_pilot_convergence.R") |>
  source(chdir = TRUE, local = new.env(parent = globalenv()))
