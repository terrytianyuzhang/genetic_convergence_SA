# Run this script from the functions directory.
rm(list = ls())

data_path <- file.path('..', 'data', 'raw')
data_file <- file.path(data_path, '10_splatter_pilot_60_120.Rdata')
effect_design_file <- file.path(data_path, '10_effect_design.csv')
dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

pert <- c('A', 'B')
n_samples <- 60
n_control <- 120
n_modules <- 6
module_sz <- 50
n_features <- n_modules * module_sz
n_signal <- 15  # 30% of genes in each modified module.

############################ Simulate baseline raw counts #####################
sce <- splatSimulate(nGenes = n_features,
  batchCells = n_control + length(pert) * n_samples,
  seed = 20261109, sparsify = FALSE, verbose = FALSE)
original_counts <- t(counts(sce))

group_labels <- c(rep(pert, each = n_samples),
                  rep('NT', n_control))

############################ Define sparse module effects #####################
module_fold_change <- list(
  A = c(1.80, 0.65, 1, 1, 1, 1),
  B = c(0.55, 1, 1.60, 1, 1, 1)
)
set.seed(20261110)
signal_positions <- lapply(seq_len(n_modules), function(module)
  sample(seq_len(module_sz), n_signal))

fold_change <- list()
for (group in pert) {
  multiplier <- rep(1, n_features)
  for (module in seq_len(n_modules)) {
    if (module_fold_change[[group]][module] == 1) next
    genes <- (module - 1) * module_sz + signal_positions[[module]]
    multiplier[genes] <- module_fold_change[[group]][module]
  }
  names(multiplier) <- colnames(original_counts)
  fold_change[[group]] <- multiplier
}

############################ Multiply raw counts and round ###################
modified_counts <- original_counts
for (group in pert) {
  rows <- which(group_labels == group)
  changed <- which(fold_change[[group]] != 1)
  modified_counts[rows, changed] <- round(sweep(
    original_counts[rows, changed, drop = FALSE], 2,
    fold_change[[group]][changed], '*'))
}

gene_info <- data.frame(gene = colnames(original_counts),
  module = rep(seq_len(n_modules), each = module_sz),
  position_in_module = rep(seq_len(module_sz), n_modules))
cell_info <- data.frame(cell = rownames(original_counts), group = group_labels)
effect_design <- do.call(rbind, lapply(pert, function(group) data.frame(
  group = group, gene = gene_info$gene, module = gene_info$module,
  fold_change = unname(fold_change[[group]]),
  is_signal = unname(fold_change[[group]]) != 1)))

sim <- list(
  control = modified_counts[group_labels == 'NT', , drop = FALSE],
  treatment = lapply(setNames(pert, pert), function(group)
    modified_counts[group_labels == group, , drop = FALSE]),
  original_counts = original_counts, cell_info = cell_info,
  gene_info = gene_info, fold_change = fold_change,
  module_fold_change = module_fold_change,
  splatter_params = metadata(sce)$Params
)

save(sim, file = data_file, compress = FALSE)
write.csv(effect_design, effect_design_file, row.names = FALSE)
