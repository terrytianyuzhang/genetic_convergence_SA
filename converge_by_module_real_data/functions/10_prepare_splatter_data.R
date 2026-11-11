# Run this script from the functions directory.
rm(list = ls())

data_path <- file.path('..', 'data', 'raw')
data_file <- file.path(data_path, '10_splatter_data.Rdata')
effect_design_file <- file.path(data_path, '10_effect_design.csv')
dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

pert <- c('A', 'B', 'C', 'D')
n_samples <- 60
n_control <- 120
n_modules <- 20
module_sz <- 50
n_features <- n_modules * module_sz
n_signal <- 10  # 20% of genes in each modified module.

############################ Simulate baseline raw counts #####################
sce <- splatSimulate(nGenes = n_features + 200,
  batchCells = n_control + length(pert) * n_samples,
  seed = 20261109, sparsify = FALSE, verbose = FALSE)
original_counts <- t(counts(sce))

group_labels <- c(rep(pert, each = n_samples),
                  rep('NT', n_control))

# Keep 1000 genes that vary in every control-plus-two-treatment comparison.
pair_groups <- combn(pert, 2, simplify = FALSE)
variable_genes <- vapply(seq_len(ncol(original_counts)), function(gene) {
  all(vapply(pair_groups, function(pair)
    sd(original_counts[group_labels %in% c('NT', pair), gene]) > 0,
    logical(1)))
}, logical(1))
stopifnot(sum(variable_genes) >= n_features)
original_counts <- original_counts[, which(variable_genes)[seq_len(n_features)],
  drop = FALSE]

############################ Define sparse module effects #####################
module_fold_change <- list(
  A = c(1.4, 1.5, 1.45, 1.4, 1.3, rep(1, n_modules - 5)),
  B = c(0.95, 0.5, 0.45, 0.55, 0.95, rep(1, n_modules - 5)),
  C = c(rep(1, 5), 1.4, 1.5, 1.45, 1.4, 1.3, rep(1, n_modules - 10)),
  D = c(0.5, 1.2, 1.5, 0.6, 1.45, rep(1, n_modules - 5))
)
set.seed(20261111)
signal_positions <- lapply(setNames(pert, pert), function(group)
  lapply(seq_len(n_modules), function(module)
    sample(seq_len(module_sz), n_signal)))

fold_change <- list()
for (group in pert) {
  multiplier <- rep(1, n_features)
  for (module in seq_len(n_modules)) {
    if (module_fold_change[[group]][module] == 1) next
    genes <- (module - 1) * module_sz + signal_positions[[group]][[module]]
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
