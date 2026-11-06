# 01: Read Xin controls, create A-F/NT, and draw gene/module effect heatmaps.
# Local paths are active; server paths are kept below as comments.
# source() executes the script directly, like the old scripts.
library(hdf5r)       # The only new package: read the h5ad file.
library(ggplot2)
library(reshape2)

rm(list = ls())

############################ Paths and parameters ############################
# Local
input_file <- '/Users/erganshang/Desktop/pert_map/Xin/data/raw/subset_adata_SCARF_10celltypes_10perturbations.h5ad'
result_path <- '/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/'
plot_path <- '/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/figures/'

# Server: to run 01 there, comment the three local lines and uncomment these.
# input_file <- '/raid6/home/eshang/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/raw/subset_adata_SCARF_10celltypes_10perturbations.h5ad'
# result_path <- '/raid6/home/eshang/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/'
# plot_path <- '/raid6/home/eshang/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/figures/'
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

cell_type <- '151 TH Prkcd Grin2c Glut'
ctrl_label <- 'Non_target'
pert <- c('A', 'B', 'C', 'D', 'E', 'F')
n_samples <- 200
ratio <- 2
n_modules <- 30
module_sz <- 50
n_features <- n_modules * module_sz
sample_seed <- 4
signal_seeds <- setNames(1:6, pert)
data_file <- paste0(result_path, 'Xin_151_200_400.Rdata')

############################ Small data-reading helpers #######################
# h5ad categorical annotations are stored as categories + zero-based codes.
read_annotation <- function(f, group, key) {
  obj <- f[[paste0(group, '/', key)]]
  if (inherits(obj, 'H5Group')) {
    codes <- obj[['codes']][]
    categories <- obj[['categories']][]
    out <- rep(NA_character_, length(codes))
    ok <- codes >= 0
    out[ok] <- as.character(categories[codes[ok] + 1L])
    out
  } else obj[]
}

# Sum all eligible ctrl cells in blocks, without loading the full h5ad matrix.
sum_ctrl_genes <- function(f, eligible, n_genes, block_size = 1024L) {
  ptr <- as.numeric(f[['X/indptr']][])
  totals <- numeric(n_genes)
  for (start in seq.int(1L, length(eligible), by = block_size)) {
    end <- min(start + block_size - 1L, length(eligible))
    if (!any(eligible[start:end])) next
    first <- ptr[start] + 1
    last <- ptr[end + 1L]
    if (first > last) next
    keep <- rep(eligible[start:end], times = diff(ptr[start:(end + 1L)]))
    values <- f[['X/data']][first:last][keep]
    columns <- f[['X/indices']][first:last][keep] + 1L
    stopifnot(all(is.finite(values)), all(values >= 0), all(values == round(values)))
    if (length(values) > 0) {
      sums <- rowsum(values, columns, reorder = FALSE)
      idx <- as.integer(rownames(sums))
      totals[idx] <- totals[idx] + sums[, 1L]
    }
  }
  totals
}

# Read the sampled cells and selected genes from CSR-encoded raw counts.
read_counts <- function(f, rows, columns, n_genes) {
  ptr <- as.numeric(f[['X/indptr']][])
  lookup <- integer(n_genes)
  lookup[columns] <- seq_along(columns)
  out <- matrix(0, length(rows), length(columns))
  for (i in seq_along(rows)) {
    first <- ptr[rows[i]] + 1
    last <- ptr[rows[i] + 1L]
    if (first > last) next
    values <- f[['X/data']][first:last]
    positions <- lookup[f[['X/indices']][first:last] + 1L]
    keep <- positions > 0L
    if (any(keep)) {
      sums <- rowsum(values[keep], positions[keep], reorder = FALSE)
      out[i, as.integer(rownames(sums))] <- sums[, 1L]
    }
  }
  stopifnot(all(is.finite(out)), all(out >= 0), all(out == round(out)))
  out
}

# Same dilution helper as the original folder, now with baseline multiplier = 1.
.expand_module_vector <- function(x, n_modules, module_size, n_features, baseline, seed = NULL) {
  if (!is.null(seed)) {
    return(withr::with_seed(seed, .expand_module_vector(
      x, n_modules, module_size, n_features, baseline, seed = NULL)))
  }
  if (length(x) == 1L) {
    rep(x, n_features)
  } else if (length(x) == n_modules) {
    out <- numeric(n_features)
    for (i in seq_len(n_modules)) {
      vals <- rep(x[i], module_size)
      n_replace <- ceiling(0.8 * module_size)
      idx_replace <- sample(seq_len(module_size), n_replace)
      vals[idx_replace] <- baseline
      start <- (i - 1) * module_size + 1
      end <- i * module_size
      out[start:end] <- vals
    }
    out
  } else if (length(x) == n_features) {
    x
  } else stop('Length must be 1, n_modules, or n_features.')
}

############################ Read the selected ctrl pool ######################
stopifnot(file.exists(input_file))
f <- H5File$new(input_file, mode = 'r')
tryCatch({
  stopifnot(identical(h5attr(f[['X']], 'encoding-type'), 'csr_matrix'))
  shape <- as.integer(h5attr(f[['X']], 'shape'))
  cell_ids <- as.character(read_annotation(f, 'obs', h5attr(f[['obs']], '_index')))
  gene_names <- as.character(read_annotation(f, 'var', h5attr(f[['var']], '_index')))
  gene_target <- read_annotation(f, 'obs', 'gene_target')
  predicted_group <- read_annotation(f, 'obs', 'predicted_group')
  batch <- read_annotation(f, 'obs', 'batch')
  stopifnot(length(cell_ids) == shape[1], length(gene_names) == shape[2],
    !anyNA(cell_ids), !anyDuplicated(cell_ids),
    !anyNA(gene_names), !anyDuplicated(gene_names))

  eligible <- !is.na(gene_target) & gene_target == ctrl_label &
    !is.na(predicted_group) & predicted_group == cell_type
  ctrl_idx <- which(eligible)
  stopifnot(length(ctrl_idx) >= n_samples * (length(pert) + ratio))
  cat('Eligible 151 ctrl cells:', length(ctrl_idx), '\n')

  ######################## Select genes BEFORE sampling cells #################
  total_gene_expression <- sum_ctrl_genes(f, eligible, shape[2])
  # Descending expression; ties follow the original gene-column order.
  gene_idx <- order(-total_gene_expression, seq_along(total_gene_expression),
    method = 'radix')[seq_len(n_features)]
  stopifnot(all(total_gene_expression[gene_idx] > 0))

  ######################## Sample and split 1600 unique cells #################
  # Fix both candidate ordering and RNG algorithm for repeatable sampling.
  ctrl_idx <- ctrl_idx[order(cell_ids[ctrl_idx], method = 'radix')]
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  set.seed(sample_seed)
  selected_idx <- sample(ctrl_idx, n_samples * (length(pert) + ratio), replace = FALSE)
  group_labels <- c(rep(pert, each = n_samples), rep('NT', n_samples * ratio))
  original_counts <- read_counts(f, selected_idx, gene_idx, shape[2])
  dimnames(original_counts) <- list(cell_ids[selected_idx], gene_names[gene_idx])
}, finally = f$close_all())

############################ Define the ORIGINAL effect multipliers ###########
# These are multipliers, not absolute expression means.
fold_change_mods_list <- list(
  A = c(1.4, 1.5, 1.45, 1.4, 1.3, rep(1, n_modules - 5)),
  B = c(0.95, 0.5, 0.45, 0.55, 0.95, rep(1, n_modules - 5)),
  C = c(rep(1, 5), 1.4, 1.5, 1.45, 1.4, 1.3, rep(1, 20)),
  D = c(0.5, 1.2, 1.5, 0.6, 1.45, rep(1, n_modules - 5)),
  E = c(1, 1.45, 1.4, 1, 1, 1, 0.45, 0.5, rep(1, n_modules - 8)),
  F = c(rep(1, 11), 1.4, 1.5, 1.45, 1.4, 1.3, rep(1, 14))
)

# Within each module, set 40/50 gene multipliers back to 1.
# A-F use the same seeds 1-6 as the original dilution helper.
fold_change_list <- list()
for (group in pert) {
  fold_change_list[[group]] <- .expand_module_vector(
    fold_change_mods_list[[group]], n_modules, module_sz, n_features,
    baseline = 1, seed = signal_seeds[[group]])
  names(fold_change_list[[group]]) <- colnames(original_counts)
}

############################ Multiply the raw counts, then ROUND ###############
modified_counts <- original_counts
for (group in pert) {
  rows <- which(group_labels == group)
  multiplier <- fold_change_list[[group]]
  changed <- which(multiplier != 1)
  modified_counts[rows, changed] <- round(
    sweep(original_counts[rows, changed, drop = FALSE], 2, multiplier[changed], '*'))
}
stopifnot(identical(original_counts[group_labels == 'NT', , drop = FALSE],
                    modified_counts[group_labels == 'NT', , drop = FALSE]))

############################ Save ONE dataset for scripts 02-05 ###############
grouping_vector <- rep(seq_len(n_modules), each = module_sz)
names(grouping_vector) <- colnames(original_counts)
cell_info <- data.frame(cell_id = cell_ids[selected_idx], source_row = selected_idx,
  group = group_labels, cell_type = predicted_group[selected_idx], batch = batch[selected_idx])
gene_info <- data.frame(gene = gene_names[gene_idx], source_column = gene_idx,
  expression_rank = seq_len(n_features), ctrl_total_expression = total_gene_expression[gene_idx],
  module = grouping_vector, position_in_module = rep(seq_len(module_sz), n_modules),
  row.names = NULL)
effect_design <- do.call(rbind, lapply(pert, function(group) data.frame(
  group = group, gene = gene_info$gene, module = gene_info$module,
  fold_change = unname(fold_change_list[[group]]),
  is_signal = unname(fold_change_list[[group]]) != 1)))

tr_sample <- lapply(pert, function(group) modified_counts[group_labels == group, , drop = FALSE])
names(tr_sample) <- pert
sim <- list(
  control = modified_counts[group_labels == 'NT', , drop = FALSE],
  treatment = tr_sample,
  original_counts = original_counts,
  cell_info = cell_info, gene_info = gene_info,
  fold_change = fold_change_list, module_fold_change = fold_change_mods_list,
  control_pool_cell_ids = cell_ids[ctrl_idx],
  source_md5 = unname(tools::md5sum(input_file)),
  seeds = list(cell_sample = sample_seed, signal_genes = signal_seeds)
)

# A rerun must reconstruct exactly the same data, cells and signal genes.
if (file.exists(data_file)) {
  previous <- new.env()
  load(data_file, envir = previous)
  if (!identical(sim, previous$sim)) stop('Recomputed data differ from the saved dataset. Check input, parameters and R version.')
}
save(sim, file = data_file, compress = FALSE)
write.csv(cell_info, paste0(result_path, 'selected_cells.csv'), row.names = FALSE)
write.csv(gene_info, paste0(result_path, 'selected_genes_and_modules.csv'), row.names = FALSE)
write.csv(effect_design, paste0(result_path, 'effect_design_all_genes.csv'), row.names = FALSE)
write.csv(effect_design[effect_design$is_signal, ], paste0(result_path, 'modified_gene_list.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(result_path, 'sessionInfo_01.txt'))

############################ Original Poisson GLM effect calculation ##########
tr_all <- do.call(rbind, sim$treatment)
Y <- rbind(tr_all, sim$control)
Y_df <- as.data.frame(Y)
Y_df$pert <- factor(c(rep(pert, each = n_samples), rep('NT', n_samples * ratio)))
Y_df$pert <- relevel(Y_df$pert, ref = 'NT')

mat_res <- matrix(0, length(pert), ncol(Y), dimnames = list(pert, colnames(Y)))
for (i in seq_len(ncol(Y))) {
  model <- glm(Y_df[, i] ~ pert, family = poisson, data = Y_df)
  for (group in pert) {
    mat_res[group, i] <- exp(coef(model)[paste0('pert', group)] + coef(model)['(Intercept)']) /
      exp(coef(model)['(Intercept)'])
  }
}
mat_res <- log2(mat_res)

program_mat <- matrix(0, length(pert), n_modules,
  dimnames = list(pert, as.character(seq_len(n_modules))))
for (i in seq_len(nrow(program_mat))) {
  for (j in seq_len(n_modules)) {
    start <- (j - 1) * module_sz + 1
    end <- j * module_sz
    program_mat[i, j] <- mean(mat_res[i, start:end])
  }
}
save(mat_res, program_mat, file = paste0(result_path, 'poisson_effects.Rdata'), compress = FALSE)
write.csv(mat_res, paste0(result_path, 'gene_log2fc.csv'))
write.csv(program_mat, paste0(result_path, 'module_mean_log2fc.csv'))

mat_long <- melt(mat_res)
p_after <- ggplot(mat_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  theme_minimal() + labs(x = 'Gene', y = 'Perturbation', fill = 'log2FC') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
ggsave(paste0(plot_path, 'logFoldChange_Xin_gene_pert_poisson_200_400.pdf'), p_after, width = 12, height = 4)

mat_long <- melt(program_mat)
p <- ggplot(mat_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = '#013e75', mid = '#FAFAFA', high = '#f5b70a', midpoint = 0) +
  theme_minimal() + labs(x = 'Program', y = 'Perturbation', fill = 'log2FC') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank(), legend.position = 'top') +
  guides(fill = guide_colorbar(
    direction = 'horizontal',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = 20,
    barheight = 0.5
  ))
ggsave(paste0(plot_path, 'logFoldChange_Xin_program_pert_poisson_200_400.pdf'), p, width = 9, height = 4)
cat('01 complete: saved the fixed dataset and two effect heatmaps.\n')
