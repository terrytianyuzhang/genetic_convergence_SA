rm(list = ls())

# INPUT
data_path <- file.path('..', 'data', 'raw')
data_file <- file.path(data_path, '10_splatter_data.Rdata')

# OUTPUT
result_path <- file.path('..', 'data', 'final', '11_pairwise_convergence')
pert <- c('A', 'B', 'C', 'D')
comb_2 <- combn(pert, 2)
pair_names <- apply(comb_2, 2, paste, collapse = '_')
result_files <- setNames(lapply(pair_names, function(pair) c(
  one = file.path(result_path, paste0('11_', pair, '_one.Rdata')),
  two = file.path(result_path, paste0('11_', pair, '_two.Rdata'))
)), pair_names)
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

load(data_file)
stopifnot(identical(names(sim$treatment), pert))
input_md5 <- unname(tools::md5sum(data_file))
ctrl <- sim$control
grouping_vector <- setNames(sim$gene_info$module, sim$gene_info$gene)
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')

for (i in seq_len(ncol(comb_2))) {
  pair <- pair_names[i]
  for (direction in c('one', 'two')) {
    tr_1 <- if (direction == 'one') sim$treatment[[comb_2[1, i]]] else
      sim$treatment[[comb_2[2, i]]]
    tr_2 <- if (direction == 'one') sim$treatment[[comb_2[2, i]]] else
      sim$treatment[[comb_2[1, i]]]
    run_seed <- 11000 + 2 * i - (if (direction == 'one') 1 else 0)
    set.seed(run_seed)
    result <- convergence_testing(
      control = ctrl, treatment1 = tr_1, treatment2 = tr_2,
      pca_method = 'dense_pca', classifier_method = 'group_lasso',
      lambda_type = 'lambda.1se', n_folds = 10,
      group = grouping_vector, verbose = FALSE, standardize_feature = TRUE
    )
    save(result, run_seed, input_md5,
      file = result_files[[pair]][[direction]], compress = FALSE)
  }
}
