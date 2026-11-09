rm(list = ls())

# INPUT
data_path <- file.path('..', 'data', 'raw')
data_file <- file.path(data_path, '10_splatter_pilot_60_120.Rdata')


# OUTPUT
result_path <- file.path('..', 'data', 'final', '11_pairwise_convergence')
result_files <- c(
  one = file.path(result_path, '11_A_B_one.Rdata'),
  two = file.path(result_path, '11_A_B_two.Rdata')
)
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

load(data_file)
input_md5 <- unname(tools::md5sum(data_file))
ctrl <- sim$control
grouping_vector <- setNames(sim$gene_info$module, sim$gene_info$gene)
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')

for (direction in c('one', 'two')) {
  tr_1 <- if (direction == 'one') sim$treatment$A else sim$treatment$B
  tr_2 <- if (direction == 'one') sim$treatment$B else sim$treatment$A
  run_seed <- if (direction == 'one') 11001 else 11002
  set.seed(run_seed)
  result <- convergence_testing(
    control = ctrl, treatment1 = tr_1, treatment2 = tr_2,
    pca_method = 'dense_pca', classifier_method = 'group_lasso',
    lambda_type = 'lambda.1se', n_folds = 10,
    group = grouping_vector, verbose = FALSE, standardize_feature = TRUE
  )
  save(result, run_seed, input_md5,
    file = result_files[[direction]], compress = FALSE)
}
