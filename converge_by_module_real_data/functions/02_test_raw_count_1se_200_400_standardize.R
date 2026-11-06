# 02: Same HMC convergence_testing() calls and parameters as the original script.
# Run on the SERVER only. This script reads the dataset saved by 01.
if (Sys.info()[['sysname']] != 'Linux') stop('Run HMC on the Linux server, not on the local computer.')
library(HMC)

rm(list = ls())

data_path <- '/raid6/home/eshang/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/'
result_path <- paste0(data_path, 'convergence/')
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
data_file <- paste0(data_path, 'Xin_151_200_400.Rdata')
load(data_file)
input_md5 <- unname(tools::md5sum(data_file))

ctrl <- sim$control
pert <- c('A', 'B', 'C', 'D', 'E', 'F')
grouping_vector <- setNames(sim$gene_info$module, sim$gene_info$gene)
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')

comb_2 <- combn(pert, 2)
seed_table <- data.frame(pair = rep(apply(comb_2, 2, paste, collapse = '_'), each = 2),
  direction = rep(c('one', 'two'), ncol(comb_2)), seed = 10000 + seq_len(2 * ncol(comb_2)))
write.csv(seed_table, paste0(data_path, 'convergence_seed_table.csv'), row.names = FALSE)

for (i in seq_len(ncol(comb_2))) {
  comb <- comb_2[, i]
  print(paste0(comb[1], ' and ', comb[2]))

  ############################ Forward direction ##############################
  tr_1 <- sim$treatment[[comb[1]]]
  tr_2 <- sim$treatment[[comb[2]]]
  run_seed <- 10000 + 2 * i - 1
  set.seed(run_seed)
  result <- convergence_testing(
    control = ctrl, treatment1 = tr_1, treatment2 = tr_2,
    pca_method = 'dense_pca', classifier_method = 'group_lasso',
    lambda_type = 'lambda.1se', n_folds = 10,
    group = grouping_vector, verbose = FALSE, standardize_feature = TRUE
  )
  result_file <- paste0(result_path, comb[1], '_', comb[2], '_one.Rdata')
  if (file.exists(result_file)) {
    previous <- new.env()
    load(result_file, envir = previous)
    if (!identical(input_md5, previous$input_md5) || !identical(result, previous$result))
      stop('Forward result/input differ from the previous run: ', comb[1], '_', comb[2])
  }
  save(result, run_seed, input_md5, file = result_file, compress = FALSE)

  ############################ Reverse direction ##############################
  tr_1 <- sim$treatment[[comb[2]]]
  tr_2 <- sim$treatment[[comb[1]]]
  run_seed <- 10000 + 2 * i
  set.seed(run_seed)
  result <- convergence_testing(
    control = ctrl, treatment1 = tr_1, treatment2 = tr_2,
    pca_method = 'dense_pca', classifier_method = 'group_lasso',
    lambda_type = 'lambda.1se', n_folds = 10,
    group = grouping_vector, verbose = FALSE, standardize_feature = TRUE
  )
  result_file <- paste0(result_path, comb[1], '_', comb[2], '_two.Rdata')
  if (file.exists(result_file)) {
    previous <- new.env()
    load(result_file, envir = previous)
    if (!identical(input_md5, previous$input_md5) || !identical(result, previous$result))
      stop('Reverse result/input differ from the previous run: ', comb[1], '_', comb[2])
  }
  save(result, run_seed, input_md5, file = result_file, compress = FALSE)
}
writeLines(capture.output(sessionInfo()), paste0(data_path, 'sessionInfo_02.txt'))
cat('02 complete: 15 perturbation pairs x 2 directions.\n')
