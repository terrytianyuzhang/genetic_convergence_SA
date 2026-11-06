# 04: Same pooled standardization and Chen test as the original script.
if (Sys.info()[['sysname']] != 'Linux') stop('Run Chen tests on the Linux server, not on the local computer.')
library(highmean)

rm(list = ls())

data_path <- '/raid6/home/eshang/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/'
result_path <- paste0(data_path, 'Chen/')
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
data_file <- paste0(data_path, 'Xin_151_200_400.Rdata')
load(data_file)
input_md5 <- unname(tools::md5sum(data_file))

ctrl <- sim$control
pert <- c('A', 'B', 'C', 'D', 'E', 'F')
n_modules <- 30
module_sz <- 50
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
p_val_mat <- matrix(0, nrow = length(pert), ncol = n_modules,
  dimnames = list(pert, as.character(seq_len(n_modules))))

for (i in seq_along(pert)) {
  pert_evaluate <- pert[i]
  tr <- sim$treatment[[pert_evaluate]]
  print(paste0('For ', pert_evaluate))
  for (j in seq_len(n_modules)) {
    start <- (j - 1) * module_sz + 1
    end <- j * module_sz
    ctrl_use <- ctrl[, start:end, drop = FALSE]
    tr_use <- tr[, start:end, drop = FALSE]

    ######################## Original pooled standardization ##################
    combined <- rbind(ctrl_use, tr_use)
    pooled_mean <- colMeans(combined)
    pooled_sd <- apply(combined, 2, sd)
    stopifnot(all(is.finite(pooled_sd)), all(pooled_sd > 0))
    normalized <- scale(combined, center = pooled_mean, scale = pooled_sd)
    ctrl_evaluate <- normalized[seq_len(nrow(ctrl_use)), , drop = FALSE]
    tr_evaluate <- normalized[nrow(ctrl_use) + seq_len(nrow(tr_use)), , drop = FALSE]

    run_seed <- 50000 + (i - 1) * n_modules + j
    set.seed(run_seed)
    result <- apval_Chen2010(ctrl_evaluate, tr_evaluate, eq.cov = FALSE)
    result_file <- paste0(result_path, pert_evaluate, '_ctrl_module_', j, '.Rdata')
    if (file.exists(result_file)) {
      previous <- new.env()
      load(result_file, envir = previous)
      if (!identical(input_md5, previous$input_md5) || !identical(result, previous$result))
        stop('Chen result/input differ from the previous run: ', pert_evaluate, ' module ', j)
    }
    save(result, run_seed, input_md5, file = result_file, compress = FALSE)
    p_val_mat[pert_evaluate, j] <- result$pval
  }
}
save(p_val_mat, input_md5, file = paste0(result_path, 'p_val_mat.Rdata'), compress = FALSE)
write.csv(p_val_mat, paste0(result_path, 'p_val_mat.csv'))
writeLines(capture.output(sessionInfo()), paste0(data_path, 'sessionInfo_04.txt'))
cat('04 complete: 6 perturbations x 30 modules = 180 Chen tests.\n')
