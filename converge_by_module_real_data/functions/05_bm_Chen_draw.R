# 05: Same shared-significant-module criterion as the original Chen plot.
library(ggplot2)

rm(list = ls())

data_path <- '/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/'
result_path <- paste0(data_path, 'Chen/')
plot_path <- '/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/figures/'
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
current_input_md5 <- unname(tools::md5sum(paste0(data_path, 'Xin_151_200_400.Rdata')))
load(paste0(result_path, 'p_val_mat.Rdata'))
stopifnot(identical(input_md5, current_input_md5))

pert <- c('A', 'B', 'C', 'D', 'E', 'F')
pair_num <- 15
n_modules <- ncol(p_val_mat)
threshold <- 0.05 / (pair_num * n_modules)
df <- data.frame(pert_1 = rep(pert, each = length(pert)),
  pert_2 = rep(pert, length(pert)), convergence = NA_character_, module_num = NA_integer_)

for (i in seq_len(nrow(df))) {
  pert_name_1 <- df$pert_1[i]
  pert_name_2 <- df$pert_2[i]
  if (pert_name_1 == pert_name_2) next
  p_mat_sub <- p_val_mat[c(pert_name_1, pert_name_2), , drop = FALSE]
  count <- 0
  for (j in seq_len(n_modules)) {
    if (p_mat_sub[pert_name_1, j] < threshold & p_mat_sub[pert_name_2, j] < threshold) {
      count <- count + 1
    }
  }
  df$module_num[i] <- count
  df$convergence[i] <- if (count > 0) 'Yes' else 'No'
}
df$convergence <- factor(df$convergence, levels = c('Yes', 'No'))
df$pert_1 <- factor(df$pert_1, levels = pert)
df$pert_2 <- factor(df$pert_2, levels = pert)
df_no_diag <- df[df$pert_1 != df$pert_2, ]
df_upper <- df_no_diag[match(df_no_diag$pert_1, pert) < match(df_no_diag$pert_2, pert), ]
df_lower <- df_no_diag[match(df_no_diag$pert_1, pert) > match(df_no_diag$pert_2, pert), ]
write.csv(df_upper, paste0(result_path, 'Chen_pair_summary.csv'), row.names = FALSE)

p <- ggplot() +
  geom_point(data = df_upper, aes(x = pert_1, y = pert_2, color = convergence), size = 6) +
  scale_color_manual(name = 'Convergence', drop = FALSE,
    values = c('Yes' = '#f5b70a', 'No' = '#EEEEEE')) +
  geom_point(data = df_lower, aes(x = pert_1, y = pert_2, size = module_num),
    alpha = 0.7, shape = 16, color = '#6a0dad') +
  scale_size_continuous(name = 'Module Number', range = c(1, 10)) +
  scale_x_discrete(limits = pert) + scale_y_discrete(limits = pert) +
  labs(title = 'Perturbation Pair Convergence', x = 'Perturbation 1', y = 'Perturbation 2') +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(plot_path, 'Chen_revise.pdf'), p, width = 5.5, height = 4)
writeLines(capture.output(sessionInfo()), paste0(data_path, 'sessionInfo_05.txt'))
cat('05 complete: saved Chen_revise.pdf and the 15-pair summary.\n')
