# 03: BOTH union and intersect plots. Paper preference: INTERSECT.
# The pair p-value is always MAX of the two directional p-values.
library(HMC)
library(ggplot2)

rm(list = ls())

data_path <- '/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/data/'
result_path <- paste0(data_path, 'convergence/')
plot_path <- '/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/functions/NT_effect/normal_artificial/redo_by_Xin_data/01/figures/'
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
data_file <- paste0(data_path, 'Xin_151_200_400.Rdata')
load(data_file)
current_input_md5 <- unname(tools::md5sum(data_file))

all_names <- c('A', 'B', 'C', 'D', 'E', 'F')
grouping_vector <- setNames(sim$gene_info$module, sim$gene_info$gene)
active_gene_num <- 5
p_thres <- 0.05
comb_2 <- combn(all_names, 2)
pair_num <- ncol(comb_2)
df <- data.frame(pert_1 = comb_2[1, ], pert_2 = comb_2[2, ],
  p_forward = NA_real_, p_reverse = NA_real_, p_value = NA_real_, TestStat = NA_real_,
  module_union = 0L, module_intersect = 0L,
  active_modules_forward = '', active_modules_reverse = '')

############################ Collect the two directions #######################
for (i in seq_len(pair_num)) {
  a <- comb_2[1, i]
  b <- comb_2[2, i]
  load(paste0(result_path, a, '_', b, '_one.Rdata'))
  stopifnot(identical(input_md5, current_input_md5))
  result1 <- result
  load(paste0(result_path, a, '_', b, '_two.Rdata'))
  stopifnot(identical(input_md5, current_input_md5))
  result2 <- result

  feature_collected_1 <- collect_active_features(result1,
    group = grouping_vector, group_threshold = active_gene_num)
  feature_collected_2 <- collect_active_features(result2,
    group = grouping_vector, group_threshold = active_gene_num)
  groups1 <- feature_collected_1$active_groups
  groups2 <- feature_collected_2$active_groups

  df$p_forward[i] <- result1$p_value
  df$p_reverse[i] <- result2$p_value
  df$p_value[i] <- max(result1$p_value, result2$p_value)
  df$TestStat[i] <- min(abs(result1$test_statistic), abs(result2$test_statistic))
  df$module_union[i] <- length(union(groups1, groups2))
  df$module_intersect[i] <- length(intersect(groups1, groups2))
  df$active_modules_forward[i] <- paste(sort(as.integer(groups1)), collapse = ';')
  df$active_modules_reverse[i] <- paste(sort(as.integer(groups2)), collapse = ';')
}
df$p_bonferroni <- pmin(1, df$p_value * pair_num)
df$Convergent <- ifelse(df$p_value <= p_thres / pair_num, 'Adj.p<0.05',
  ifelse(df$p_value >= p_thres, 'p>0.05', 'p<0.05'))
write.csv(df, paste0(data_path, 'convergence_pair_summary.csv'), row.names = FALSE)

############################ Plot both count methods ##########################
for (count_method in c('union', 'intersect')) {
  df_upper <- df
  df_lower <- df
  df_lower$pert_1 <- df$pert_2
  df_lower$pert_2 <- df$pert_1
  df_lower$module_num <- df[[paste0('module_', count_method)]]
  df_upper$Convergent <- factor(df_upper$Convergent,
    levels = c('Adj.p<0.05', 'p>0.05', 'p<0.05'))

  p <- ggplot() +
    geom_point(data = df_upper, aes(x = pert_1, y = pert_2, color = Convergent), size = 8) +
    scale_color_manual(name = 'Convergence', drop = FALSE,
      breaks = c('Adj.p<0.05', 'p>0.05'),
      values = c('Adj.p<0.05' = '#f5b70a', 'p>0.05' = '#EEEEEE', 'p<0.05' = '#013e75')) +
    geom_point(data = df_lower, aes(x = pert_1, y = pert_2, size = module_num),
      alpha = 0.7, shape = 16, color = '#663B8C') +
    scale_size_continuous(name = 'Module Number', range = c(0, 10)) +
    scale_x_discrete(limits = all_names) + scale_y_discrete(limits = all_names) +
    labs(title = 'Perturbation Pair Convergence',
      x = 'Perturbation 1', y = 'Perturbation 2') +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(plot_path, 'Xin_NT_3_color_raw_grlasso_1se_200_400_', count_method, '.pdf'),
    p, width = if (count_method == 'union') 10 else 5.5,
    height = if (count_method == 'union') 8 else 4)
}
writeLines(capture.output(sessionInfo()), paste0(data_path, 'sessionInfo_03.txt'))
cat('03 complete: max p-value, union and intersect figures. Paper preference: intersect.\n')
