rm(list = ls())

data_path <- file.path('..', 'data', 'raw')
result_path <- file.path('..', 'data', 'final', '11_pairwise_convergence')
plot_path <- file.path('..', 'figures', '13_splatter_pilot')
data_file <- file.path(data_path, '10_splatter_pilot_60_120.Rdata')
result_one_file <- file.path(result_path, '11_A_B_one.Rdata')
result_two_file <- file.path(result_path, '11_A_B_two.Rdata')
module_log2fc_file <- file.path(data_path, '12_module_mean_log2fc.csv')
pair_summary_file <- file.path(result_path, '13_convergence_pair_summary.csv')
combined_files <- c(
  union = file.path(plot_path, '13_splatter_pilot_combined_union.pdf'),
  intersect = file.path(plot_path, '13_splatter_pilot_combined_intersect.pdf')
)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
load(data_file)
current_input_md5 <- unname(tools::md5sum(data_file))

############################ Recreate the heatmap from script 12 data #########
module_mat <- as.matrix(read.csv(module_log2fc_file,
  row.names = 1, check.names = FALSE))
module_data <- as.data.frame(as.table(module_mat))
names(module_data) <- c('Perturbation', 'Module', 'log2FC')
p_heatmap <- ggplot(module_data,
  aes(x = Module, y = Perturbation, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = '#013e75', mid = '#FAFAFA', high = '#f5b70a',
    midpoint = 0) +
  labs(x = 'Gene module', y = 'Perturbation', fill = 'log2FC') +
  theme_minimal() + theme(panel.grid = element_blank())

############################ Collect the two directions #######################
grouping_vector <- setNames(sim$gene_info$module, sim$gene_info$gene)
active_gene_num <- 5
p_thres <- 0.05

load(result_one_file)
stopifnot(identical(input_md5, current_input_md5))
result1 <- result
load(result_two_file)
stopifnot(identical(input_md5, current_input_md5))
result2 <- result

features1 <- collect_active_features(result1,
  group = grouping_vector, group_threshold = active_gene_num)
features2 <- collect_active_features(result2,
  group = grouping_vector, group_threshold = active_gene_num)
groups1 <- features1$active_groups
groups2 <- features2$active_groups

df <- data.frame(pert_1 = 'A', pert_2 = 'B',
  p_forward = result1$p_value, p_reverse = result2$p_value,
  p_value = max(result1$p_value, result2$p_value),
  TestStat = min(abs(result1$test_statistic), abs(result2$test_statistic)),
  module_union = length(union(groups1, groups2)),
  module_intersect = length(intersect(groups1, groups2)),
  active_modules_forward = paste(sort(as.integer(groups1)), collapse = ';'),
  active_modules_reverse = paste(sort(as.integer(groups2)), collapse = ';'))
df$p_bonferroni <- pmin(1, df$p_value)  # One pair in this pilot.
df$Convergent <- ifelse(df$p_value <= p_thres, 'Adj.p<0.05', 'p>0.05')
write.csv(df, pair_summary_file, row.names = FALSE)

############################ Combined heatmap and convergence plots ###########
for (count_method in c('union', 'intersect')) {
  df_upper <- df
  df_lower <- df
  df_lower$pert_1 <- 'B'
  df_lower$pert_2 <- 'A'
  df_lower$module_num <- df[[paste0('module_', count_method)]]
  df_upper$Convergent <- factor(df_upper$Convergent,
    levels = c('Adj.p<0.05', 'p>0.05'))

  p_convergence <- ggplot() +
    geom_point(data = df_upper,
      aes(x = pert_1, y = pert_2, color = Convergent), size = 8) +
    scale_color_manual(name = 'Convergence', drop = FALSE,
      values = c('Adj.p<0.05' = '#f5b70a', 'p>0.05' = '#EEEEEE')) +
    geom_point(data = df_lower,
      aes(x = pert_1, y = pert_2, size = module_num),
      alpha = 0.7, shape = 16, color = '#663B8C') +
    scale_size_continuous(name = 'Module Number', range = c(0, 10)) +
    scale_x_discrete(limits = c('A', 'B')) +
    scale_y_discrete(limits = c('A', 'B')) +
    labs(title = 'Perturbation Pair Convergence',
      x = 'Perturbation 1', y = 'Perturbation 2') +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

  grDevices::pdf(combined_files[[count_method]], width = 10.5, height = 4)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2,
    widths = grid::unit(c(5.5, 5), 'null'))))
  print(p_heatmap, vp = grid::viewport(layout.pos.row = 1,
    layout.pos.col = 1), newpage = FALSE)
  print(p_convergence, vp = grid::viewport(layout.pos.row = 1,
    layout.pos.col = 2), newpage = FALSE)
  grDevices::dev.off()
}
