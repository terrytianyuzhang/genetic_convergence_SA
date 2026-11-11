rm(list = ls())

# INPUT
data_path <- file.path('..', 'data', 'raw')
data_file <- file.path(data_path, '10_splatter_data.Rdata')

convergence_path <- file.path('..', 'data', 'final', '11_pairwise_convergence')
pert <- c('A', 'B', 'C', 'D')
comb_2 <- combn(pert, 2)
pair_names <- apply(comb_2, 2, paste, collapse = '_')
result_files <- setNames(lapply(pair_names, function(pair) c(
  one = file.path(convergence_path, paste0('11_', pair, '_one.Rdata')),
  two = file.path(convergence_path, paste0('11_', pair, '_two.Rdata'))
)), pair_names)

module_log2fc_file <- file.path('..', 'data', 'final', '12_heatmap', '12_module_mean_log2fc.csv')

# OUTPUT
pair_summary_file <- file.path('..', 'data', 'final', '13_convergence_pair_summary.csv')

plot_path <- file.path('..', 'figures', '13_splatter_pilot')
combined_files <- c(
  union = file.path(plot_path, '13_splatter_pilot_combined_union.pdf'),
  intersect = file.path(plot_path, '13_splatter_pilot_combined_intersect.pdf')
)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
load(data_file)
stopifnot(identical(names(sim$treatment), pert))
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
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank(), legend.position = 'top') +
  guides(fill = guide_colorbar(direction = 'horizontal',
    title.position = 'top', title.hjust = 0.5,
    barwidth = 20, barheight = 0.5))

############################ Collect the two directions #######################
grouping_vector <- setNames(sim$gene_info$module, sim$gene_info$gene)
active_gene_num <- 5
p_thres <- 0.05

pair_num <- ncol(comb_2)
df <- data.frame(pert_1 = comb_2[1, ], pert_2 = comb_2[2, ],
  p_forward = NA_real_, p_reverse = NA_real_, p_value = NA_real_,
  TestStat = NA_real_, module_union = 0L, module_intersect = 0L,
  active_modules_forward = '', active_modules_reverse = '')

for (i in seq_len(pair_num)) {
  pair <- pair_names[i]
  load(result_files[[pair]][['one']])
  stopifnot(identical(input_md5, current_input_md5))
  result1 <- result
  load(result_files[[pair]][['two']])
  stopifnot(identical(input_md5, current_input_md5))
  result2 <- result

  features1 <- collect_active_features(result1,
    group = grouping_vector, group_threshold = active_gene_num)
  features2 <- collect_active_features(result2,
    group = grouping_vector, group_threshold = active_gene_num)
  groups1 <- features1$active_groups
  groups2 <- features2$active_groups

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
write.csv(df, pair_summary_file, row.names = FALSE)

############################ Combined heatmap and convergence plots ###########
for (count_method in c('union', 'intersect')) {
  df_upper <- df
  df_lower <- df
  df_lower$pert_1 <- df$pert_2
  df_lower$pert_2 <- df$pert_1
  df_lower$module_num <- df[[paste0('module_', count_method)]]
  df_upper$Convergent <- factor(df_upper$Convergent,
    levels = c('Adj.p<0.05', 'p>0.05', 'p<0.05'))

  p_convergence <- ggplot() +
    geom_point(data = df_upper,
      aes(x = pert_1, y = pert_2, color = Convergent), size = 8) +
    scale_color_manual(name = 'Convergence', drop = FALSE,
      values = c('Adj.p<0.05' = '#f5b70a', 'p>0.05' = '#EEEEEE',
        'p<0.05' = '#013e75')) +
    geom_point(data = df_lower,
      aes(x = pert_1, y = pert_2, size = module_num),
      alpha = 0.7, shape = 16, color = '#663B8C') +
    scale_size_continuous(name = 'Module Number', range = c(0, 10)) +
    scale_x_discrete(limits = pert) +
    scale_y_discrete(limits = pert) +
    coord_fixed(ratio = 1) +
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
