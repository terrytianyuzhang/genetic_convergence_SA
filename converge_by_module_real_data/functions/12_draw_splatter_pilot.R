rm(list = ls())

data_path <- file.path('..', 'data', 'raw')
plot_path <- file.path('..', 'figures', 'splatter_pilot')
data_file <- file.path(data_path, '10_splatter_pilot_60_120.Rdata')
module_log2fc_file <- file.path(data_path, '12_module_mean_log2fc.csv')
heatmap_file <- file.path(plot_path, '12_splatter_pilot_module_heatmap.pdf')
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
load(data_file)

############################ Plot observed module log2 fold changes ###########
# A 0.5 pseudocount keeps the display finite when a gene is all zero in a group.
pert <- c('A', 'B')
n_modules <- length(unique(sim$gene_info$module))
control_mean <- colMeans(sim$control)
module_mat <- matrix(0, length(pert), n_modules,
  dimnames = list(pert, as.character(seq_len(n_modules))))
for (group in pert) {
  gene_log2fc <- log2((colMeans(sim$treatment[[group]]) + 0.5) /
                     (control_mean + 0.5))
  for (module in seq_len(n_modules)) {
    genes <- which(sim$gene_info$module == module)
    module_mat[group, module] <- mean(gene_log2fc[genes])
  }
}
write.csv(module_mat, module_log2fc_file)
module_data <- as.data.frame(as.table(module_mat))
names(module_data) <- c('Perturbation', 'Module', 'log2FC')
p_heatmap <- ggplot(module_data,
  aes(x = Module, y = Perturbation, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = '#013e75', mid = '#FAFAFA', high = '#f5b70a',
    midpoint = 0) +
  labs(x = 'Gene module', y = 'Perturbation', fill = 'log2FC') +
  theme_minimal() + theme(panel.grid = element_blank())
ggsave(heatmap_file, p_heatmap, width = 5.5, height = 2.8)
