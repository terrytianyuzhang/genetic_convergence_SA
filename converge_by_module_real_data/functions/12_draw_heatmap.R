rm(list = ls())

# INPUT
data_path <- file.path('..', 'data', 'raw')
data_file <- file.path(data_path, '10_splatter_data.Rdata')
# OUTPUT
result_path <- file.path('..', 'data', 'final', '12_heatmap')
module_log2fc_file <- file.path(result_path, '12_module_mean_log2fc.csv')

plot_path <- file.path('..', 'figures', '12_splatter_pilot')
heatmap_file <- file.path(plot_path, '12_splatter_pilot_module_heatmap.pdf')
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)


############################ Plot observed module log2 fold changes ###########
# Match 06: remove the group library-size offset measured before perturbation.
load(data_file)
pert <- names(sim$treatment)
n_modules <- length(unique(sim$gene_info$module))
original_library_size <- rowSums(sim$original_counts)
stopifnot(all(original_library_size > 0))
mean_library_size <- tapply(original_library_size, sim$cell_info$group, mean)
library_offset <- log2(mean_library_size[pert] / mean_library_size['NT'])
control_mean <- colMeans(sim$control)
module_mat <- matrix(0, length(pert), n_modules,
  dimnames = list(pert, as.character(seq_len(n_modules))))
for (group in pert) {
  gene_log2fc <- log2((colMeans(sim$treatment[[group]]) + 0.5) /
                     (control_mean + 0.5))
  for (module in seq_len(n_modules)) {
    genes <- which(sim$gene_info$module == module)
    module_mat[group, module] <- mean(gene_log2fc[genes]) - library_offset[group]
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
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank(), legend.position = 'top') +
  guides(fill = guide_colorbar(direction = 'horizontal',
    title.position = 'top', title.hjust = 0.5,
    barwidth = 20, barheight = 0.5))
ggsave(heatmap_file, p_heatmap, width = 9, height = 4)
