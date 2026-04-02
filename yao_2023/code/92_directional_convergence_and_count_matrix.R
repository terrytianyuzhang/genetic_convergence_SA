
batch_name <- '51_pairwise'
processed_results <- fread(glue(work_directory, '/yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))

processed_results <- processed_results[, .(comparison, p_value, active_group)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

threshold <- 0.025
processed_results[, is_significant := 0]
processed_results[, is_significant := is_significant + (p_value < threshold)]
processed_results[, is_significant := is_significant + (p_value < threshold / 2450)]
processed_results[, is_significant := factor(is_significant, levels = c(0, 1, 2))]

processed_results[, relevant_module := 0L]
valid_active_group <- !is.na(processed_results$active_group) & trimws(processed_results$active_group) != ""
processed_results[valid_active_group,
                  relevant_module := vapply(strsplit(active_group, "/", fixed = TRUE), function(x) {
                    module_ids <- suppressWarnings(as.integer(trimws(x)))
                    length(unique(module_ids[!is.na(module_ids)]))
                  }, integer(1))]

gene_levels <- sort(unique(c(processed_results$gene1, processed_results$gene2)))

significance_plot <- ggplot(processed_results, aes(x = gene1, y = gene2, color = is_significant)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_discrete(limits = gene_levels) +
  scale_y_discrete(limits = gene_levels) +
  scale_color_manual(
    values = c("0" = "#EEEEEE", "1" = "#013e75", "2" = "#f5b70a"),
    labels = c("Adj.p<0.05", "p<0.05", "p>0.05"),
    breaks = c("2", "1", "0"),
    name = "Convergence"
  ) +
  theme_minimal() +
  labs(x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

count_plot <- ggplot(processed_results, aes(x = gene1, y = gene2, color = relevant_module, size = relevant_module)) +
  geom_point(alpha = 0.7) +
  scale_x_discrete(limits = gene_levels) +
  scale_y_discrete(limits = gene_levels) +
  scale_color_gradient(
    low = "#EEEEEE",
    high = "#663B8C",
    name = "# of relevant modules"
  ) +
  scale_size_continuous(guide = "none") +
  theme_minimal() +
  labs(x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dir.create(file.path(work_directory, "paper_plots", "data"), recursive = TRUE)
dir.create(file.path(work_directory, "yao_2023", "report"), recursive = TRUE)

plot_path <- file.path(
  work_directory,
  "paper_plots",
  "data",
  paste0("92_", batch_name, "_directional_significance.pdf")
)
ggsave(filename = plot_path, plot = significance_plot, width = 9, height = 6 * 9 / 8, bg = "transparent")

plot_path <- file.path(
  work_directory,
  "yao_2023",
  "report",
  paste0("92_", batch_name, "_directional_significance.pdf")
)
ggsave(filename = plot_path, plot = significance_plot, width = 9, height = 6 * 9 / 8, bg = "transparent")

plot_path <- file.path(
  work_directory,
  "paper_plots",
  "data",
  paste0("92_", batch_name, "_directional_count.pdf")
)
ggsave(filename = plot_path, plot = count_plot, width = 9, height = 6 * 9 / 8, bg = "transparent")

plot_path <- file.path(
  work_directory,
  "yao_2023",
  "report",
  paste0("92_", batch_name, "_directional_count.pdf")
)
ggsave(filename = plot_path, plot = count_plot, width = 9, height = 6 * 9 / 8, bg = "transparent")
