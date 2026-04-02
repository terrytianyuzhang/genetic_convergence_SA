
batch_name <- '51_pairwise'
processed_results <- fread(glue(work_directory, '/yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))

processed_results <- processed_results[, .(comparison, p_value, active_group)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

threshold <- 0.025
processed_results[, is_significant := 0]
processed_results[, is_significant := is_significant + (p_value < threshold)]
processed_results[, is_significant := is_significant + (p_value < threshold / 2450)]
processed_results[, is_significant := factor(is_significant, levels = c(0, 1, 2))]

parse_nums <- function(x) {
  if (is.na(x) || trimws(x) == "") return(integer(0))
  parts <- unlist(strsplit(x, "/", fixed = TRUE))
  nums <- suppressWarnings(as.integer(trimws(parts)))
  nums <- nums[!is.na(nums)]
  unique(nums)
}

processed_results[, relevant_module := vapply(active_group, function(x) length(parse_nums(x)), integer(1))]

gene_levels <- sort(unique(c(processed_results$gene1, processed_results$gene2)))

base_theme <- theme_minimal() +
  labs(x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_directional_significance <- function(plot_dt) {
  ggplot(plot_dt, aes(x = gene1, y = gene2, color = is_significant)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_x_discrete(limits = gene_levels) +
    scale_y_discrete(limits = rev(gene_levels)) +
    scale_color_manual(
      values = c("0" = "#EEEEEE", "1" = "#013e75", "2" = "#f5b70a"),
      labels = c("Adj.p<0.05", "p<0.05", "p>0.05"),
      breaks = c("2", "1", "0"),
      name = "Convergence"
    ) +
    base_theme
}

plot_directional_count <- function(plot_dt) {
  ggplot(plot_dt, aes(x = gene1, y = gene2, color = relevant_module, size = relevant_module)) +
    geom_point(alpha = 0.7) +
    scale_x_discrete(limits = gene_levels) +
    scale_y_discrete(limits = rev(gene_levels)) +
    scale_color_gradient(
      low = "#EEEEEE",
      high = "#663B8C",
      name = "# of relevant modules"
    ) +
    scale_size_continuous(guide = "none") +
    base_theme
}

significance_plot <- plot_directional_significance(processed_results)
count_plot <- plot_directional_count(processed_results)

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
