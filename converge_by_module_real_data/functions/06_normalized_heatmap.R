# Run this script from the functions directory after 01_prepare_Xin_ctrl_and_effect.R.
data_file <- file.path(
  "..", "data", "Xin_151_200_400.Rdata"
)
effect_file <- file.path(
  "..", "data", "poisson_effects.Rdata"
)
pair_file <- file.path(
  "..", "data", "convergence_pair_summary.csv"
)
module_pdf <- file.path(
  "..", "figures", "06_normalized_module_heatmap.pdf"
)
combined_pdf <- file.path(
  "..", "figures", "06_normalized_convergence_figure.pdf"
)

load(data_file)
load(effect_file)

groups <- c("A", "B", "C", "D", "E", "F")
stopifnot(
  identical(rownames(program_mat), groups),
  ncol(program_mat) == 30,
  all(c(groups, "NT") %in% sim$cell_info$group)
)

# Calculate library sizes before injecting effects. The group-level offset
# leaves the original counts and saved Poisson results unchanged.
original_library_size <- rowSums(sim$original_counts)
stopifnot(all(original_library_size > 0))
mean_library_size <- tapply(
  original_library_size,
  sim$cell_info$group,
  mean
)
library_offset <- log2(
  mean_library_size[groups] / mean_library_size["NT"]
)
module_log2fc <- sweep(program_mat, 1, library_offset, "-")

module_data <- as.data.frame(as.table(module_log2fc))
module_plot <- ggplot2::ggplot(
  module_data,
  ggplot2::aes(x = Var2, y = Var1, fill = Freq)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#013e75", mid = "#FAFAFA", high = "#f5b70a", midpoint = 0
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Gene module", y = "Perturbation", fill = "log2FC"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    legend.position = "top"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 20,
      barheight = 0.5
    )
  )

ggplot2::ggsave(module_pdf, module_plot, width = 9, height = 4)

# Pairwise convergence panel, matching the intersection-based Figure 3 style.
pair_data <- read.csv(pair_file)
stopifnot(
  nrow(pair_data) == choose(length(groups), 2),
  all(c("pert_1", "pert_2", "Convergent", "module_intersect") %in%
    names(pair_data)),
  all(pair_data$pert_1 %in% groups),
  all(pair_data$pert_2 %in% groups)
)
upper_data <- pair_data
upper_data$Convergent <- factor(
  upper_data$Convergent,
  levels = c("Adj.p<0.05", "p>0.05", "p<0.05")
)
lower_data <- pair_data
lower_data$pert_1 <- pair_data$pert_2
lower_data$pert_2 <- pair_data$pert_1

pair_plot <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = upper_data,
    ggplot2::aes(x = pert_1, y = pert_2, color = Convergent),
    size = 8
  ) +
  ggplot2::scale_color_manual(
    name = "Convergence",
    breaks = c("Adj.p<0.05", "p>0.05"),
    values = c(
      "Adj.p<0.05" = "#f5b70a",
      "p>0.05" = "#EEEEEE",
      "p<0.05" = "#013e75"
    )
  ) +
  ggplot2::geom_point(
    data = lower_data,
    ggplot2::aes(x = pert_1, y = pert_2, size = module_intersect),
    alpha = 0.7,
    shape = 16,
    color = "#663B8C"
  ) +
  ggplot2::scale_size_continuous(
    name = "Module Number", range = c(0, 10),
    breaks = 0:max(pair_data$module_intersect)
  ) +
  ggplot2::scale_x_discrete(limits = groups) +
  ggplot2::scale_y_discrete(limits = groups) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::labs(
    title = "Gene Pair Significance",
    x = "Perturbation 1", y = "Perturbation 2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14, face = "bold"),
    axis.title = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 12)
  )

module_panel <- module_plot +
  ggplot2::scale_x_discrete(breaks = c("1", "10", "20", "30")) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1),
    legend.title = ggplot2::element_text(size = 13),
    legend.text = ggplot2::element_text(size = 12)
  )
combined_plot <- cowplot::plot_grid(
  module_panel, pair_plot,
  nrow = 1,
  rel_widths = c(1.5, 1.5),
  labels = c("a", "b"),
  label_size = 20,
  label_fontface = "bold"
)

ggplot2::ggsave(combined_pdf, combined_plot, width = 13.5, height = 5.05)
