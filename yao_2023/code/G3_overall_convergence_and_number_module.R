library(glue)
library(data.table)
library(ggplot2)
library(ggnewscale)
library(here)

work_directory <- here::here()

batch_name <- "G1_pairwise"
processed_results <- fread(file.path(
  work_directory,
  "yao_2023", "data", "intermediate_data", batch_name, "G2_processed_results.csv"
))

processed_results <- processed_results[, .(comparison, p_value, active_genes)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(
  processed_results,
  processed_results,
  by.x = c("gene1", "gene2"),
  by.y = c("gene2", "gene1")
)
merged_results[, p_max := pmax(p_value.x, p_value.y)]

threshold <- 0.025
merged_results[, is_significant := 0]
merged_results[, is_significant := is_significant + (p_max < threshold)]
merged_results[, is_significant := is_significant + (p_max < threshold / 2450)]
merged_results[, is_significant := factor(is_significant, levels = c(0, 1, 2))]

p_value_dt <- merged_results[gene1 < gene2, .(gene1, gene2, is_significant)]

parse_genes <- function(x) {
  if (is.na(x) || trimws(x) == "") return(character(0))
  parts <- unlist(strsplit(x, "/", fixed = TRUE))
  genes <- trimws(parts)
  unique(genes[genes != ""])
}

count_number_strings <- function(s1, s2, mode = c("union", "intersection")) {
  mode <- match.arg(mode)
  a <- parse_genes(s1)
  b <- parse_genes(s2)
  if (mode == "union") {
    length(union(a, b))
  } else {
    length(intersect(a, b))
  }
}

merged_results[, relevant_gene :=
                 mapply(count_number_strings, active_genes.x, active_genes.y,
                        MoreArgs = list(mode = "intersection"))]
gene_num_dt <- merged_results[gene1 > gene2, .(gene1, gene2, relevant_gene)]

p <- ggplot() +
  geom_point(
    data = gene_num_dt,
    aes(x = gene1, y = gene2, color = relevant_gene, size = relevant_gene),
    alpha = 0.7
  ) +
  scale_color_gradient(
    low = "#EEEEEE",
    high = "#663B8C",
    name = "# of relevant genes",
    guide = guide_colorbar(order = 2)
  ) +
  scale_size_continuous(guide = "none") +
  theme_minimal() +
  labs(x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggnewscale::new_scale_color() +
  geom_point(
    data = p_value_dt,
    aes(x = gene1, y = gene2, color = is_significant),
    size = 3,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("0" = "#EEEEEE", "1" = "#013e75", "2" = "#f5b70a"),
    labels = c("Adj.p<0.05", "p<0.05", "p>0.05"),
    breaks = c("2", "1", "0"),
    name = "Convergence",
    guide = guide_legend(order = 1)
  )

plot(p)

dir.create(file.path(work_directory, "paper_plots", "data"), recursive = TRUE)
plot_path <- file.path(
  work_directory,
  "paper_plots", "data", paste0("G3_", batch_name, "_overall_significance.pdf")
)
ggsave(filename = plot_path, plot = p, width = 9, height = 6 * 9 / 8, bg = "transparent")

dir.create(file.path(work_directory, "yao_2023", "report"), recursive = TRUE)
plot_path <- file.path(
  work_directory,
  "yao_2023", "report", paste0("G3_", batch_name, "_overall_significance.pdf")
)
ggsave(filename = plot_path, plot = p, width = 9, height = 6 * 9 / 8, bg = "transparent")
