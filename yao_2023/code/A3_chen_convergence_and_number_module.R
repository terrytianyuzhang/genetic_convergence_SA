
out_dir <- file.path(work_directory, "yao_2023", "data", "intermediate_data", "A1_chen_method")
rds_path <- file.path(out_dir, "processed_data.rds")
processed_results <- readRDS(file = rds_path)

processed_results <- processed_results[, .(comparison, p_value, active_group)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(processed_results, processed_results, 
                        by.x = c("gene1", "gene2"), by.y = c("gene2", "gene1"))
# merged_results[, p_min := pmin(p_value.x, p_value.y)]
merged_results[, p_max := pmax(p_value.x, p_value.y)]

threshold <- 0.025
merged_results[, is_significant := 0]
merged_results[, is_significant := is_significant + 2*(p_max < threshold)]
merged_results[, is_significant := factor(is_significant, levels = c(0, 1, 2))]

p_value_dt <- merged_results[gene1 < gene2, .(gene1, gene2, is_significant)]


parse_nums <- function(x) {
  if (is.na(x) || trimws(x) == "") return(integer(0))
  parts <- unlist(strsplit(x, "/", fixed = TRUE))
  nums  <- suppressWarnings(as.integer(trimws(parts)))
  nums  <- nums[!is.na(nums)]
  unique(nums)
}

count_number_strings <- function(s1, s2, mode = c("union", "intersection")) {
  mode <- match.arg(mode)
  a <- parse_nums(s1)
  b <- parse_nums(s2)
  if (mode == "union") {
    length(union(a, b))
  } else {
    length(intersect(a, b))
  }
}

merged_results[, relevant_module :=
                 mapply(count_number_strings, active_group.x, active_group.y,
                        MoreArgs = list(mode = "intersection"))]
module_num_dt <- merged_results[gene1 > gene2, .(gene1, gene2, relevant_module)]

library(ggplot2)
library(ggnewscale)

p <- ggplot() +
  # modules: continuous color (will be 2nd legend)
  geom_point(data = module_num_dt,
             aes(x = gene1, y = gene2, color = relevant_module, size = relevant_module),
             alpha = 0.7) +
  scale_color_gradient(low = "#EEEEEE", high = "#663B8C",
                       name = "# of relevant modules",
                       guide = guide_colorbar(order = 2)) +
  scale_size_continuous(guide = "none") +
  theme_minimal() +
  labs(x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  
  ggnewscale::new_scale_color() +
  
  # significance: discrete colors (will be 1st legend, on top)
  geom_point(data = p_value_dt,
             aes(x = gene1, y = gene2, color = is_significant),
             size = 3, alpha = 0.7) +
  scale_color_manual(
    values = c("0" = "#EEEEEE", "2" = "#f5b70a"),
    labels = c("Yes", "No"),
    breaks = c("2", "0"),
    name   = "Convergence",
    guide  = guide_legend(order = 1)
  )

plot(p)

dir.create(file.path(work_directory, "paper_plots", "data"), 
           recursive = TRUE)
ggsave(p, file = here(work_directory, "paper_plots/data/A3_chen_convergence.pdf"),
       width = 9, height = 6*9/8, bg = "transparent")

ggsave(p, file = here(work_directory, "yao_2023/report/A3_chen_convergence.pdf"),
       width = 9, height = 6*9/8, bg = "transparent")

summary(p_value_dt$is_significant)
