plot_path<-file.path(work_directory, 'yao_2023', 'report', 'C3_hist_LOEUF')
dir.create(plot_path, recursive = TRUE)

conv_path<-file.path(work_directory, 'yao_2023', 'data')

supp_df<-read_excel(paste0(conv_path, "/raw_data/NIHMS1844621-supplement-Supplementary_Table.xlsx"), sheet="Supplementary Table 5")
supp_df <- supp_df[!is.na(supp_df$gene), ]
supp_df$pLI <- as.numeric(supp_df$pLI)

KO_conventional_effect_size<-read.csv(paste0(conv_path, '/raw_data/GSM6858447_KO_conventional_FRPerturb_effect_sizes.csv.gz'))
KO_matrix<-acast(KO_conventional_effect_size, Downstream_gene ~ Perturbed_gene, value.var = "qvalue", fill = 1)

thres<-0.05/nrow(KO_matrix)
target_MAPK1_name<-names(which(KO_matrix[, 'MAPK1']<thres))

LOEUF<-supp_df$LOEUF
target <- supp_df$LOEUF[supp_df$gene == 'MAPK1'][1]
LOEUF_target<-supp_df$LOEUF[supp_df$gene %in% target_MAPK1_name]

LOEUF_nontarget<-supp_df$LOEUF[!(supp_df$gene %in% target_MAPK1_name)]

df_a<-data.frame(
  score=c(LOEUF_target, LOEUF_nontarget),
  group=c(rep('target', length(LOEUF_target)), rep('non-target', length(LOEUF_nontarget)))
                 )

df_a$group<-as.factor(df_a$group)

p2_density<-ggplot(df_a, aes(x = score, fill = group)) +
  geom_histogram(
    aes(y = after_stat(density)),     # 用密度显示更平滑（可改成Count）
    position = "identity",    # 允许重叠
    alpha = 0.5,              # 半透明
    bins = 50,                # 直方图分箱数
    color = "white"           # 每个柱子加白边区分
  ) +
  scale_fill_manual(values = c("target" = "#f5b70a", "non-target" = "#663B8C"))+
  theme_bw() +
  labs(
    title = "Distribution of LOEUF",
    x = "LOEUF Score",
    y = "Density"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank(),   # 移除外框
    legend.position = "top"
  )
print(p2_density)

ggsave(paste0(plot_path, '/LOEUF.pdf'), p2_density, width = 10, height = 8) # 5.5, 4







