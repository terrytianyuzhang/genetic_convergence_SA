plot_path<-file.path(work_directory, 'simulation', 'report')
dir.create(plot_path, recursive = TRUE)

data_path<-file.path(work_directory, 'simulation', 'data', 'final_data')

CI<-read.csv(paste0(data_path, '/CI.csv'))
CI$gene <- rownames(CI)
CI$gene <- factor(CI$gene, levels = CI$gene)  

label_genes <- CI$gene[seq(1, nrow(CI), by = 50)]
label_genes <- label_genes[-1]
CI_sub <- CI[CI$gene %in% label_genes, ]
CI_sub$gene <- droplevels(CI_sub$gene)  # 去除无关 levels

dodge <- position_dodge(width = 0.2)

p1 <- ggplot() +
  geom_point(data = CI_sub, aes(x = gene, y = diff_mean, color = "Generated Data"), 
             size = 2, position = dodge) +
  # geom_errorbar(data = CI_sub, 
  #               aes(x = gene, ymin = diff_mean - 2 * diff_sd, ymax = diff_mean + 2 * diff_sd, color = "Generated Data"), 
  #               width = 0.5, position = dodge) +
  geom_point(data = CI_sub, aes(x = gene, y = observed_mean, color = "Observed Data"), 
             size = 2, position = dodge) +
  # geom_errorbar(data = CI_sub, 
  #               aes(x = gene, ymin = observed_mean - 2 * observed_sd, ymax = observed_mean + 2 * observed_sd, color = "Observed Data"), 
  #               width = 0.5, position = dodge) +
  scale_x_discrete(breaks = label_genes) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Gene Index", y = "Value", color = "Type", 
       title = "Mean of Observed Data and Generated data") +
  scale_color_manual(values = c("Generated Data" = "#013e75", "Observed Data" = "#f5b70a")) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),    # 图例标题字号
    legend.text = element_text(size = 14)      # 图例内容字号
  )

print(p1)

################################## QQ ########################################
load(paste0(data_path, '/test_stat_one_null_sample_200.Rdata'))
n <- length(test_stat_one)
observed <- sort(test_stat_one)
theoretical <- qnorm(ppoints(n))  # 正态分布的理论分位数

df_qq <- data.frame(theoretical = theoretical, observed = observed)

p2<-ggplot(df_qq, aes(x = theoretical, y = observed)) +
  geom_point(color='#f5b70a') +
  geom_abline(slope = 1, intercept = 0, color = "#222222") +  # 45度参考线
  labs(title = "QQ Plot for Test Stat Under the Null",
       x = "Theoretical Quantiles (Normal)",
       y = "Observed Quantiles") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),       # 标题字号和加粗
    axis.title.x = element_text(size = 16),                    # x轴标签字号
    axis.title.y = element_text(size = 16),                    # y轴标签字号
    axis.text.x = element_text(size = 12),                     # x轴刻度文字字号
    axis.text.y = element_text(size = 12)                      # y轴刻度文字字号
  )
print(p2)
################################## NULL ########################################
thres<-0.05
sample_use<-c(100, 200, 500, 800, 1000)
type_I_max<-numeric(length(sample_use))
type_I_min<-numeric(length(sample_use))
type_1<-numeric(length(sample_use))
type_2<-numeric(length(sample_use))
type_ave<-numeric(length(sample_use))
for(i in 1:length(sample_use)){
  n<-sample_use[i]
  load(paste0(data_path, '/p_vec_one_null_sample_', n, '.Rdata'))
  load(paste0(data_path, '/p_vec_two_null_sample_', n, '.Rdata'))
  type_1[i]<-sum(p_vec_one<thres)/length(p_vec_one)
  type_2[i]<-sum(p_vec_two<thres)/length(p_vec_two)
  p_min<-pmin(p_vec_one, p_vec_two)
  p_max<-pmax(p_vec_one, p_vec_two)
  p_ave<-(p_vec_one+p_vec_two)/2
  type_I_max[i]<-sum(p_max<thres)/length(p_max)
  type_I_min[i]<-sum(p_min<thres)/length(p_min)
  type_ave[i]<-sum(p_ave<0.5*thres)/length(p_ave)
}
df <- data.frame(
  index = sample_use,
  value_max = type_I_max,
  value_min = type_I_min,
  value_1 = type_1,
  value_2 = type_2,
  value_ave = type_ave
)

df_long <- df %>%
  pivot_longer(
    cols = c(value_max, value_min, value_1, value_2),
    names_to = "type",
    values_to = "value"
  )

# 2. 设置每种 type 对应的颜色和线型
color_map <- c(
  "value_max" = "#f5b70a",
  "value_min" = "#013e75", 
  "value_1"   = "#D9D9D9",
  "value_2"   = "#4D4D4D" 
)

linetype_map <- c(
  "value_max" = "solid",
  "value_min" = "solid",   # 改成点线
  "value_1"   = "dotdash",  # 改成点划线
  "value_2"   = "dotdash"  # 改成长虚线
)

label_map <- c(
  "value_max" = 'Pooled max',
  "value_min" = 'Pooled min',
  "value_1"   = 'Indv. 1to2',
  "value_2"   = 'Indv. 2to1'
)

shape_map <- c(
  "value_max" = 15,  # 方块
  "value_min" = 18,  # 菱形
  "value_1"   = 16,  # 圆点
  "value_2"   = 17   # 三角
)

p3 <- ggplot(df_long, aes(x = index, y = value,
                          color = type, linetype = type, shape = type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = thres,
             color = "black", linetype = "dashed") +
  scale_color_manual(values = color_map,
                     labels = label_map,
                     name = NULL) +
  scale_linetype_manual(values = linetype_map,
                        labels = label_map,
                        name = NULL) +
  scale_shape_manual(values = shape_map,
                     labels = label_map,
                     name = NULL) +
  labs(
    title = "Type-I Error under Different Sample Sizes under the Null",
    x = "Sample Size",
    y = "Type-I Error"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.position = "top"
  ) +
  guides(
    color = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

print(p3)


################################## ALTER ########################################
sample_use<-c(10, 30, 50, 80, 100, 300)
power_max<-numeric(length(sample_use))
power_min<-numeric(length(sample_use))
power_ave<-numeric(length(sample_use))
for(i in 1:length(sample_use)){
  n<-sample_use[i]
  load(paste0(data_path, '/p_vec_one_alter_sample_', n, '.Rdata'))
  load(paste0(data_path, '/p_vec_two_alter_sample_', n, '.Rdata'))
  p_min<-pmin(p_vec_one, p_vec_two)
  p_max<-pmax(p_vec_one, p_vec_two)
  power_max[i]<-sum(p_max<thres)/length(p_max)
  power_min[i]<-sum(p_min<thres)/length(p_min)
  p_ave<-(p_vec_one+p_vec_two)/2
  power_ave[i]<-sum(p_ave<0.5*thres)/length(p_ave)
}

df <- data.frame(
  index = sample_use,
  value_max = power_max,
  value_min = power_min,
  value_ave = power_ave
)

df_long <- df %>%
  pivot_longer(
    cols = c(value_max, value_min, value_ave),
    names_to = "type",
    values_to = "value"
  )

color_map <- c(
  "value_max" = "#f5b70a",
  "value_min" = "#013e75"
  # "value_ave"   = "#6F42C0"
)

linetype_map <- c(
  "value_max" = "solid",
  "value_min" = "solid"   # 改成点线
  # "value_ave"   = "dotdash"  # 改成点划线
)

label_map <- c(
  "value_max" = 'Pooled max',
  "value_min" = 'Pooled min'
  #"value_ave"   = expression(Power['ave'])
)

shape_map <- c(
  "value_max" = 15,  # 方块
  "value_min" = 18  # 菱形
  # "value_ave"   = 16  # 圆点
)

p4 <- ggplot(df_long, aes(x = index, y = value,
                          color = type, linetype = type, shape = type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = color_map,
                     labels = label_map,
                     name = NULL) +
  scale_linetype_manual(values = linetype_map,
                        labels = label_map,
                        name = NULL) +
  scale_shape_manual(values = shape_map,
                     labels = label_map,
                     name = NULL) +
  labs(
    title = "Power under Different Sample Sizes under the Alternative",
    x = "Sample Size",
    y = "Power"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.position = "top"
  ) +
  guides(
    color = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

print(p4)

ggsave(paste0(plot_path, "/sim_a.pdf"), p1, width = 6.5, height = 4.5)
ggsave(paste0(plot_path, "/sim_b.pdf"), p2, width = 6.5, height = 4.5)
ggsave(paste0(plot_path, "/sim_c.pdf"), p3, width = 6.5, height = 4.5)
ggsave(paste0(plot_path, "/sim_d.pdf"), p4, width = 6.5, height = 4.5)

