result_path<-file.path(work_directory, 'converge_by_module', 'data', 'intermediate_data', 'Dset')
plot_path<-file.path(work_directory, 'converge_by_module', 'report', 'Dset_result')
dir.create(plot_path, recursive = TRUE)
load(paste0(result_path, '/p_val_mat.Rdata'))


pert<-c('A', 'B', 'C', 'D', 'E', 'F')

pair_num<-15
n_modules<-ncol(p_val_mat)
threshold<-0.05/(pair_num*n_modules)

df<-data.frame(pert_1=rep(pert, each=length(pert)),
               pert_2=rep(pert, length(pert)),
               convergence=rep(NA, length(pert)^2),
               module_num=rep(NA, length(pert)^2))

for(i in 1:nrow(df)){
  pert_name_1<-df$pert_1[i]
  pert_name_2<-df$pert_2[i]
  if(pert_name_1==pert_name_2){
    next
  }else{
    p_mat_sub<-p_val_mat[c(pert_name_1, pert_name_2), ]
    count<-0
    for(j in 1:ncol(p_mat_sub)){
      if(p_mat_sub[pert_name_1, j]<threshold & p_mat_sub[pert_name_2, j]<threshold){
        count<-count+1
      }
    }
    df$module_num[i]<-count
    if(count>0){
      df$convergence[i]<-'Yes'
    }else{
      df$convergence[i]<-'No'
    }
  }
}
df$convergence<-factor(df$convergence)
df$pert_1 <- factor(df$pert_1, levels = pert)
df$pert_2 <- factor(df$pert_2, levels = pert)
df_no_diag <- df[df$pert_1 != df$pert_2, ]

df_upper <- df_no_diag[match(df_no_diag$pert_1, pert) < match(df_no_diag$pert_2, pert), ]
df_lower <- df_no_diag[match(df_no_diag$pert_1, pert) > match(df_no_diag$pert_2, pert), ]

p <- ggplot() +
  # 左上角：按 convergence 着色
  geom_point(
    data = df_upper,
    aes(x = pert_1, y = pert_2, color = convergence),
    size = 6
  ) +
  scale_color_manual(
    name = "Convergence",
    values = c("Yes" = "#f5b70a", "No" = "#EEEEEE")
  ) +
  
  new_scale_color() +
  
  # 右下角：点大小表示 module_num 
  geom_point(
    data = df_lower,
    aes(x = pert_1, y = pert_2, size = module_num),
    alpha = 0.7,
    shape = 16,
    color = "#6a0dad"
  ) +
  scale_size_continuous(
    name = "Module Number",
    range = c(1, 10)   # 控制点大小范围
  ) +
  
  scale_x_discrete(limits = pert) +
  scale_y_discrete(limits = pert) +
  labs(
    title = "Perturbation Pair Convergence",
    x = "Perturbation 1", y = "Perturbation 2"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p)

ggsave(paste0(plot_path, '/Chen_revise.pdf'), p, width = 5.5, height = 4)


