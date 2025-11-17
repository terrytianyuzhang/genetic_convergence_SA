result_path<-file.path(work_directory, 'converge_by_module', 'data', 'intermediate_data', 'XConTest')
plot_path<-file.path(work_directory, 'converge_by_module', 'report', 'XConTest_result')
dir.create(plot_path, recursive = TRUE)

amazing_scatterplot_asym<-function(all_names, data_path, store_path, png_name, scale_divider, grouping_vector, active_gene_num, count_method, p_thres=0.05){
  # method is "shared" or "refit"
  all_names<<-all_names
  data_path<<-data_path
  scale_divider<<-scale_divider
  p_thres<<-p_thres
  grouping_vector<<-grouping_vector
  
  num_all_eGene<-length(all_names)
  df<-data.frame(
    eGene_1=rep(all_names, each=length(all_names)),
    eGene_2=rep(all_names, times=length(all_names)),
    Convergent=rep(NA, length(all_names)^2),
    TestStat=rep(NA, length(all_names)^2),
    p_value=rep(0, length(all_names)^2)
  )
  for(i in 1:nrow(df)){
    gene_1<-df$eGene_1[i]
    gene_2<-df$eGene_2[i]
    if(gene_1==gene_2){
      df$TestStat[i]<-0
    }else{
      result_name<-paste0(data_path, '/', gene_1, '_', gene_2, '_', 'one', '.Rdata')
      if(!file.exists(result_name)){
        result_name<-paste0(data_path, '/', gene_2, '_', gene_1, '_', 'one', '.Rdata')
      }
      load(result_name)
      result1<-result
      p_val_1<-result$p_value
      stat_1<-result$test_statistic
      
      result_name<-paste0(data_path, '/', gene_1, '_', gene_2, '_', 'two', '.Rdata')
      if(!file.exists(result_name)){
        result_name<-paste0(data_path, '/', gene_2, '_', gene_1, '_', 'two', '.Rdata')
      }
      load(result_name)
      result2<-result
      p_val_2<-result$p_value
      stat_2<-result$test_statistic
      fold_num<-length(result$fold_data)
      feature_collected_1<-collect_active_features(result1, 
                                                   group = grouping_vector, 
                                                   group_threshold = active_gene_num) 
      feature_collected_2<-collect_active_features(result2, 
                                                   group = grouping_vector, 
                                                   group_threshold = active_gene_num) 
      if(count_method=='union'){
        df$module_num[i] <- length(union(feature_collected_1$active_groups, feature_collected_2$active_groups))
      }else if(count_method=='intersect'){
        df$module_num[i] <- length(intersect(feature_collected_1$active_groups, feature_collected_2$active_groups))
      }else{
        print('Please choose count_method from union or intersect!')
      }
      if(which(all_names==gene_1)>which(all_names==gene_2)){ # bottom right -> number of modules
        # df$TestStat[i]<-max(abs(stat_1), abs(stat_2))
        df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
        df$p_value[i]<-max(p_val_1, p_val_2)
        df$PlotValue[i]<-df$module_num[i]
      }else{ # which(all_names==gene_1)<which(all_names==gene_2) # upper left -> convergent
        # df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
        df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
        df$p_value[i]<-max(p_val_1, p_val_2)
        df$PlotValue[i] <- df$p_value[i] 
      }
      if(df$p_value[i]<=p_thres/(0.5*num_all_eGene*(num_all_eGene-1))){
        df$Convergent[i]<-"Adj.p<0.05"
      }else if(df$p_value[i]>=p_thres){
        df$Convergent[i]<-"p>0.05"
      }else if(df$p_value[i]<p_thres && df$p_value[i]>p_thres/(0.5*num_all_eGene*(num_all_eGene-1))){
        df$Convergent[i]<-"p<0.05"
      }
    } 
  }
  cat('The biggest test statistic is ', max(df$TestStat))
  df$eGene_1 <- factor(df$eGene_1, levels = all_names)
  df$eGene_2 <- factor(df$eGene_2, levels = all_names)
  df_no_na<-df[!is.na(df$Convergent), ]
  
  df_upper <- df_no_na[match(df_no_na$eGene_1, all_names) < match(df_no_na$eGene_2, all_names), ]
  df_upper$Convergent <- factor(df_upper$Convergent,
                                levels = c("Adj.p<0.05", "p>0.05", "p<0.05"))
  df_lower <- df_no_na[match(df_no_na$eGene_1, all_names) > match(df_no_na$eGene_2, all_names), ]
  
  p <- ggplot() +
    # 左上角：按 Convergent 着色
    geom_point(
      data = df_upper,
      aes(x = eGene_1, y = eGene_2, 
          # size = TestStat, 
          color = Convergent),
      size = 8
    ) +
    scale_color_manual(
      name = "Convergence",
      values = c("Adj.p<0.05"="#f5b70a",
                 "p>0.05"="#EEEEEE",
                 "p<0.05"="#013e75")
    ) +
    scale_size_continuous(name = "Test Statistic") +
    
    new_scale_color() +
    
    # 右下角：固定大小
    geom_point(
      data = df_lower,
      aes(
        x = eGene_1,
        y = eGene_2,
        size = module_num  # 按 module_num 决定点大小
      ),
      alpha=0.7,
      shape = 16,          # 圆形
      color = "#663B8C"     # 圆圈颜色
      # color = NULL      # 边框颜色
    ) +
    scale_size_continuous(
      name = "Module Number",
      range = c(0, 10)     # 控制点大小范围
    ) +
    scale_x_discrete(limits = all_names) +  
    scale_y_discrete(limits = all_names) +
    labs(
      title = "Gene Pair Significance",
      x = "Perturbation 1", y = "Perturbation 2"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 14)
    )
  
  print(p)
  
  return(p)
}

# amazing_scatterplot_asym_lasso<-function(all_names, data_path, store_path, png_name, scale_divider, count_method, p_thres=0.05){
#   # method is "shared" or "refit"
#   all_names<<-all_names
#   data_path<<-data_path
#   scale_divider<<-scale_divider
#   p_thres<<-p_thres
#   grouping_vector<<-grouping_vector
#   
#   num_all_eGene<-length(all_names)
#   df<-data.frame(
#     eGene_1=rep(all_names, each=length(all_names)),
#     eGene_2=rep(all_names, times=length(all_names)),
#     Convergent=rep(NA, length(all_names)^2),
#     TestStat=rep(NA, length(all_names)^2),
#     p_value=rep(0, length(all_names)^2)
#   )
#   for(i in 1:nrow(df)){
#     gene_1<-df$eGene_1[i]
#     gene_2<-df$eGene_2[i]
#     if(gene_1==gene_2){
#       df$TestStat[i]<-0
#     }else{
#       result_name<-paste0(data_path, gene_1, '_', gene_2, '_', 'one_lasso', '.Rdata')
#       if(!file.exists(result_name)){
#         result_name<-paste0(data_path, gene_2, '_', gene_1, '_', 'one_lasso', '.Rdata')
#       }
#       load(result_name)
#       result1<-result
#       p_val_1<-result$p_value
#       stat_1<-result$test_statistic
#       
#       result_name<-paste0(data_path, gene_1, '_', gene_2, '_', 'two_lasso', '.Rdata')
#       if(!file.exists(result_name)){
#         result_name<-paste0(data_path, gene_2, '_', gene_1, '_', 'two_lasso', '.Rdata')
#       }
#       load(result_name)
#       result2<-result
#       p_val_2<-result$p_value
#       stat_2<-result$test_statistic
#       fold_num<-length(result$fold_data)
#       genes_1<-c()
#       genes_2<-c()
#       for(j in 1:fold_num){
#         genes_1<-union(genes_1, names(result1$fold_data[[j]]$final_beta[result1$fold_data[[j]]$final_beta!=0]))
#       }
#       for(j in 1:fold_num){
#         genes_2<-union(genes_2, names(result2$fold_data[[j]]$final_beta[result2$fold_data[[j]]$final_beta!=0]))
#       }
#       if(count_method=='union'){
#         df$genes_len[i] <- length(union(genes_1, genes_2))
#       }else if(count_method=='intersect'){
#         df$genes_len[i] <- length(intersect(genes_1, genes_2))
#       }else{
#         print('Please choose count_method from union or intersect!')
#       }
#       if(which(all_names==gene_1)>which(all_names==gene_2)){ # bottom right -> number of modules
#         # df$TestStat[i]<-max(abs(stat_1), abs(stat_2))
#         df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
#         df$p_value[i]<-max(p_val_1, p_val_2)
#         df$PlotValue[i]<-df$genes_len[i]
#       }else{ # which(all_names==gene_1)<which(all_names==gene_2) # upper left -> convergent
#         # df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
#         df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
#         df$p_value[i]<-max(p_val_1, p_val_2)
#         df$PlotValue[i] <- df$p_value[i] 
#       }
#       if(df$p_value[i]<=p_thres/(0.5*num_all_eGene*(num_all_eGene-1))){
#         df$Convergent[i]<-"Adj.p<0.05"
#       }else if(df$p_value[i]>=p_thres){
#         df$Convergent[i]<-"p>0.05"
#       }else if(df$p_value[i]<p_thres && df$p_value[i]>p_thres/(0.5*num_all_eGene*(num_all_eGene-1))){
#         df$Convergent[i]<-"p<0.05"
#       }
#     } 
#   }
#   cat('The biggest test statistic is ', max(df$TestStat))
#   df$eGene_1 <- factor(df$eGene_1, levels = all_names)
#   df$eGene_2 <- factor(df$eGene_2, levels = all_names)
#   df_no_na<-df[!is.na(df$Convergent), ]
#   
#   df_upper <- df_no_na[match(df_no_na$eGene_1, all_names) < match(df_no_na$eGene_2, all_names), ]
#   df_upper$Convergent <- factor(df_upper$Convergent,
#                                 levels = c("Adj.p<0.05", "p>0.05", "p<0.05"))
#   df_lower <- df_no_na[match(df_no_na$eGene_1, all_names) > match(df_no_na$eGene_2, all_names), ]
#   
#   p <- ggplot() +
#     # 左上角：按 Convergent 着色
#     geom_point(
#       data = df_upper,
#       aes(x = eGene_1, y = eGene_2, 
#           # size = TestStat, 
#           color = Convergent),
#       size = 8
#     ) +
#     scale_color_manual(
#       name = "Convergence",
#       values = c("Adj.p<0.05"="#f5b70a",
#                  "p>0.05"="#EEEEEE",
#                  "p<0.05"="#013e75")
#     ) +
#     # scale_size_continuous(name = "Test Statistic") +
#     
#     new_scale_color() +
#     
#     # 右下角：固定大小
#     geom_point(data = df_lower,
#                aes(x = eGene_1, y = eGene_2,
#                    size = genes_len,
#                    color = genes_len),
#                alpha=0.7,
#                show.legend = c(size = FALSE, color = TRUE)) +  # 关闭 size 图例 
#     scale_color_gradient(name = "Genes Count",
#                          low = "#f3e5f5", high = "#6a0dad") +
#     scale_x_discrete(limits = all_names) +  
#     scale_y_discrete(limits = all_names) +
#     labs(
#       title = "Gene Pair Significance",
#       x = "Perturbation 1", y = "Perturbation 2"
#     ) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(size = 14, face = "bold"),
#       axis.title.x = element_text(size = 16),
#       axis.title.y = element_text(size = 16),
#       axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
#       axis.text.y = element_text(size = 12),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 14)
#     )
#   
#   print(p)
#   
#   return(p)
# }

.feature_names <- function(n_modules, module_size, prefix = "feature") {
  mods <- rep(seq_len(n_modules), each = module_size)
  ords <- rep(seq_len(module_size), times = n_modules)
  paste0(prefix, "_", mods, "_", ords)
}

n_samples<-200
n_modules<-30
module_sz<-50
cn<-.feature_names(n_modules, module_sz, prefix = 'feature')
grouping_vector <- rep(1:n_modules, each = module_sz)
names(grouping_vector) <- cn


scale<-638
all_names<-c('A', "B", 'C', 'D', 'E', 'F')
data_path<-result_path
store_path<-plot_path
png_name<-"artificial_NT_3color"
active_gene_num<-5

count_method<-'intersect'

p1<-amazing_scatterplot_asym(all_names, data_path, store_path, png_name, scale, grouping_vector, active_gene_num, count_method)
ggsave(paste0(plot_path, "/artificial_NT_3_color_raw_grlasso_1se_200_400_intersect.pdf"), p1, width = 5.5, height = 4)





