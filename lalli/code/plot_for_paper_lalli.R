plot_path<-file.path(work_directory, 'lalli', 'report')
dir.create(plot_path, recursive = TRUE)

data_path<-file.path(work_directory, 'lalli', 'data', 'intermediate_data')

amazing_scatterplot_asym<-function(all_names, method, data_path, postfix, store_path, png_name, scale_divider, count_method, p_thres=0.05){
  # method is "shared" or "refit"
  all_names<<-all_names
  method<<-method
  data_path<<-data_path
  postfix<<-postfix
  scale_divider<<-scale_divider
  p_thres<<-p_thres
  
  all_names <- gsub(" late", "", all_names)
  
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
      result_name<-paste0(data_path, '/', gene_1, ' late', '-', gene_2, ' late', '-', method, '-', postfix, '.Rdata')
      if(!file.exists(result_name)){
        result_name<-paste0(data_path, '/', gene_2, ' late', '-', gene_1, ' late', '-', method, '-', postfix, '.Rdata')
      }
      load(result_name)
      p_val_1<-result$p_value
      stat_1<-result$test_statistic
      result1<-result
      result_name<-paste0(data_path, '/', gene_1, ' late', '-', gene_2, ' late', '-', method, '-', '2big-', postfix, '.Rdata')
      if(!file.exists(result_name)){
        result_name<-paste0(data_path, '/', gene_2, ' late', '-', gene_1, ' late', '-', method, '-', '2big-', postfix, '.Rdata')
      }
      load(result_name)
      result2<-result
      p_val_2<-result$p_value
      stat_2<-result$test_statistic
      genes_1<-c()
      genes_2<-c()
      fold_num<-length(result1$split_data)
      for(j in 1:fold_num){
        genes_1<-union(genes_1, names(result1$split_data[[j]]$final_beta[result1$split_data[[j]]$final_beta!=0]))
      }
      for(j in 1:fold_num){
        genes_2<-union(genes_2, names(result2$split_data[[j]]$final_beta[result2$split_data[[j]]$final_beta!=0]))
      }
      if(count_method=='union'){
        df$genes_len[i] <- length(union(genes_1, genes_2))
      }else if(count_method=='intersect'){
        df$genes_len[i] <- length(intersect(genes_1, genes_2))
      }else{
        print('Please choose count_method from union or intersect!')
      }
      ################### We count numbers of genes for the result with a bigger p value ###########
      if(which(all_names==gene_1)>which(all_names==gene_2)){ # bottom right -> number of genes
        df$TestStat[i]<-min(abs(stat_1), abs(stat_2))
        df$p_value[i]<-max(p_val_1, p_val_2)
        df$PlotValue[i] <- df$genes_len[i]
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
  # print(levels(df_upper$eGene_1))
  # print(levels(df_lower$eGene_1))
  # print(levels(df_upper$eGene_2))
  # print(levels(df_lower$eGene_2))
  
  p <- ggplot() +
    geom_point(data = df_upper,
               aes(x = eGene_1, y = eGene_2,
                   # size = TestStat, 
                   color = Convergent),
               size = 6) +
    scale_color_manual(name = "Convergence",
                       values = c("Adj.p<0.05"="#f5b70a",
                                  "p<0.05"="#013e75",
                                  "p>0.05"="#EEEEEE"
                                  ),
                       breaks = c("Adj.p<0.05", "p<0.05", "p>0.05"),
                       guide = guide_legend(order = 1))+   
    new_scale_color() +  
    geom_point(data = df_lower,
               aes(x = eGene_1, y = eGene_2,
                   size = genes_len,
                   color = genes_len),
                    alpha=0.7,
                   show.legend = c(size = FALSE, color = TRUE)) +  # 关闭 size 图例) +
    scale_color_gradient(name = "Genes Count",
                         low = "#f3e5f5", high = "#6a0dad",
                         guide = guide_colorbar(order = 2)) +
    # scale_size_continuous(name = "Genes Count",
    #                       # range = c(0, max(df$genes_len, na.rm = TRUE))  
    #                       range = c(2, 8) # 控制点大小范围
    #                       ) +
    scale_x_discrete(limits = all_names) +
    scale_y_discrete(limits = all_names) +
    labs(title = "Gene Pair Significance",
         x = "Gene 1", y = "Gene 2") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) 
  
  
  print(p)
  return(p)
}

scale<-9.6
all_names<-c("CHD2 late", "DYRK1A late", "HDAC5 late", "PTEN late", "CHD8 late", "ADNP late", "POGZ late",
         "RELN late", "MECP2 late", "MYT1L late", "CTNND2 late", "SETD5 late", "ASH1L late", "ARID1B late")
all_names<-c("CHD2 late", "PTEN late", "CHD8 late", "ADNP late", "RELN late", "MECP2 late", "CTNND2 late", "SETD5 late")
method<-"refit-lognorm-Lallilate"
postfix<-"seq_result_indivdrop"
store_path<-'/Users/erganshang/Desktop/diffusion_and_protein/diff_conv/eff_emb/plot/for_paper/'
png_name<-"lognorm_indivdrop_refit_14_3color"

count_method<-'intersect'

p1<-amazing_scatterplot_asym(all_names, method, data_path, postfix, store_path, png_name, scale, count_method)
ggsave(paste0(plot_path, "/lalli_dot_plot_intersect.pdf"), p1, width = 5.5, height = 4)
################################################################################################
load(paste0(data_path, "/normal_drop_t_lallilate.Rdata"))
all_genes<-colnames(normal_drop_t_lallilate[[1]])

load(paste0(data_path, "/ADNP late-PTEN late-refit-lognorm-Lallilate-2big-seq_result_indivdrop.Rdata"))

gene_name<-c()
for(i in 1:10){
  gene_name<-union(gene_name, names(result$split_data[[i]]$final_beta[result$split_data[[i]]$final_beta!=0]))
}

load(paste0(data_path, "/ADNP late-PTEN late-refit-lognorm-Lallilate-seq_result_indivdrop.Rdata"))

for(i in 1:10){
  gene_name<-union(gene_name, names(result$split_data[[i]]$final_beta[result$split_data[[i]]$final_beta!=0]))
}

GO_result <- enrichGO(gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP", # try "BP", "MF", "CC" instead of "ALL"
                      pvalueCutoff = 0.05)
GO_result_df <- as.data.frame(GO_result@result)

n_use<-15
df_top <- GO_result_df %>%
  arrange(pvalue) %>%
  slice_head(n = n_use)

# 画图
p <- ggplot(df_top, aes(x = reorder(Description, -pvalue),
                          y = Count,
                          fill = -log10(pvalue))) +
  geom_col() +
  coord_flip() +   # 横着的柱状图
  scale_fill_gradient(
    low = "grey90",
    high = "#663B8C"
  ) +
  labs(x = "", y = "Count", title = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11),
    legend.position = "top",          # <--- 把图例放在上面
    legend.direction = "horizontal"   # <--- 横向色条
  )+
  guides(
    fill = guide_colorbar(
      title.position = "top",   # 标题在色条上方
      title.hjust = 0.5,        # 标题居中
      barwidth = 20,            # 色条变长 (可调大)
      barheight = 0.5             # 色条高度 (可适当调小)
    )
  )
print(p)

ggsave(filename = paste0(plot_path, "/GO_barplot.pdf"), plot = p, width = 8, height = 4.5, dpi = 300)

write.csv(GO_result_df, file=paste0(plot_path, "/ADNP-PTEN-refit-1_and_2big-log-indivdrop-GO-term.csv"))









