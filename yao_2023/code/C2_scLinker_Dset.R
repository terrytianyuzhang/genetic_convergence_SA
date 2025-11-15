plot_path<-file.path(work_directory, 'yao_2023', 'report', 'C2_scLinker_Dset')
dir.create(plot_path, recursive = TRUE)

conv_path<-file.path(work_directory, 'yao_2023', 'data')

scLinker_df<-read_excel(paste0(conv_path, "/raw_data/supp_yao.xlsx"), sheet="Supplementary Table 4b")
colnames(scLinker_df)<-scLinker_df[1, ]
scLinker_df<-scLinker_df[-1,]
scLinker_df$Enrichment_score_pvalue <- as.numeric(as.character(scLinker_df$Enrichment_score_pvalue))

direct_df<-readRDS(paste0(conv_path, '/intermediate_data/A1_chen_method/processed_data.rds'))
all_pert<-unique(direct_df$treatment1)

scLinker_pheno_conv<-function(scLinker_df_fun, conv_df_fun, phenotype, knock_way, p_method, pheno_thres){
  print(phenotype)
  print(knock_way)
  # phenotype could be 'Inflammatory bowel disease'
  # knock_way (a set) could be (pos)/(neg) + Knock-down/Knock-out, say (pos) Knock-down
  # p_method could be 'max' or 'min'
  # phenotype<-'Primary biliary cirrhosis'
  # knock_way<-c("(pos) Knock-down", "(pos) Knock-out", "(neg) Knock-down", "(neg) Knock-out")
  # scLinker_df_fun<-scLinker_df
  # conv_df_fun<-conv_df
  scLinker_use<-scLinker_df_fun[scLinker_df$Trait==phenotype, ]
  scLinker_use <- scLinker_use[
    sapply(scLinker_use$Geneset, function(x) any(endsWith(x, knock_way))),
  ]
  # ACTL6A
  scLinker_use$base_geneset <- sub(" \\([^)]*\\).*", "", scLinker_use$Geneset)
  sel_idx <- unlist(tapply(seq_len(nrow(scLinker_use)), scLinker_use$base_geneset, function(idxs) {
    idxs[which.min(scLinker_use$Enrichment_score_pvalue[idxs])]
  }))
  scLinker_use <- scLinker_use[sel_idx, , drop = FALSE]
  scLinker_use$Geneset <- scLinker_use$base_geneset
  scLinker_use$base_geneset <- NULL
  
  # print(class(scLinker_use$Enrichment_score_pvalue))
  # print(any(is.na(scLinker_use$Enrichment_score_pvalue)))
  
  print(all(all_pert %in% scLinker_use$Geneset))
  
  all_gene_focus<-intersect(all_pert, scLinker_use$Geneset)
  conv_thres<-0.05/(0.5*length(all_gene_focus)*(length(all_gene_focus)-1))
  
  conv_df_fun$conv_p<-rep(NA, nrow(conv_df_fun))
  conv_df_fun$conv<-rep(NA, nrow(conv_df_fun))
  conv_df_fun$pheno<-rep(NA, nrow(conv_df_fun))
  for(i in 1:nrow(conv_df_fun)){
    pert_1<-conv_df_fun$treatment1[i]
    pert_2<-conv_df_fun$treatment2[i]
    if(!(pert_1 %in% all_gene_focus) || !(pert_2 %in% all_gene_focus)) {
      next
    }
    rowA<-which(conv_df_fun$treatment1 == pert_1 & conv_df_fun$treatment2 == pert_2)
    rowB<-which(conv_df_fun$treatment1 == pert_2 & conv_df_fun$treatment2 == pert_1)
    pA<-conv_df_fun$p_value[rowA]
    pB<-conv_df_fun$p_value[rowB]
    if(p_method=='max'){
      conv_df_fun$conv_p[i]<-max(pA, pB)
    }else if(p_method=='min'){
      conv_df_fun$conv_p[i]<-min(pA, pB)
    }
    
    if(scLinker_use[scLinker_use$Geneset==pert_1, 'Enrichment_score_pvalue']<pheno_thres
       & scLinker_use[scLinker_use$Geneset==pert_2, 'Enrichment_score_pvalue']<pheno_thres){
      conv_df_fun$pheno[i]<-'Significant'
    }else{
      conv_df_fun$pheno[i]<-'Not_Significant'
    }
  }
  conv_df_fun <- conv_df_fun[!is.na(conv_df_fun$conv_p), ]
  conv_df_fun$conv<-ifelse(conv_df_fun$conv_p<conv_thres, 'Significant', 'Not_Significant')
  conv_df_fun$combined <- ifelse(conv_df_fun$conv == "Significant" & conv_df_fun$pheno == "Significant",
                                 "Significant",
                                 "Not_Significant")
  return(list(conv_df_return=conv_df_fun, all_gene_focus_return=all_gene_focus))
}

pheno_thres<-0.05

#####################################################################################################################
pheno = "Eczema" 
knock_way<-c('(pos) Knock-out')
combo_clean <- gsub("[()]", "", knock_way)       # 去掉括号
combo_clean <- gsub(" ", "-", combo_clean)   # 空格换成 -
combo_str <- paste(combo_clean, collapse = "_")  # 用 _ 连接
combo_str
result<-scLinker_pheno_conv(scLinker_df_fun=scLinker_df, conv_df_fun=direct_df, phenotype=pheno, knock_way=knock_way, p_method='max', pheno_thres=pheno_thres)

df_use<-result$conv_df_return
write.csv(df_use, paste0(plot_path, '/', pheno, "_df_direct.csv"))
all_gene_focus<-sort(result$all_gene_focus_return)

df_use$treatment1 <- factor(df_use$treatment1, levels = all_gene_focus)
df_use$treatment2 <- factor(df_use$treatment2, levels = all_gene_focus)
df_no_na<-df_use[!is.na(df_use$conv), ]

df_upper <- df_no_na[match(df_no_na$treatment1, all_gene_focus) < match(df_no_na$treatment2, all_gene_focus), ]
df_upper$conv <- factor(df_upper$conv)
df_upper$combined <- factor(df_upper$combined)
df_lower <- df_no_na[match(df_no_na$treatment1, all_gene_focus) > match(df_no_na$treatment2, all_gene_focus), ]
df_lower$conv <- factor(df_lower$conv)
df_lower$combined <- factor(df_lower$combined)

df_upper$color_group <- ifelse(df_upper$combined == "Significant", "Both",
                               ifelse(df_upper$pheno == "Significant", "Pheno_only", "None"))


df_lower$color_group <- ifelse(df_lower$combined == "Significant", "Both",
                               ifelse(df_lower$conv == "Significant", "Direct_only", "None"))

p1 <- ggplot() +
  geom_point(data = df_upper,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  geom_point(data = df_lower,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  scale_color_manual(
    name = "Convergence",
    values = c(
      "Both" = "#f5b70a",
      "Pheno_only" = "#013e75",
      "Direct_only" = "#663B8C",
      "None" = alpha("#EEEEEE", 0.4)
    ),
    breaks = c("Both", "Pheno_only", "Direct_only", "None"),
    drop = FALSE
  )+
  scale_x_discrete(limits = all_gene_focus) +
  scale_y_discrete(limits = all_gene_focus) +
  labs(title = pheno,
       x = "Gene 1", y = "Gene 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank()   # 移除外框
  )
print(p1)
#####################################################################################################################
pheno = "Primary biliary cirrhosis" 
knock_way<-c('(pos) Knock-out')
combo_clean <- gsub("[()]", "", knock_way)       # 去掉括号
combo_clean <- gsub(" ", "-", combo_clean)   # 空格换成 -
combo_str <- paste(combo_clean, collapse = "_")  # 用 _ 连接
combo_str
result<-scLinker_pheno_conv(scLinker_df_fun=scLinker_df, conv_df_fun=direct_df, phenotype=pheno, knock_way=knock_way, p_method='max', pheno_thres=pheno_thres)

df_use<-result$conv_df_return
write.csv(df_use, paste0(plot_path, '/', pheno, "_df_direct.csv"))
all_gene_focus<-sort(result$all_gene_focus_return)

df_use$treatment1 <- factor(df_use$treatment1, levels = all_gene_focus)
df_use$treatment2 <- factor(df_use$treatment2, levels = all_gene_focus)
df_no_na<-df_use[!is.na(df_use$conv), ]

df_upper <- df_no_na[match(df_no_na$treatment1, all_gene_focus) < match(df_no_na$treatment2, all_gene_focus), ]
df_upper$conv <- factor(df_upper$conv)
df_upper$combined <- factor(df_upper$combined)
df_lower <- df_no_na[match(df_no_na$treatment1, all_gene_focus) > match(df_no_na$treatment2, all_gene_focus), ]
df_lower$conv <- factor(df_lower$conv)
df_lower$combined <- factor(df_lower$combined)

df_upper$color_group <- ifelse(df_upper$combined == "Significant", "Both",
                               ifelse(df_upper$pheno == "Significant", "Pheno_only", "None"))


df_lower$color_group <- ifelse(df_lower$combined == "Significant", "Both",
                               ifelse(df_lower$conv == "Significant", "Direct_only", "None"))

p2 <- ggplot() +
  geom_point(data = df_upper,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  geom_point(data = df_lower,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  scale_color_manual(
    name = "Convergence",
    values = c(
      "Both" = "#f5b70a",
      "Pheno_only" = "#013e75",
      "Direct_only" = "#663B8C",
      "None" = alpha("#EEEEEE", 0.4)
    ),
    breaks = c("Both", "Pheno_only", "Direct_only", "None"),
    drop = FALSE
  )+
  scale_x_discrete(limits = all_gene_focus) +
  scale_y_discrete(limits = all_gene_focus) +
  labs(title = pheno,
       x = "Gene 1", y = "Gene 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank()   # 移除外框
  )
print(p2)
#####################################################################################################################
pheno = "Eosinophil percentage" 
knock_way<-c('(pos) Knock-out')
combo_clean <- gsub("[()]", "", knock_way)       # 去掉括号
combo_clean <- gsub(" ", "-", combo_clean)   # 空格换成 -
combo_str <- paste(combo_clean, collapse = "_")  # 用 _ 连接
combo_str
result<-scLinker_pheno_conv(scLinker_df_fun=scLinker_df, conv_df_fun=direct_df, phenotype=pheno, knock_way=knock_way, p_method='max', pheno_thres=pheno_thres)

df_use<-result$conv_df_return
write.csv(df_use, paste0(plot_path, '/', pheno, "_df_direct.csv"))
all_gene_focus<-sort(result$all_gene_focus_return)

df_use$treatment1 <- factor(df_use$treatment1, levels = all_gene_focus)
df_use$treatment2 <- factor(df_use$treatment2, levels = all_gene_focus)
df_no_na<-df_use[!is.na(df_use$conv), ]

df_upper <- df_no_na[match(df_no_na$treatment1, all_gene_focus) < match(df_no_na$treatment2, all_gene_focus), ]
df_upper$conv <- factor(df_upper$conv)
df_upper$combined <- factor(df_upper$combined)
df_lower <- df_no_na[match(df_no_na$treatment1, all_gene_focus) > match(df_no_na$treatment2, all_gene_focus), ]
df_lower$conv <- factor(df_lower$conv)
df_lower$combined <- factor(df_lower$combined)

df_upper$color_group <- ifelse(df_upper$combined == "Significant", "Both",
                               ifelse(df_upper$pheno == "Significant", "Pheno_only", "None"))


df_lower$color_group <- ifelse(df_lower$combined == "Significant", "Both",
                               ifelse(df_lower$conv == "Significant", "Direct_only", "None"))

p3 <- ggplot() +
  geom_point(data = df_upper,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  geom_point(data = df_lower,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  scale_color_manual(
    name = "Convergence",
    values = c(
      "Both" = "#f5b70a",
      "Pheno_only" = "#013e75",
      "Direct_only" = "#663B8C",
      "None" = alpha("#EEEEEE", 0.4)
    ),
    breaks = c("Both", "Pheno_only", "Direct_only", "None"),
    drop = FALSE
  )+
  scale_x_discrete(limits = all_gene_focus) +
  scale_y_discrete(limits = all_gene_focus) +
  labs(title = pheno,
       x = "Gene 1", y = "Gene 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank()   # 移除外框
  )
print(p3)
#####################################################################################################################
pheno = "Inflammatory bowel disease" 
knock_way<-c('(pos) Knock-out')
combo_clean <- gsub("[()]", "", knock_way)       # 去掉括号
combo_clean <- gsub(" ", "-", combo_clean)   # 空格换成 -
combo_str <- paste(combo_clean, collapse = "_")  # 用 _ 连接
combo_str
result<-scLinker_pheno_conv(scLinker_df_fun=scLinker_df, conv_df_fun=direct_df, phenotype=pheno, knock_way=knock_way, p_method='max', pheno_thres=pheno_thres)

df_use<-result$conv_df_return
write.csv(df_use, paste0(plot_path, '/', pheno, "_df_direct.csv"))
all_gene_focus<-sort(result$all_gene_focus_return)

df_use$treatment1 <- factor(df_use$treatment1, levels = all_gene_focus)
df_use$treatment2 <- factor(df_use$treatment2, levels = all_gene_focus)
df_no_na<-df_use[!is.na(df_use$conv), ]

df_upper <- df_no_na[match(df_no_na$treatment1, all_gene_focus) < match(df_no_na$treatment2, all_gene_focus), ]
df_upper$conv <- factor(df_upper$conv)
df_upper$combined <- factor(df_upper$combined)
df_lower <- df_no_na[match(df_no_na$treatment1, all_gene_focus) > match(df_no_na$treatment2, all_gene_focus), ]
df_lower$conv <- factor(df_lower$conv)
df_lower$combined <- factor(df_lower$combined)

df_upper$color_group <- ifelse(df_upper$combined == "Significant", "Both",
                               ifelse(df_upper$pheno == "Significant", "Pheno_only", "None"))


df_lower$color_group <- ifelse(df_lower$combined == "Significant", "Both",
                               ifelse(df_lower$conv == "Significant", "Direct_only", "None"))

p4 <- ggplot() +
  geom_point(data = df_upper,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  geom_point(data = df_lower,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  scale_color_manual(
    name = "Convergence",
    values = c(
      "Both" = "#f5b70a",
      "Pheno_only" = "#013e75",
      "Direct_only" = "#663B8C",
      "None" = alpha("#EEEEEE", 0.4)
    ),
    breaks = c("Both", "Pheno_only", "Direct_only", "None"),
    drop = FALSE
  )+
  scale_x_discrete(limits = all_gene_focus) +
  scale_y_discrete(limits = all_gene_focus) +
  labs(title = pheno,
       x = "Gene 1", y = "Gene 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank()   # 移除外框
  )
print(p4)
#####################################################################################################################
pheno = "Rheumatoid arthritis" 
knock_way<-c('(pos) Knock-out')
combo_clean <- gsub("[()]", "", knock_way)       # 去掉括号
combo_clean <- gsub(" ", "-", combo_clean)   # 空格换成 -
combo_str <- paste(combo_clean, collapse = "_")  # 用 _ 连接
combo_str
result<-scLinker_pheno_conv(scLinker_df_fun=scLinker_df, conv_df_fun=direct_df, phenotype=pheno, knock_way=knock_way, p_method='max', pheno_thres=pheno_thres)

df_use<-result$conv_df_return
write.csv(df_use, paste0(plot_path, '/', pheno, "_df_direct.csv"))
all_gene_focus<-sort(result$all_gene_focus_return)

df_use$treatment1 <- factor(df_use$treatment1, levels = all_gene_focus)
df_use$treatment2 <- factor(df_use$treatment2, levels = all_gene_focus)
df_no_na<-df_use[!is.na(df_use$conv), ]

df_upper <- df_no_na[match(df_no_na$treatment1, all_gene_focus) < match(df_no_na$treatment2, all_gene_focus), ]
df_upper$conv <- factor(df_upper$conv)
df_upper$combined <- factor(df_upper$combined)
df_lower <- df_no_na[match(df_no_na$treatment1, all_gene_focus) > match(df_no_na$treatment2, all_gene_focus), ]
df_lower$conv <- factor(df_lower$conv)
df_lower$combined <- factor(df_lower$combined)

df_upper$color_group <- ifelse(df_upper$combined == "Significant", "Both",
                               ifelse(df_upper$pheno == "Significant", "Pheno_only", "None"))


df_lower$color_group <- ifelse(df_lower$combined == "Significant", "Both",
                               ifelse(df_lower$conv == "Significant", "Direct_only", "None"))

p5 <- ggplot() +
  geom_point(data = df_upper,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  geom_point(data = df_lower,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  scale_color_manual(
    name = "Convergence",
    values = c(
      "Both" = "#f5b70a",
      "Pheno_only" = "#013e75",
      "Direct_only" = "#663B8C",
      "None" = alpha("#EEEEEE", 0.4)
    ),
    breaks = c("Both", "Pheno_only", "Direct_only", "None"),
    drop = FALSE
  )+
  scale_x_discrete(limits = all_gene_focus) +
  scale_y_discrete(limits = all_gene_focus) +
  labs(title = pheno,
       x = "Gene 1", y = "Gene 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank()   # 移除外框
  )
print(p5)

##############################################################################################
phenos<-c("Eczema", "Primary biliary cirrhosis", "Eosinophil percentage", "Inflammatory bowel disease", "Rheumatoid arthritis")
knock_way<-c('(pos) Knock-out')
combo_clean <- gsub("[()]", "", knock_way)       # 去掉括号
combo_clean <- gsub(" ", "-", combo_clean)   # 空格换成 -
combo_str <- paste(combo_clean, collapse = "_")  # 用 _ 连接
combo_str
results<-vector("list", length(phenos))
for(i in 1:length(phenos)){
  pheno<-phenos[i]
  result<-scLinker_pheno_conv(scLinker_df_fun=scLinker_df, conv_df_fun=direct_df, phenotype=pheno, knock_way=knock_way, p_method='max', pheno_thres=pheno_thres)
  results[[i]]<-result
}

at_least_one_df<-results[[1]]$conv_df_return
for(i in 1:nrow(at_least_one_df)){
  t1 <- at_least_one_df$treatment1[i]
  t2 <- at_least_one_df$treatment2[i]
  
  subset_list <- lapply(results, function(res) {
    df <- res$conv_df_return
    df[df$treatment1 == t1 & df$treatment2 == t2, c("conv", "pheno", "combined")]
  })
  
  combined_rows <- do.call(rbind, subset_list)
  
  for (col in c("conv", "pheno", "combined")) {
    if (any(combined_rows[[col]] == "Significant")) {
      at_least_one_df[i, col] <- "Significant"
    } else {
      at_least_one_df[i, col] <- "Not_Significant"
    }
  }
}

write.csv(at_least_one_df, paste0(plot_path, "/all_combined_at_least_one_df_direct.csv"))
all_gene_focus<-sort(results[[1]]$all_gene_focus_return)

at_least_one_df$treatment1 <- factor(at_least_one_df$treatment1, levels = all_gene_focus)
at_least_one_df$treatment2 <- factor(at_least_one_df$treatment2, levels = all_gene_focus)
at_least_one_df_no_na<-at_least_one_df[!is.na(at_least_one_df$conv), ]

df_upper <- at_least_one_df_no_na[match(at_least_one_df_no_na$treatment1, all_gene_focus) < match(at_least_one_df_no_na$treatment2, all_gene_focus), ]
df_upper$conv <- factor(df_upper$conv)
df_upper$combined <- factor(df_upper$combined)
df_lower <- at_least_one_df_no_na[match(at_least_one_df_no_na$treatment1, all_gene_focus) > match(at_least_one_df_no_na$treatment2, all_gene_focus), ]
df_lower$conv <- factor(df_lower$conv)
df_lower$combined <- factor(df_lower$combined)

df_upper$color_group <- ifelse(df_upper$combined == "Significant", "Both",
                               ifelse(df_upper$pheno == "Significant", "Pheno_only", "None"))


df_lower$color_group <- ifelse(df_lower$combined == "Significant", "Both",
                               ifelse(df_lower$conv == "Significant", "Direct_only", "None"))

p <- ggplot() +
  geom_point(data = df_upper,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  geom_point(data = df_lower,
             aes(x = treatment1, y = treatment2, color = color_group),
             size = 3) +
  scale_color_manual(
    name = "Convergence",
    values = c(
      "Both" = "#f5b70a",
      "Pheno_only" = "#013e75",
      "Direct_only" = "#663B8C",
      "None" = alpha("#EEEEEE", 0.4)
    ),
    breaks = c("Both", "Pheno_only", "Direct_only", "None"),
    drop = FALSE
  )+
  scale_x_discrete(limits = all_gene_focus) +
  scale_y_discrete(limits = all_gene_focus) +
  labs(title = 'Combined',
       x = "Gene 1", y = "Gene 2") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.border = element_blank()   # 移除外框
  )
print(p)

# ggsave(paste0(plot_path, "combined_at_least_one.pdf"), p, width = 10, height = 8) # 5.5, 4

combined_plot <- p1 + p2 + p3 + p4 + p5 + p+
  plot_layout(ncol=2, guides = "collect") & 
  theme(legend.position = "bottom")

combined_plot <- combined_plot +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(
      face = "bold",   # 粗体
      size = 25        # 字号（可调整）
    )
  )

print(combined_plot)

ggsave(paste0(plot_path, "/all_combined_and_combined.pdf"), combined_plot, width = 20, height = 25) # 5.5, 4






