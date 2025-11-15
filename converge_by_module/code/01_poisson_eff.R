save_path<-file.path(work_directory, 'converge_by_module', 'report', 'heatmap')
dir.create(save_path, recursive = TRUE)

# Helper: expand module-level specs into a per-feature vector
.expand_module_vector <- function(x, n_modules, module_size, n_features, baseline, seed=NULL) {
  if (!is.null(seed)) {
    # 用局部随机种子，保证 reproducibility 且不影响全局随机流
    return(withr::with_seed(seed, .expand_module_vector(x, n_modules, module_size, n_features, baseline, seed = NULL)))
  }
  if (length(x) == 1L) {
    rep(x, n_features)
  } else if (length(x) == n_modules) {
    # rep(x, each = module_size)
    out <- numeric(n_modules * module_size)
    for (i in seq_len(n_modules)) {
      vals <- rep(x[i], module_size)
      n_replace <- ceiling(0.8 * module_size)
      idx_replace <- sample(seq_len(module_size), n_replace)
      vals[idx_replace] <- baseline
      start <- (i - 1) * module_size + 1
      end <- i * module_size
      out[start:end] <- vals
    }
    out
  } else if (length(x) == n_features) {
    x
  } else {
    stop("Length of vector must be 1, n_modules, or n_features.")
  }
}

# Helper: build block-diagonal covariance from per-module SDs or a single SD
.build_block_cov <- function(sd_within, n_modules, module_size, n_features) {
  sds <- rep(sd_within, each = module_size) # .expand_module_vector(sd_within, n_modules, module_size, n_features)
  diag(sds^2)
}

# Helper: feature names as feature_(moduleIndex)_(orderWithinModule)
.feature_names <- function(n_modules, module_size, prefix = "feature") {
  mods <- rep(seq_len(n_modules), each = module_size)
  ords <- rep(seq_len(module_size), times = n_modules)
  paste0(prefix, "_", mods, "_", ords)
}

# -------- Modular simulator --------
simulate_modular_data <- function(
    n_samples   = 500,
    n_modules   = 7,            # number of modules
    module_size = 30,           # dimensions per module
    ratio       = 3,            # control:treatment sample ratio
    # Per-module mean shifts (can be length 1, n_modules, or n_features)
    mu_control  = 0,
    mu_tr_list,
    # Per-module within-feature SDs (same length rules as above)
    sd_within   = 1,
    # Alternatively supply full covariance matrices (overrides sd_within if given)
    Sigma_control = NULL,
    Sigma_tr1     = NULL,
    Sigma_tr2     = NULL,
    feature_prefix = "feature"
) {
  n_features <- n_modules * module_size
  
  # Expand means
  mu_expand_list<-list()
  
  mu_con <- rep(mu_control, each = module_size) #.expand_module_vector(mu_control, n_modules, module_size, n_features)
  for(i in 1:length(mu_tr_list)){
    mu_expand_list[[i]]<-.expand_module_vector(mu_tr_list[[i]], n_modules, module_size, n_features, baseline=mu_control[1], seed=i)
  }
  names(mu_expand_list)<-names(mu_tr_list)
  
  # Build covariance
  sd_expand_list<-list()
  if (is.null(Sigma_control)) Sigma_control <- .build_block_cov(sd_within, n_modules, module_size, n_features)
  
  for(i in 1:length(mu_tr_list)){
    sd_expand_list[[i]]<-.build_block_cov(sd_within, n_modules, module_size, n_features)
  }
  names(sd_expand_list)<-names(mu_tr_list)
  
  # Draw
  tr_sample<-list()
  control    <- round(MASS::mvrnorm(n_samples * ratio, mu = mu_con, Sigma = Sigma_control))
  for(i in 1:length(mu_tr_list)){
    tr_sample[[i]]<-round(MASS::mvrnorm(n_samples, mu = mu_expand_list[[i]],  Sigma = sd_expand_list[[i]]))
  }
  names(tr_sample)<-names(mu_tr_list)
  # print(Sigma_control)
  # Names
  rn_c  <- paste0("Sample_", seq_len(n_samples * ratio))
  cn    <- .feature_names(n_modules, module_size, prefix = feature_prefix)
  for(i in 1:length(tr_sample)){
    rownames(tr_sample[[i]])<-paste0("Sample_", seq_len(n_samples))
    colnames(tr_sample[[i]])<-cn
  }
  rownames(control)    <- rn_c
  colnames(control)    <- cn
  
  list(
    control = control,
    treatment = tr_sample,
    mean = mu_expand_list,
    meta = list(
      n_samples = n_samples,
      n_modules = n_modules,
      module_size = module_size,
      n_features = n_features
    )
  )
}

set.seed(4)
n_samples<-200
n_modules<-30
module_sz<-50


mu_control_mods <- rep(10, n_modules)   # length n_modules

mu_tr_mods_list<-list('A'=10*c(1.4, 1.5, 1.45, 1.4, 1.3, rep(1, n_modules-5)),
                      'B'=10*c(0.95, 0.5, 0.45, 0.55, 0.95, rep(1, n_modules-5)),
                      'C'=10*c(rep(1, 5), 1.4, 1.5, 1.45, 1.4, 1.3, rep(1, 20)),
                      'D'=10*c(0.5, 1.2, 1.5, 0.6, 1.45, rep(1, n_modules-5)),
                      'E'=10*c(1, 1.45, 1.4, 1, 1, 1, 0.45, 0.5, rep(1, n_modules-8)),
                      'F'=10*c(rep(1, 11), 1.4, 1.5, 1.45, 1.4, 1.3, rep(1, 14))
)

sim <- simulate_modular_data(
  n_samples   = n_samples,
  n_modules   = n_modules,
  module_size = module_sz,
  ratio       = 2,
  mu_control  = mu_control_mods,               # scalar => 0 everywhere
  mu_tr_mods_list,
  sd_within = 1*rep(1,n_modules)
)
ratio<-2

################################ eff calculated ################################################################
tr_all<-do.call(rbind, sim$treatment)
Y<-rbind(tr_all, sim$control)
Y_df<-as.data.frame(Y)
Y_df$pert<-c(rep("A", n_samples), rep("B", n_samples), rep("C", n_samples), rep("D", n_samples),
             rep("E", n_samples), rep("F", n_samples), rep("NT", n_samples*ratio))
Y_df$pert<-factor(Y_df$pert)
Y_df$pert<-relevel(Y_df$pert, ref='NT')

mat_res<-matrix(0, 6, ncol(Y))
rownames(mat_res)<-c('A', 'B', 'C', 'D', 'E', 'F')
colnames(mat_res)<-colnames(tr_all)

for(i in 1:(ncol(Y_df)-1)){
  model<-glm(Y_df[, i]~pert, family=poisson, data=Y_df)
  for(j in 1:length(rownames(mat_res))){
    letter<-rownames(mat_res)[j]
    mat_res[letter, i]<-exp(coef(model)[paste0('pert', letter)]+coef(model)['(Intercept)'])/exp(coef(model)['(Intercept)'])
  }
}
mat_res<-log2(mat_res)


program_mat<-matrix(0, nrow(mat_res), ncol(mat_res)/module_sz)
rownames(program_mat)<-c("A", "B", "C", 'D', 'E', 'F')
for(i in 1:nrow(program_mat)){
  count<-0
  for(j in 1:ncol(program_mat)){
    # gene_name<-old_gene_name[(count+1):(count+gene_num)]
    sub_eff_mat<-mat_res[i, (count+1):(count+module_sz)]
    program_mat[i, j]<-mean(sub_eff_mat)
    count<-count+module_sz
  }
}



mat_long <- melt(program_mat)
p <- ggplot(mat_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#013e75", mid = "#fAfAfA", high = "#f5b70a", midpoint = 0) +
  theme_minimal() +
  labs(x = "Program", y = "Perturbation", fill = "log2FC") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12),   # x 轴标题大小
    axis.title.y = element_text(size = 12),   # y 轴标题大小
    axis.text = element_text(size = 12),      # 坐标刻度文字大小
    legend.title = element_text(size = 10),   # 图例标题大小
    legend.text = element_text(size = 12),    # 图例刻度文字大小
    legend.position = "top"                   # 图例放在上方
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",   # 标题在色条上方
      title.hjust = 0.5,        # 标题居中
      barwidth = 20,            # 色条变长 (可调大)
      barheight = 0.5             # 色条高度 (可适当调小)
    )
  )

print(p)


ggsave(paste0(save_path, "/logFoldChange_artificial_program_pert_poisson_200_400.pdf"), plot = p)






