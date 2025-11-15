library(parallel)
library(data.table)
library(HMC, lib.loc = "~/Rlibs")

data_path<-'/home/eshang/diffusion_and_protein/diff_conv/eff_emb/data/sampled_data_csv/'
result_path<-'/home/eshang/diffusion_and_protein/diff_conv/eff_emb/data/result/alter_increase_1_15_and_1_to_5_fac1_095_fac2_09/'
sample_use<-c(10, 30, 50, 80, 100, 300, 500)
fac1<-0.95
fac2<-0.9
count<-1

num_cores<-20

run_single_seed<-function(i, n, data_path, count, fac1, fac2, result_path){
  tryCatch({
    cat(paste0("For sample size ", n, " and seed ", i, "\n"))
    set.seed(i)
    
    pert1<-as.matrix(read.csv(paste0(data_path, 'non-targeting_v3_eta1.5_5000gen_acc_wo_ZC3H13_seed_', i, '.csv')))[count:(count+n-1), ]
    pert2<-as.matrix(read.csv(paste0(data_path, 'fake_pert_v3_eta1.5_5000gen_acc_wo_ZC3H13_seed_', i+2000, '.csv')))[count:(count+n-1), ]
    ctrl<-as.matrix(read.csv(paste0(data_path, 'non-targeting_v3_eta1.5_5000gen_acc_wo_ZC3H13_seed_', i+4000, '.csv')))[count:(count+n-1), ]
    
    
    pert1[, 1:15]<-pert1[, 1:15] * fac1
    pert2[, 1:5]<-pert2[, 1:5] * fac2
    
    cat("for one\n")
    test_result_one<-convergence_testing(
      control = ctrl,
      treatment1 = pert1,
      treatment2 = pert2,
      pca_method = "dense_pca",
      classifier_method = "lasso",
      lambda_type = "lambda.min",
      n_folds = 5,
      verbose = FALSE
    )
    formatted_fac1 <- sprintf("%03d", as.integer(fac1 * 100))
    formatted_fac2 <- sprintf("%02d", as.integer(fac2 * 10))
    save(test_result_one, file = paste0(result_path, 'non_tar_1_15_increase_and_1_5_increase_seed_', i, '_sample_', n, 'fac1_', formatted_fac1, 'fac2_', formatted_fac2, '_one.Rdata'))
    
    cat("for two\n")
    test_result_two<-convergence_testing(
      control = ctrl,
      treatment1 = pert2,
      treatment2 = pert1,
      pca_method = "dense_pca",
      classifier_method = "lasso",
      lambda_type = "lambda.min",
      n_folds = 5,
      verbose = FALSE
    )
    save(test_result_two, file = paste0(result_path, 'non_tar_1_15_increase_and_1_5_increase_seed_', i, '_sample_', n, 'fac1_', formatted_fac1, 'fac2_', formatted_fac2, '_two.Rdata'))
    
    return(NULL)
  }, error = function(e) {
    cat(paste0("Error at seed ", i, ": ", e$message, "\n"))
    return(NULL)
  })
}

for(n in sample_use){
  mclapply(
    1:2000,
    run_single_seed,
    n = n,
    data_path = data_path,
    count = count,
    fac1 = fac1,
    fac2 = fac2,
    result_path = result_path,
    mc.cores = num_cores
  )
  count<-count+n
}








