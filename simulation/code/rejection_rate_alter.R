result_path<-'/home/eshang/diffusion_and_protein/diff_conv/eff_emb/data/result/alter_increase_1_15_and_1_to_5_fac1_095_fac2_09/'
sample_use<-c(10, 30, 50, 80, 100, 300, 500)
seed_num<-2000
fac1<-0.95
fac2<-0.9
formatted_fac1 <- sprintf("%03d", as.integer(fac1 * 100))
formatted_fac2 <- sprintf("%02d", as.integer(fac2 * 10))
save_paper_path<-'/home/eshang/diffusion_and_protein/diff_conv/eff_emb/data/for_paper/'
for(n in sample_use){
  print(paste0('The sample size is ', n))
  p_vec_one<-numeric(seed_num)
  p_vec_two<-numeric(seed_num)
  p_vec_max<-numeric(seed_num)
  p_vec_min<-numeric(seed_num)
  for(i in 1:seed_num){
    load(paste0(result_path, 'non_tar_1_15_increase_and_1_5_increase_seed_', i, '_sample_', n, 'fac1_', formatted_fac1, 'fac2_', formatted_fac2, '_one.Rdata'))
    p_vec_one[i]<-test_result_one$p_value
    load(paste0(result_path, 'non_tar_1_15_increase_and_1_5_increase_seed_', i, '_sample_', n, 'fac1_', formatted_fac1, 'fac2_', formatted_fac2, '_two.Rdata'))
    p_vec_two[i]<-test_result_two$p_value
    
    p_vec_max[i]<-max(p_vec_one[i], p_vec_two[i])
    p_vec_min[i]<-min(p_vec_one[i], p_vec_two[i])
  }
  print(paste0('Power(max_p) is ', sum(p_vec_max<0.05)/length(p_vec_max), ' when sample size is ', n))
  print(paste0('Power(min_p) is ', sum(p_vec_min<0.05)/length(p_vec_min), ' when sample size is ', n))
  save(p_vec_one, file=paste0(save_paper_path, 'p_vec_one_alter_sample_', n, '.Rdata'))
  save(p_vec_two, file=paste0(save_paper_path, 'p_vec_two_alter_sample_', n, '.Rdata'))
}





