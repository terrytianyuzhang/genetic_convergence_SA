result_path<-'/home/eshang/diffusion_and_protein/diff_conv/eff_emb/data/result/increase_1_30_and_increase_31_60/'
sample_use<-c(100, 200, 500, 800, 1000)
seed_num<-2000
fac<-5
save_paper_path<-'/home/eshang/diffusion_and_protein/diff_conv/eff_emb/data/for_paper/'
for(n in sample_use){
  print(paste0('The sample size is ', n))
  p_vec_one<-numeric(seed_num)
  p_vec_two<-numeric(seed_num)
  p_vec_max<-numeric(seed_num)
  p_vec_min<-numeric(seed_num)
  for(i in 1:seed_num){
    load(paste0(result_path, 'non_tar_1_30_increase_and_31_60_increase_seed_', i, '_sample_', n, 'fac_', fac, '_one.Rdata'))
    p_vec_one[i]<-test_result_one$p_value
    load(paste0(result_path, 'non_tar_1_30_increase_and_31_60_increase_seed_', i, '_sample_', n, 'fac_', fac, '_two.Rdata'))
    p_vec_two[i]<-test_result_two$p_value
    
    p_vec_max[i]<-max(p_vec_one[i], p_vec_two[i])
    p_vec_min[i]<-min(p_vec_one[i], p_vec_two[i])
  }
  print(paste0('Rejection rate(max_p) is ', sum(p_vec_max<0.05)/length(p_vec_max), ' when sample size is ', n))
  print(paste0('Rejection rate(min_p) is ', sum(p_vec_min<0.05)/length(p_vec_min), ' when sample size is ', n))
  save(p_vec_one, file=paste0(save_paper_path, 'p_vec_one_null_sample_', n, '.Rdata'))
  save(p_vec_two, file=paste0(save_paper_path, 'p_vec_two_null_sample_', n, '.Rdata'))
}
  

