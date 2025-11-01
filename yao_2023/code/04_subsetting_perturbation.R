 
residual_matrix <- readRDS(file = file.path(work_directory, 
                                            'yao_2023/data/intermediate_data/residual_matrix.rds'))

perturbation_interest <- c('non-targeting', 'IFNAR1','TYK2','STAT1','STAT2','YEATS4','MAP2K1','MAP2K2','AHR','MAPK14',
                            'PRDM1','ADO','JUNB','RAB5C','MAPK1','EHMT2','ATXN7L3','DNTTIP1','CYLD',
                            'TXNDC17','ARPC4','TNFAIP3','PICALM','CEBPG','SEPHS2','STK11','FLI1',
                            'TRIB1','MEF2C','KIDINS220','XPR1','SBDS','PPP2R1A','PGM3','MIDN',
                            'CSNK1A1','RNF31','TAB1','TLR1','IKBKB','TIRAP','ALG2','RELA','IKBKG',
                            'MAP3K7','TLR2','MYD88','HSP90B1','IRAK4','IRAK1','TRAF6')
length(perturbation_interest)

residual_subset <- residual_matrix[residual_matrix$Guides_collapsed_by_gene %in% perturbation_interest, ]
length(table(residual_subset$Guides_collapsed_by_gene))

perturbation_interest[!perturbation_interest %in% names(table(residual_subset$Guides_collapsed_by_gene))]
saveRDS(residual_subset, file.path(work_directory, 
                                   'yao_2023/data/intermediate_data/residual_matrix_all_in_paper.rds'))

