library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/')

saige_NK <- fread('Saige/step3/outputs/NK_NI_best_eQTLs_no_perm_withqvalues.txt')
saige_CD4T <- fread('Saige/step3/outputs/CD4T_NI_best_eQTLs_no_perm_withqvalues.txt')

matrixeqtl_NK <- fread('HALEYs/matrixEQTL_results/NI_NK_adj_3PCs_best_cisQTL_withqvalue_sumstats.txt')
matrixeqtl_CD4T <- fread('HALEYs/matrixEQTL_results/NI_T-CD4_adj_4PCs_best_cisQTL_withqvalue_sumstats.txt')

NK_joint <- inner_join(matrixeqtl_NK, saige_NK, by=c('gene'))
NK_joint_summary <- data.frame('matrix'=sum(NK_joint$qvals.x<0.05),
                               'saigeqtl'=sum(NK_joint$qvals.y<0.05))

CD4T_joint <- inner_join(matrixeqtl_CD4T, saige_CD4T, by=c('gene'))
CD4T_joint_summary <- data.frame('matrix'=sum(CD4T_joint$qvals.x<0.05),
                               'saigeqtl'=sum(CD4T_joint$qvals.y<0.05))
