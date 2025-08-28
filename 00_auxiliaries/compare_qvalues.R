library(tidyverse)
library(data.table)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/')

### COMPILE SAIGE-QTL RESULTS
# NK
for (perm in c('no_perm', 'perm1', 'perm2', 'perm3', 'perm4', 'perm5')){
  print(perm)
  if (perm=='no_perm'){
    saige_files <- list.files('Saige/step3/outputs/NK_NI/perms/no_perm/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) 
      
      # combine into a sigle dataframe
      if (exists('best_NK_s')){
        best_NK_s <- rbind(best_NK_s, tmp)
      } else {best_NK_s <- tmp}
    }
  } else {
    saige_files <- list.files('Saige/step3/outputs/NK_NI/perms/'%&%perm%&%'/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) %>% mutate(permutation=perm)
      
      # combine into a sigle dataframe
      if (exists('perms_NK_s')){
        perms_NK_s <- rbind(perms_NK_s, tmp)
      } else {perms_NK_s <- tmp}
    }
  }
}
#fwrite(best_NK_s, 'Saige/step3/outputs/NK_NI_best_eQTLs_no_perm.txt')
#fwrite(perms_NK_s, 'Saige/step3/outputs/NK_NI_best_eQTLs_perm.txt')
best_NK_s <- fread('Saige/step3/outputs/NK_NI_best_eQTLs_no_perm.txt')
perms_NK_s <- fread('Saige/step3/outputs/NK_NI_best_eQTLs_perm.txt')
best_NK_s <- best_NK_s %>% arrange(gene)
perms_NK_s <- perms_NK_s %>% arrange(permutation, gene)

# get top_pval for each permutation in wide format
perms_NK_s <- perms_NK_s %>%
  select(gene, permutation, top_pval) %>%
  pivot_wider(
    names_from = permutation,
    values_from = top_pval,
    names_prefix = 'perm_')

# CD4 T
for (perm in c('no_perm', 'perm1', 'perm2', 'perm3', 'perm4', 'perm5')){
  print(perm)
  if (perm=='no_perm'){
    saige_files <- list.files('Saige/step3/outputs/CD4_T_NI/perms/no_perm/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) 
      
      # combine into a sigle dataframe
      if (exists('best_CD4T_s')){
        best_CD4T_s <- rbind(best_CD4T_s, tmp)
      } else {best_CD4T_s <- tmp}
    }
  } else {
    saige_files <- list.files('Saige/step3/outputs/CD4_T_NI/perms/'%&%perm%&%'/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) %>% mutate(permutation=perm)
      
      # combine into a sigle dataframe
      if (exists('perms_CD4T_s')){
        perms_CD4T_s <- rbind(perms_CD4T_s, tmp)
      } else {perms_CD4T_s <- tmp}
    }
  }
}
#fwrite(best_CD4T_s, 'Saige/step3/outputs/CD4T_NI_best_eQTLs_no_perm.txt')
#fwrite(perms_CD4T_s, 'Saige/step3/outputs/CD4T_NI_best_eQTLs_perm.txt')
best_CD4T_s <- fread('Saige/step3/outputs/CD4T_NI_best_eQTLs_no_perm.txt')
perms_CD4T_s <- fread('Saige/step3/outputs/CD4T_NI_best_eQTLs_perm.txt')
best_CD4T_s <- best_CD4T_s %>% arrange(gene)
perms_CD4T_s <- perms_CD4T_s %>% arrange(permutation, gene)

# get top_pval for each permutation in wide format
perms_CD4T_s <- perms_CD4T_s %>%
  select(gene, permutation, top_pval) %>%
  pivot_wider(
    names_from = permutation,
    values_from = top_pval,
    names_prefix = 'perm_')

### COMPILE MATRIXEQTL RESULTS
# NK
# select top SNP for each gene in non permuted file
best_NK_m <- fread('HALEYs/matrixEQTL_results/NI_NK_adj_3PCs_cisQTL_sumstats.txt') %>% group_by(gene) %>% 
  slice_min(pvalue, with_ties=FALSE) %>% arrange(gene)

# select top SNP for each gene in permuted files
for (perm in seq(1:5)){
  print(perm)
  tmp <- fread('HALEYs/matrixEQTL_results/NI_NK_Perm'%&%as.character(perm)%&%'_adj_3PCs_cisQTL_sumstats.txt') %>% 
    group_by(gene) %>% slice_min(pvalue, with_ties=FALSE) %>% ungroup() %>% 
    select(gene, pvalue) %>% arrange(gene)
  
  if (exists('perms_NK_m')){
    perms_NK_m <- inner_join(perms_NK_m, tmp, by=c('gene'='gene'))
  } else {perms_NK_m <- tmp}
}

# CD4 T
# select top SNP for each gene in non permuted file
best_CD4T_m <- fread('HALEYs/matrixEQTL_results/NI_T-CD4_adj_4PCs_cisQTL_sumstats.txtt') %>% group_by(gene) %>% 
  slice_min(pvalue, with_ties=FALSE) %>% arrange(gene)

# select top SNP for each gene in permuted files
for (perm in seq(1:5)){
  print(perm)
  tmp <- fread('HALEYs/matrixEQTL_results/NI_T-CD4_Perm'%&%as.character(perm)%&%'_adj_4PCs_cisQTL_sumstats.txt') %>% 
    group_by(gene) %>% slice_min(pvalue, with_ties=FALSE) %>% ungroup() %>% 
    select(gene, pvalue) %>% arrange(gene)
  
  if (exists('perms_CD4T_m')){
    perms_CD4T_m <- inner_join(perms_CD4T_m, tmp, by=c('gene'='gene'))
  } else {perms_CD4T_m <- tmp}
}

### compute qvalues for NK
nk_intersection <- intersect(best_NK_s$gene, best_NK_m$gene)
best_NK_s <- best_NK_s %>% filter(gene %in% nk_intersection)
perms_NK_s <- perms_NK_s %>% filter(gene %in% nk_intersection)
best_NK_m <- best_NK_m %>% filter(gene %in% nk_intersection)
perms_NK_m <- perms_NK_m %>% filter(gene %in% nk_intersection)

empP <- empPvals(stat=-log10(best_NK_s$top_pval), stat0=-log10(as.matrix(perms_NK_s[,2:6])), pool=TRUE)
best_NK_s$qvals <- qvalue(empP)$qvalue
fwrite(best_NK_s, 'Saige/step3/outputs/NK_NI_best_eQTLs_no_perm_withqvalues.txt', sep=' ')

empP <- empPvals(stat=-log10(best_NK_m$pvalue), stat0=-log10(as.matrix(perms_NK_m[,2:6])), pool=TRUE)
best_NK_m$qvals <- qvalue(empP)$qvalue
fwrite(best_NK_m, 'HALEYs/matrixEQTL_results/NI_NK_adj_3PCs_best_cisQTL_withqvalue_sumstats.txt', sep=' ')

### compute qvalues for CD4T
cd4t_intersection <- intersect(best_CD4T_s$gene, best_CD4T_m$gene)
best_CD4T_s <- best_CD4T_s %>% filter(gene %in% cd4t_intersection)
perms_CD4T_s <- perms_CD4T_s %>% filter(gene %in% cd4t_intersection)
best_CD4T_m <- best_CD4T_m %>% filter(gene %in% cd4t_intersection)
perms_CD4T_m <- perms_CD4T_m %>% filter(gene %in% cd4t_intersection)

empP <- empPvals(stat=-log10(best_CD4T_s$top_pval), stat0=-log10(as.matrix(perms_CD4T_s[,2:6])), pool=TRUE)
best_CD4T_s$qvals <- qvalue(empP)$qvalue
fwrite(best_CD4T_s, 'Saige/step3/outputs/CD4T_NI_best_eQTLs_no_perm_withqvalues.txt', sep=' ')

empP <- empPvals(stat=-log10(best_CD4T_m$pvalue), stat0=-log10(as.matrix(perms_CD4T_m[,2:6])), pool=TRUE)
best_CD4T_m$qvals <- qvalue(empP)$qvalue
fwrite(best_CD4T_m, 'HALEYs/matrixEQTL_results/NI_CD4T_adj_4PCs_best_cisQTL_withqvalue_sumstats.txt', sep=' ')

### summarise and compare 
NK_joint <- inner_join(best_NK_m, best_NK_s, by=c('gene'))
NK_joint_summary <- data.frame('matrix'=sum(NK_joint$qvals.x<0.05),
                               'saigeqtl'=sum(NK_joint$qvals.y<0.05))

CD4T_joint <- inner_join(best_CD4T_m, best_CD4T_s, by=c('gene'))
CD4T_joint_summary <- data.frame('matrix'=sum(CD4T_joint$qvals.x<0.05),
                                 'saigeqtl'=sum(CD4T_joint$qvals.y<0.05))