library(tidyverse)
library(data.table)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/')

# compile best snp for all saige-qtl outputs
for (perm in c('no_perm', 'perm1', 'perm2', 'perm3', 'perm4', 'perm5')){
  print(perm)
  if (perm=='no_perm'){
    saige_files <- list.files('NK_NI/perms/no_perm/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) 
      
      # combine into a sigle dataframe
      if (exists('best_NK')){
        best_NK <- rbind(best_NK, tmp)
      } else {best_NK <- tmp}
    }
  } else {
    saige_files <- list.files('NK_NI/perms/'%&%perm%&%'/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) %>% mutate(permutation=perm)
      
      # combine into a sigle dataframe
      if (exists('perms_NK')){
        perms_NK <- rbind(perms_NK, tmp)
      } else {perms_NK <- tmp}
    }
  }
}
#best_NK <- fread('NK_NI_best_eQTLs_no_perm.txt')
#perms_NK <- fread('NK_NI_best_eQTLs_perm.txt')
best_NK <- best_NK %>% arrange(gene)
perms_NK <- perms_NK %>% arrange(permutation, gene)

# get top_pval for each permutation in wide format
perms_NK <- perms_NK %>%
  select(gene, permutation, top_pval) %>%
  pivot_wider(
    names_from = permutation,
    values_from = top_pval,
    names_prefix = 'perm_')

# compute qvalues with the top SNPs
empP <- empPvals(stat=-log10(best_NK$top_pval), stat0=-log10(as.matrix(perms_NK[,2:6])), pool=TRUE)
best_NK$qvals <- qvalue(empP)$qvalue

fwrite(best_NK, 'NK_NI_best_eQTLs_no_perm_withqvalues.txt', sep=' ')


####


# compile best snp for all saige-qtl outputs
for (perm in c('no_perm', 'perm1', 'perm2', 'perm3', 'perm4', 'perm5')){
  print(perm)
  if (perm=='no_perm'){
    saige_files <- list.files('CD4_T_NI/perms/no_perm/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) 
      
      # combine into a sigle dataframe
      if (exists('best_CD4T')){
        best_CD4T <- rbind(best_CD4T, tmp)
      } else {best_CD4T <- tmp}
    }
  } else {
    saige_files <- list.files('CD4_T_NI/perms/'%&%perm%&%'/', full.names=T)
    for (f in saige_files){
      tmp <- fread(f) %>% mutate(permutation=perm)
      
      # combine into a sigle dataframe
      if (exists('perms_CD4T')){
        perms_CD4T <- rbind(perms_CD4T, tmp)
      } else {perms_CD4T <- tmp}
    }
  }
}
#best_CD4T <- fread('CD4T_NI_best_eQTLs_no_perm.txt')
#perms_CD4T <- fread('CD4T_NI_best_eQTLs_perm.txt')
best_CD4T <- best_CD4T %>% arrange(gene)
perms_CD4T <- perms_CD4T %>% arrange(permutation, gene)

# get top_pval for each permutation in wide format
perms_CD4T <- perms_CD4T %>%
  select(gene, permutation, top_pval) %>%
  pivot_wider(
    names_from = permutation,
    values_from = top_pval,
    names_prefix = 'perm_')

# compute qvalues with the top SNPs
empP <- empPvals(stat=-log10(best_CD4T$top_pval), stat0=-log10(as.matrix(perms_CD4T[,2:6])), pool=TRUE)
best_CD4T$qvals <- qvalue(empP)$qvalue

fwrite(best_CD4T, 'CD4T_NI_best_eQTLs_no_perm_withqvalues.txt', sep=' ')