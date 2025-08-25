library(tidyverse)
library(data.table)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/HALEYs/matrixEQTL_results')

# select top SNP for each gene in non permuted file
best_NK <- fread('NI_NK_adj_3PCs_cisQTL_sumstats.txt') %>% group_by(gene) %>% 
  slice_min(pvalue, with_ties=FALSE) %>% arrange(gene)

# select top SNP for each gene in permuted files
for (perm in seq(1:5)){
  print(perm)
  tmp <- fread('NI_NK_Perm'%&%as.character(perm)%&%'_adj_3PCs_cisQTL_sumstats.txt') %>% 
    group_by(gene) %>% slice_min(pvalue, with_ties=FALSE) %>% ungroup() %>% 
    select(gene, pvalue) %>% arrange(gene)

  if (exists('compiled.perm.NK')){
    compiled.perm.NK <- inner_join(compiled.perm.NK, tmp, by=c('gene'='gene'))
  } else {compiled.perm.NK <- tmp}
}

# compute qvalues with the top SNPs
empP <- empPvals(stat=-log10(best_NK$pvalue), stat0=-log10(as.matrix(compiled.perm.NK[,2:6])), pool=TRUE)
best_NK$qvals <- qvalue(empP)$qvalue
fwrite(best_NK, 'NI_NK_adj_3PCs_best_cisQTL_withqvalue_sumstats.txt', sep=' ')

# select top SNP for each gene in non permuted file
best_TCD4 <- fread('NI_T-CD4_adj_4PCs_cisQTL_sumstats.txtt') %>% group_by(gene) %>% 
  slice_min(pvalue, with_ties=FALSE) %>% arrange(gene)

# select top SNP for each gene in permuted files
for (perm in seq(1:5)){
  print(perm)
  tmp <- fread('NI_T-CD4_Perm'%&%as.character(perm)%&%'_adj_4PCs_cisQTL_sumstats.txt') %>% 
    group_by(gene) %>% slice_min(pvalue, with_ties=FALSE) %>% ungroup() %>% 
    select(gene, pvalue) %>% arrange(gene)
  
  if (exists('compiled.perm.TCD4')){
    compiled.perm.TCD4 <- inner_join(compiled.perm.TCD4, tmp, by=c('gene'='gene'))
  } else {compiled.perm.TCD4 <- tmp}
}

# compute qvalues with the top SNPs
empP <- empPvals(stat=-log10(best_TCD4$pvalue), stat0=-log10(as.matrix(compiled.perm.TCD4[,2:6])), pool=TRUE)
best_TCD4$qvals <- qvalue(empP)$qvalue
fwrite(best_TCD4, 'NI_T-CD4_adj_4PCs_best_cisQTL_withqvalue_sumstats.txt', sep=' ')
