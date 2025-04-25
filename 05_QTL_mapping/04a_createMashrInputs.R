library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')

# get gene intersection
for (cond in conditions){
  for (ctype in celltypes){
    tmp <- fread('../'%&%cond%&%'_'%&%ctype%&%'_elbowPCs.txt') %>% pull(GENES)
    if (exists('shared_genes')){
      shared_genes <- intersect(shared_genes, tmp)
    } else {shared_genes <- tmp}
  }
}

# retrieve cis-snps pairs that are shared across all conditions/celltypes
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    tmp <- fread(cond%&%'_'%&%ctype%&%'_elbowPCs_cisQTL_sumstats.txt') %>% 
      filter(gene %in% shared_genes) %>% select(gene, snps)
      
    if (exists('combined_snps')){
      combined_snps <- intersect(combined_snps, tmp)
    } else {combined_snps <- tmp}
  }
}
fwrite(combined_snps, 'intersection_cisQTL.txt', sep='\t', quote=F)
rm(shared_genes, tmp)

# get the most significant (lowest qvalue) cis-snps across all conditions/celltypes
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    tmp <- fread(cond%&%'_'%&%ctype%&%'_best_cisQTL_sumstats.txt') %>% 
      select(gene, snps, qvals) 
    
    if (exists('top_pairs')){
      # combine both data frames
      top_pairs <- rbind(top_pairs, tmp)
      
      # a gene should only be in the df once (lowest qval)
      top_pairs <- top_pairs %>% group_by(gene) %>%
        slice_min(qvals, with_ties=F) %>% ungroup()
      
    } else {top_pairs <- tmp}
  }
}
rm(tmp)
top_pairs <- top_pairs %>% select(gene, snps)

# get betas/SEs for each random cis-snps pair and for each top cis-snps pair for every conditions/celltypes
combined_snps <- combined_snps %>% slice_sample(n=200000)
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    
    tmp <- fread(cond%&%'_'%&%ctype%&%'_elbowPCs_cisQTL_sumstats.txt') %>% 
      select(gene, snps, beta, SE) 
    
    tmp_random <- tmp %>% right_join(combined_snps, by=c('gene', 'snps')) %>%
      group_by(gene, snps) %>% slice_head(n=1)
    colnames(tmp_random)[3:4] <- c(cond%&%'_'%&%ctype%&%'_beta', cond%&%'_'%&%ctype%&%'_SE')
    
    tmp_strong <- tmp %>% right_join(top_pairs, by=c('gene', 'snps'))
    colnames(tmp_strong)[3:4] <- c(cond%&%'_'%&%ctype%&%'_beta', cond%&%'_'%&%ctype%&%'_SE')
    
    if (exists('random_df')){
      random_df <- left_join(random_df, tmp_random, by=c('gene', 'snps'))
    } else {random_df <- tmp_random}
    
    if (exists('strong_df')){
      strong_df <- left_join(strong_df, tmp_strong, by=c('gene', 'snps'))
    } else {strong_df <- tmp_strong}
    
  }
}
fwrite(random_df, '../mashr/mashr_in_random_df.txt', quote=F, sep='\t', na='NA')
fwrite(strong_df, '../mashr/mashr_in_strong_df.txt', quote=F, sep='\t', na='NA')
