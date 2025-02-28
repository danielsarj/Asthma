library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

# get gene intersection
for (cond in conditions){
  for (ctype in celltypes){
    tmp <- fread('../'%&%cond%&%'_'%&%ctype%&%'_rinResiduals.txt') %>% pull(GENES)
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
    tmp <- fread(cond%&%'_'%&%ctype%&%'_cisQTL_sumstats.txt.gz') %>% 
      filter(gene %in% shared_genes) %>% select(gene, snps)
      
    if (exists('combined_snps')){
      combined_snps <- intersect(combined_snps, tmp)
    } else {combined_snps <- tmp}
  }
}
rm(shared_genes, tmp)

# get the most significant (lowest FDR) cis-snps across all conditions/celltypes
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    tmp <- fread(cond%&%'_'%&%ctype%&%'_cisQTL_sumstats.txt.gz') %>% 
      select(gene, snps, FDR) %>% right_join(combined_snps, by=c('gene', 'snps')) %>% 
      group_by(gene) %>% slice_min(FDR, with_ties=F)
    
    if (exists('top_pairs')){
      top_pairs <- rbind(top_pairs, tmp) %>% group_by(gene) %>% slice_min(FDR, with_ties=F)
    } else {top_pairs <- tmp}
  }
}
rm(tmp)
top_pairs <- top_pairs %>% select(-FDR)

# get betas/SEs for random cis-snps pairs for every conditions/celltypes
combined_snps <- combined_snps %>% slice_sample(n=200000)
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    tmp <- fread(cond%&%'_'%&%ctype%&%'_cisQTL_sumstats.txt.gz') %>% 
      select(gene, snps, beta, SE) %>% right_join(combined_snps, by=c('gene', 'snps')) %>%
      group_by(gene, snps) %>% slice_head(n=1)
    colnames(tmp)[3:4] <- c(cond%&%'_'%&%ctype%&%'_beta', cond%&%'_'%&%ctype%&%'_SE')
    
    if (exists('random_df')){
      random_df <- left_join(random_df, tmp, by=c('gene', 'snps'))
    } else {random_df <- tmp}
  }
}
fwrite(random_df, '../mashr_in_random_df.txt', sep=' ', quote=F)
rm(random_df, tmp, combined_snps)
    
# finally, get betas/SEs for each top cis-snps pair for every conditions/celltypes
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    tmp <- fread(cond%&%'_'%&%ctype%&%'_cisQTL_sumstats.txt.gz') %>% 
      select(gene, snps, beta, SE) %>% right_join(top_pairs, by=c('gene', 'snps'))
    colnames(tmp)[3:4] <- c(cond%&%'_'%&%ctype%&%'_beta', cond%&%'_'%&%ctype%&%'_SE')
    
    if (exists('strong_df')){
      strong_df <- left_join(strong_df, tmp, by=c('gene', 'snps'))
    } else {strong_df <- tmp}
  }
}
fwrite(strong_df, '../mashr_in_strong_df.txt', sep=' ', quote=F)
rm(strong_df, tmp)
