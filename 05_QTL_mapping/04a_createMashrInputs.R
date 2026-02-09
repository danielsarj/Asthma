library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')
input_prefix <- c('IVA_B_4',
                  'NI_B_5',
                  'RV_B_4',
                  'IVA_CD4-T_0',
                  'NI_CD4-T_0',
                  'RV_CD4-T_1',
                  'IVA_CD8-T_1',
                  'NI_CD8-T_0',
                  'RV_CD8-T_2',
                  'IVA_Mono_14',
                  'NI_Mono_19',
                  'RV_Mono_0',
                  'IVA_NK_0',
                  'NI_NK_2',
                  'RV_NK_0')

# get gene intersection
for (f in input_prefix){
  tmp <- fread('../'%&%f%&%'PCs_new.txt') %>% pull(GENES)
  if (exists('shared_genes')){
    shared_genes <- intersect(shared_genes, tmp)
  } else {shared_genes <- tmp}
}

# retrieve cis-snps pairs that are shared across all conditions/celltypes
for (f in input_prefix){
  tmp <- fread(f%&%'PCs_cisQTL_sumstats_new.txt') %>% 
    filter(gene %in% shared_genes) %>% select(gene, snps)
    
  if (exists('combined_snps')){
    combined_snps <- intersect(combined_snps, tmp)
  } else {combined_snps <- tmp}
}
fwrite(combined_snps, 'intersection_cisQTL.txt', sep='\t', quote=F)
rm(tmp)

# get the most significant (lowest qvalue) cis-snps across all conditions/celltypes
for (f in input_prefix){
  tmp <- fread(f%&%'_best_cisQTL_sumstats_new.txt') %>% 
    filter(gene %in% shared_genes) %>% select(gene, snps, qvals) 
    
  if (exists('top_pairs')){
    # combine both data frames
    top_pairs <- rbind(top_pairs, tmp)
      
    # a gene should only be in the df once (lowest qval)
    top_pairs <- top_pairs %>% group_by(gene) %>%
      slice_min(qvals, with_ties=F) %>% ungroup()
    
  } else {top_pairs <- tmp}
}
rm(shared_genes, tmp)
top_pairs <- top_pairs %>% select(gene, snps)

# get betas/SEs for each random cis-snps pair and for each top cis-snps pair for every conditions/celltypes
combined_snps <- combined_snps %>% slice_sample(n=200000)
for (f in input_prefix){
  base <- sub('_[^_]+$', '', f)
  
  tmp <- fread(f%&%'PCs_cisQTL_sumstats_new.txt') %>% 
    select(gene, snps, beta, SE) 
    
  tmp_random <- tmp %>% right_join(combined_snps, by=c('gene', 'snps')) 
  colnames(tmp_random)[3:4] <- c(base%&%'_beta', base%&%'_SE')

  tmp_strong <- tmp %>% right_join(top_pairs, by=c('gene', 'snps'))
  colnames(tmp_strong)[3:4] <- c(base%&%'_beta', base%&%'_SE')
    
  if (exists('random_df')){
    random_df <- left_join(random_df, tmp_random, by=c('gene', 'snps'))
  } else {random_df <- tmp_random}
  
  if (exists('strong_df')){
    strong_df <- left_join(strong_df, tmp_strong, by=c('gene', 'snps'))
  } else {strong_df <- tmp_strong}
}
random_df <- random_df %>% select(gene, snps, contains('NI'), contains('RV'), contains('IVA'))
strong_df <- strong_df %>% select(gene, snps, contains('NI'), contains('RV'), contains('IVA'))

fwrite(random_df, '../mashr/mashr_in_random_df_new.txt', quote=F, sep='\t', na='NA')
fwrite(strong_df, '../mashr/mashr_in_strong_df_new.txt', quote=F, sep='\t', na='NA')