library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/')
celltypes <- c('NK_NI')
pcs <- c('PCs1', 'PCs2', 'PCs3', 'PCs4', 'PCs5','PCs6','PCs7','PCs8','PCs9','PCs10', 'PCs11')
#pcs <- c('PCs11')

final_df <- fread('compiled_step3_outputs.txt')
for (ct in celltypes){
  print(ct)
  for (pc in pcs){
    print(pc)
    # find all SAIGE-QTL step3 outputs
    saige_files <- list.files(ct%&%'/exp_'%&%pc, full.names=T)
  
    for (f in saige_files){
      # read file
      tmp <- fread(f) %>% mutate(celltype=ct, exp_pcs=pc)
    
      # combine into a sigle dataframe
      if (exists('final_df')){
      final_df <- rbind(final_df, tmp)
      } else {final_df <- tmp}
    }
  }
}
fwrite(final_df, 'compiled_step3_outputs.txt', sep=' ', col.names=T)

# find eGenes with pval < 1e-5
final_df$exp_pcs <- factor(final_df$exp_pcs, levels=pcs)
summary_df <- final_df %>% group_by(celltype, exp_pcs) %>% summarise(n_eGenes=sum(top_pval<1e-5), all_genes=n(),
                                                                     prop=n_eGenes/all_genes)

# summarise 
ggplot(summary_df, aes(x=exp_pcs, y=n_eGenes, group=celltype)) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(~celltype)
ggplot(summary_df, aes(x=exp_pcs, y=prop, group=celltype)) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(~celltype)

# read Haleys results
haleys <- fread('../../../HALEYs/matrixEQTL_results/NI_NK_cisQTL_sumstats.txt') %>% group_by(gene) %>% 
  slice_min(pvalue, with_ties=F) %>% ungroup()
haleys_summary <- haleys %>% group_by(celltype) %>% summarise(n_eGenes=sum(pvalue<1e-5), all_genes=n(),
                                                                           prop=n_eGenes/all_genes)
