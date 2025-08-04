library(tidyverse)
library(data.table)
library(ggrepel)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('asthma','asthma_alb','income')

### INTEGRATE LIMMA RESULTS AND MAKE VOLCANO PLOT
for (int in interactions){
  print(int)
  for (i in 1:length(conditions)){
    print(conditions[i])
    for (ctype in cells_seurat){
      print(ctype)
    
      # read results per celltype/condition
      results <- fread('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_results_mod.txt') %>% 
        filter(condition!='NI') %>% mutate(celltype=ctype)
      
      # volcano plot
      ggplot(results) + geom_point(aes(logFC, -log10(P.Value)), size=0.5, alpha=0.5) +
        theme_bw() + ylab('-log10(pvalue)') + ggtitle(conditions[i]%&%' - '%&%ctype) +
        geom_text_repel(aes(logFC, -log10(P.Value), label=ifelse(P.Value<0.00001,gene, '')), 
                    colour='red', size=3) 
      ggsave('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_volcanoplot_mod.pdf', height=3, width=4)
    
      # combine dataframes
      if (exists('full_results')){
        full_results <- rbind(full_results, results)
      } else {full_results <- results}
    }
  }
  fwrite(full_results, 'NI_IVAxRV_'%&%int%&%'_limma_results_mod.txt', sep=' ')
  
  # histogram of pvalues
  ggplot(full_results, aes(x=P.Value)) + geom_histogram() + theme_bw() +
    facet_grid(rows=vars(condition), cols=vars(celltype), scales='free')
  ggsave('NI_IVAxRV_'%&%int%&%'_pval_histogram_mod.pdf', height=4, width=10)
  
  # volcano plot
  ggplot(full_results) + geom_point(aes(logFC, -log10(P.Value)), size=0.5, alpha=0.5) +
    theme_bw() + ylab('-log10(pvalue)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
    geom_hline(yintercept=1.30103, color='red')
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot_mod.pdf', height=5, width=8)
  rm(full_results)
}
