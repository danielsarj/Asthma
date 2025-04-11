library(tidyverse)
library(data.table)
library(ggrepel)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('asthma', 'income')

# define minimum logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','T-CD4','T-CD8','Mono','NK',
                                            'B','T-CD4','T-CD8','Mono','NK'),
                                 threshold=c(6.0,1.9,0.9,3.7,5.2,
                                             5.1,0.1,1.6,3.7,5.8),
                                 condition=c(rep('IVA',5),rep('RV',5)))

# get avg logCPM per gene
logCPM <- fread('genes_avglogCPM.txt')

### INTEGRATE LIMMA RESULTS AND MAKE VOLCANO PLOT
for (int in interactions){
  for (i in 1:length(conditions)){
    for (ctype in cells_seurat){
      
      # retrieve logCPM threshold
      tmp_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
    
      # read results per celltype/condition
      results <- fread('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_results.txt') %>% 
        filter(condition!='NI')
      
      if (nrow(results)>0){
        # merge with avg logCPM values
        subset_logCPM <- logCPM %>% filter(celltype==ctype, condition==conditions[i])
        results <- inner_join(results, subset_logCPM, by=c('gene'='Gene', 'condition'))

        # filter results
        results <- results %>% filter(AVG_logCPM>=tmp_threshold) %>% 
          mutate(direction=ifelse(logFC>0, 'UP', 'DOWN'))
        
        # compute new adjusted FDR
        results$adj.P.Val <- p.adjust(results$P.Value, method='BH', n=nrow(results))

        # volcano plot
        ggplot(results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
          theme_bw() + ylab('-log10(FDR)') + ggtitle(conditions[i]%&%' - '%&%ctype) +
          geom_text_repel(aes(logFC, -log10(adj.P.Val), label=ifelse(adj.P.Val<0.05,gene, '')), 
                      colour='red', size=3) 
        ggsave('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_volcanoplot.pdf', height=3, width=4)
    
        # combine dataframes
        if (exists('full_results')){
          full_results <- rbind(full_results, results)
        } else {full_results <- results}
      }
      }
  }
  fwrite(full_results, 'NI_IVAxRV_'%&%int%&%'_limma_results_avglogCPM.filtered.txt', sep=' ')
  
  ggplot(full_results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
    theme_bw() + ylab('-log10(FDR)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
    geom_hline(yintercept=1.30103, color='red')
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot.pdf', height=5, width=8)
  rm(full_results)
}
