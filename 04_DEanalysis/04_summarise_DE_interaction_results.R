library(Hmisc)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(viridis)
library(ggpointdensity)
library(UpSetR)
library(grid)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','CD4-T','CD8-T','DC','Mono','NK')
interactions <- c('asthma', 'income')

# define minimum logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','CD4-T','CD8-T','DC','Mono','NK',
                                            'B','CD4-T','CD8-T','DC','Mono','NK'),
                                 threshold=c(1.4,0.1,2.7,2.4,0.4,1.7,
                                             2.9,1.2,1.4,0.4,0.1,3.2),
                                 condition=c(rep('IVA',6),rep('RV',6)))

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
      results <- fread('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_results.txt', fill=T) %>% 
        filter(condition!='NI') %>% drop_na()
      
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
        ggsave('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_volcanoplot.pdf', height=6, width=8)
    
        # combine dataframes
        if (exists('full_results')){
          full_results <- rbind(full_results, results)
        } else {full_results <- results}
      }
      }
    }
  ggplot(full_results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
    theme_bw() + ylab('-log10(adjusted p-value)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
    geom_hline(yintercept=1.30103, color='red')
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot.pdf', height=6, width=10)
  rm(full_results)
}
