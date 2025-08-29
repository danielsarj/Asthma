library(tidyverse)
library(data.table)
library(ggrepel)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('asthma','asthma_alb','income')

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
  print(int)
  for (i in 1:length(conditions)){
    print(conditions[i])
    for (ctype in cells_seurat){
      print(ctype)
      
      # retrieve logCPM threshold
      tmp_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
    
      # read results per celltype/condition
      results <- fread('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_results_v2design.txt') %>% 
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
        ggplot(results) + geom_point(aes(logFC, -log10(P.Value)), size=0.5, alpha=0.5) +
          theme_bw() + ylab('-log10(pvalue)') + ggtitle(conditions[i]%&%' - '%&%ctype) +
          geom_text_repel(aes(logFC, -log10(P.Value), label=ifelse(P.Value<0.00001,gene, '')), 
                      colour='red', size=3) 
        ggsave('NI_'%&%conditions[i]%&%'_'%&%int%&%'_limma_'%&%ctype%&%'_volcanoplot_v2design.pdf', height=3, width=4)
    
        # combine dataframes
        if (exists('full_results')){
          full_results <- rbind(full_results, results)
        } else {full_results <- results}
      }
      }
  }
  fwrite(full_results, 'NI_IVAxRV_'%&%int%&%'_limma_results_avglogCPM.filtered_v2design.txt', sep=' ')
  
  # histogram of pvalues
  ggplot(full_results, aes(x=P.Value)) + geom_histogram() + theme_bw() +
    facet_grid(rows=vars(condition), cols=vars(celltype), scales='free')
  ggsave('NI_IVAxRV_'%&%int%&%'_pval_histogram_v2design.pdf', height=4, width=10)
  
  # volcano plot
  ggplot(full_results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
    theme_bw() + ylab('-log10(adjusted p-value)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
    geom_hline(yintercept=1.30103, color='red')
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot_v2design.pdf', height=5, width=8)
  
  # get summary
  summary_results <- full_results %>% filter(abs(logFC)>=1, adj.P.Val<0.05) %>% 
    group_by(celltype, condition, direction) %>% summarise(n_genes=n())
  
  if (nrow(summary_results)>0){
    # make bar plot
    ggplot(summary_results) + geom_col(aes(x=celltype, y=n_genes, fill=direction), position='dodge') +
      geom_text(aes(x=celltype, y=n_genes, label=n_genes, group=direction), position=position_dodge(width=0.9),
                vjust=-0.5, size=4) + theme_bw() + facet_wrap(~condition)
    ggsave('NI_IVAxRV_'%&%int%&%'_NumOfDEgenes_barplot_v2design.pdf', height=4, width=6)
  }
  rm(full_results)
}
