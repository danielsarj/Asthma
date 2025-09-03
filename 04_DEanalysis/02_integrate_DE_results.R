library(tidyverse)
library(data.table)
library(ggrepel)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('none','asthma_alb','income')

# integrate DE limma results
for (int in interactions){
  print(int)
  for (i in 1:length(conditions)){
    print(conditions[i])
    for (ctype in cells_seurat){
      print(ctype)
      
      if (int=='none'){
        results <- fread('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_limma_results_wqvals.txt') %>%
          mutate(interaction='none', celltype=ctype)
      } else {
        results <- fread('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_'%&%int%&%'_limma_results_wqvals.txt') %>%
          mutate(interaction=int, celltype=ctype)
      }
    
      # combine dataframes
      if (exists('full_results')){
        full_results <- rbind(full_results, results)
      } else {full_results <- results}
    }
  }
}
full_results <- full_results %>% mutate(direction=ifelse(logFC>0, 'UP', 'DOWN'))
full_results$interaction <- gsub('asthma_alb', 'asthma', full_results$interaction)
fwrite(full_results, 'NI_IVAxRV_integrated_limma_results.txt', sep=' ')

# volcano plots of qvalues
for (int in unique(full_results$interaction)){
  tmp <- full_results %>% filter(interaction==int)
  
  ggplot(tmp) + geom_point(aes(logFC, -log10(qvals)), size=0.5, alpha=0.5) +
    theme_bw() + ylab('-log10(qvalue)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
    geom_hline(yintercept=1.30103, color='red') + ggtitle(int)
  
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot.pdf', height=4, width=8)
}

# histograms of [unadjusted] pvalues
for (int in unique(full_results$interaction)){
  tmp <- full_results %>% filter(interaction==int)
  
  ggplot(tmp, aes(x=P.Value)) + geom_histogram(binwidth=0.01, boundary=0) +
    facet_wrap(~condition+celltype, scales='free_y', ncol=length(unique(tmp$celltype))) +
    theme_bw() + ggtitle(int) + xlab('Unadjusted p-values')
  
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_pval_histogram.pdf', height=4, width=8)
}

# filter significant results based on arbitrary thresholds
filtered_results <- full_results %>% filter(abs(logFC)>0.5, qvals<0.05)

# get summary across interactions, conditions, and celltypes
summary_results <- filtered_results %>% group_by(interaction, celltype, condition, direction) %>%
  summarise(n_genes=n())

# bar plots of significant DE genes
for (int in unique(full_results$interaction)){
  tmp <- summary_results %>% filter(interaction==int)
  
  if (nrow(tmp)>0){
    ggplot(summary_results) + geom_col(aes(x=celltype, y=n_genes, fill=direction), position='dodge') +
      geom_text(aes(x=celltype, y=n_genes, label=n_genes, group=direction), position=position_dodge(width=0.9),
                vjust=-0.5, size=4) + theme_bw() + facet_wrap(~condition) + ggtitle(int)
    
    ggsave('NI_IVAxRV_'%&%int%&%'_NumOfDEgenes_barplot.pdf', height=4, width=7)
  }
}
