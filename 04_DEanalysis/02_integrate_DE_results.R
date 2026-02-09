library(tidyverse)
library(data.table)
library(ggrepel)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','CD4-T','CD8-T','Mono','NK')
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
full_results <- full_results %>% mutate(direction=ifelse(logFC>0, 'UP', 'DOWN'),
                                        sig=ifelse(qvals<0.05, TRUE, FALSE))
full_results$interaction <- gsub('asthma_alb', 'asthma', full_results$interaction)
fwrite(full_results, 'NI_IVAxRV_integrated_limma_results_new.txt', sep=' ')

# volcano plots and histograms of pvalues
for (int in unique(full_results$interaction)){
  tmp <- full_results %>% filter(interaction==int)
  
  # volcano plot
  ggplot(tmp) + geom_point(aes(logFC, -log10(P.Value), color=sig), size=0.5, alpha=0.5) +
    theme_bw() + ylab('-log10(qvalue)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
    ggtitle(int) + scale_color_manual(values=c('TRUE'='red', 'FALSE'='black')) +
    theme(legend.position='none')
  
  #ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot.pdf', height=4, width=8)
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_facetgrid_volcanoplot_new.png', height=4, width=8)
  
  # histogram 
  ggplot(tmp, aes(x=P.Value)) + geom_histogram(binwidth=0.05, boundary=0) +
    facet_wrap(~condition+celltype, scales='free_y', ncol=length(unique(tmp$celltype))) +
    theme_bw() + ggtitle(int) + xlab('Unadjusted p-values')
  
  #ggsave('NI_IVAxRV_'%&%int%&%'_limma_pval_histogram.pdf', height=4, width=8)
  ggsave('NI_IVAxRV_'%&%int%&%'_limma_pval_histogram_new.png', height=4, width=8)
  
  # compare logFCs of significant DE genes
  tmp <- tmp %>% filter(sig==TRUE) 
  if (nrow(tmp)>0){
    tmp <- tmp %>% select(Gene, logFC, condition, interaction, celltype, sig) %>%
      pivot_wider(names_from=condition, values_from=logFC) %>% drop_na()
    
    ggplot(tmp, aes(x=IVA, y=RV)) + geom_point() + theme_bw() + facet_wrap(~celltype) +
      ylab('RV LogFC') + xlab('IVA LogFC') + geom_abline(slope=1, color='red') +
      geom_hline(yintercept=0, color='blue') + geom_vline(xintercept=0, color='blue')
    
    #ggsave('NI_IVAxRV_'%&%int%&%'_ComparelogFCs_barplot.pdf', height=4, width=4)
    ggsave('NI_IVAxRV_'%&%int%&%'_ComparelogFCs_barplot_new.png', height=4, width=4)
    
  }
}

# filter significant results
filtered_results <- full_results %>% filter(sig==TRUE)

# get summary across interactions, conditions, and celltypes
summary_results <- filtered_results %>% group_by(interaction, celltype, condition, direction) %>%
  summarise(n_genes=n(), .groups='drop') %>% mutate(n_genes_signed=ifelse(direction=='DOWN', -n_genes, n_genes))

# bar plots of significant DE genes
for (int in unique(full_results$interaction)){
  tmp <- summary_results %>% filter(interaction==int)
  
  if (nrow(tmp)>0){
    ggplot(tmp, aes(x=celltype, y=n_genes_signed, fill=direction)) + geom_col() +
      geom_text(aes(label=abs(n_genes), group=direction), 
                vjust=ifelse(summary_results$n_genes_signed>0, -0.5, 1.2), size=4) + 
      theme_bw() + facet_wrap(~condition) + ggtitle(int) + geom_hline(yintercept=0, color='black') +
      labs(y='Number of DE genes (Up vs Down)') + theme(legend.position='none')
    
    #ggsave('NI_IVAxRV_'%&%int%&%'_NumOfDEgenes_barplot.pdf', height=4, width=7)
    ggsave('NI_IVAxRV_'%&%int%&%'_NumOfDEgenes_barplot_new.png', height=4, width=7)
    
  }
}
