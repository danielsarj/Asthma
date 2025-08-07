library(tidyverse)
library(data.table)
library(reshape2)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('asthma','asthma_alb','income', 'none')

# load pseudobulk object
obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# merge metadata
mdata <- obj@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
obj@meta.data <- mdata

# read results
for (int in interactions){
  print(int)
  if (int=='none'){
    DE_results <- fread('NI_IVAxRV_limma_results_avglogCPM.filtered.txt') %>% rename(gene=Gene)
  } else {
    DE_results <- fread('NI_IVAxRV_'%&%int%&%'_limma_results_avglogCPM.filtered.txt')
  }
  
  # find top 5 genes
  top_genes <- DE_results %>% group_by(celltype, condition) %>% slice_min(P.Value, n=5, with_ties=F) %>% 
    select(gene, condition, celltype, logFC, P.Value) %>% ungroup()
  
  # extract count matrix and keep only top genes
  count <- obj@assays$RNA$counts
  count <- count[rownames(count) %in% top_genes$gene,] %>% as.data.frame()

  for (cond in unique(top_genes$condition)){
    for (ct in unique(top_genes$celltype)){
      genes_to_look <- top_genes %>% filter(celltype==ct, condition==cond) %>% pull(gene)
      sub_count <- count %>% select(contains(cond%&%'_'%&%ct), contains('NI_'%&%ct)) %>% rownames_to_column('gene') %>%
        filter(gene %in% genes_to_look) %>% melt() %>% separate(variable, into=c('ID', 'condition', 'celltype'), sep='_')
      
      # merge with metadata
      sub_count <- inner_join(sub_count, sample_m, by=c('ID'='ID'))
      
      # make boxplot
      if (int=='none'){
        ggplot(sub_count, aes(x=condition, y=value)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group=condition), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')
        
      } else if (int=='asthma'){
        ggplot(sub_count, aes(x=condition, y=value, color=asthma)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')   
          
      } else if (int=='asthma_alb'){
        ggplot(sub_count, aes(x=condition, y=value, color=asthma)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')         
        
      } else {
        ggplot(sub_count, aes(x=condition, y=value, color=income)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')
        
      }
    }
  }
}
