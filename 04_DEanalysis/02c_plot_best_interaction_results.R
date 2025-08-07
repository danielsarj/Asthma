library(tidyverse)
library(data.table)
library(edgeR)
library(limma)
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
  
  for (cond in unique(top_genes$condition)){
    for (ct in unique(top_genes$celltype)){
      # subset count object
      sub_count <- count[,(grepl(cond, colnames(count)) | grepl('NI', colnames(count))) & grepl(ct, colnames(count))] 
      genes_to_look <- top_genes %>% filter(celltype==ct, condition==cond) %>% pull(gene)
      
      # transform count into dge object
      sub_count <- DGEList(counts=sub_count)
      sub_count <- calcNormFactors(sub_count)
      sub_count <- sub_count[genes_to_look, , keep.lib.sizes=TRUE]
      
      # subset metadata
      sub_mdata <- mdata %>% rownames_to_column() %>% filter(rowname %in% rownames(sub_count$samples))
      
      # create design matrix (without effect of interest)
      if (int=='none'){
        design <- model.matrix(~IDs+batch+age+gender+n+avg_mt, data=sub_mdata)
      } else if (int=='asthma_alb'){
        design <- model.matrix(~age+gender+n+avg_mt+albuterol+condition, data=sub_mdata)
      } else {
        design <- model.matrix(~age+gender+n+avg_mt+condition, data=sub_mdata)
      }
      
      # normalize and compute precision weights
      v <- voom(sub_count, design)
      
      # access normalized expression
      voom_logCPM <- v$E %>% as.data.frame()
      
      # filter by genes of interest
      sub_count <- voom_logCPM %>% select(contains(cond%&%'_'%&%ct), contains('NI_'%&%ct)) %>% rownames_to_column('gene') %>%
        filter(gene %in% genes_to_look) %>% melt() %>% separate(variable, into=c('ID', 'condition', 'celltype'), sep='_')
      
      # merge with metadata
      sub_count <- inner_join(sub_count, sample_m, by=c('ID'='ID'))
      
      # make boxplot
      if (int=='none'){
        ggplot(sub_count, aes(x=condition, y=value)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group=condition), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')
        ggsave('boxplots/top5_DE_'%&%int%&%'_'%&%cond%&%'_'%&%ct%&%'_boxplots.pdf', height=6, width=9)
        
      } else if (int=='income'){
        ggplot(sub_count, aes(x=condition, y=value, color=income)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group = income), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')
        ggsave('boxplots/top5_DE_'%&%int%&%'_'%&%cond%&%'_'%&%ct%&%'_boxplots.pdf', height=6, width=9)
        
      } else {
        ggplot(sub_count, aes(x=condition, y=value, color=asthma)) + geom_boxplot() + theme_bw() +
          stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
          facet_wrap(~gene, scales='free')
        ggsave('boxplots/top5_DE_'%&%int%&%'_'%&%cond%&%'_'%&%ct%&%'_boxplots.pdf', height=6, width=9)
      }
    }
  }
}
