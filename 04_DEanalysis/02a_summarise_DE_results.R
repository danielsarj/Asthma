library(Seurat)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(viridis)
library(ggpointdensity)
library(UpSetR)
library(grid)
library(edgeR)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')

### INTEGRATE LIMMA RESULTS AND MAKE VOLCANO PLOT
for (i in 1:length(conditions)){
  for (ctype in cells_seurat){
    
    # read results per celltype/condition, add metadata 
    results <- fread('NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_results.txt') %>% 
      mutate(celltype=ctype, condition=conditions[i], direction=ifelse(logFC>0, 'UP', 'DOWN'))
    
    # volcano plot
    ggplot(results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
      theme_bw() + ylab('-log10(FDR)') + ggtitle(conditions[i]%&%' - '%&%ctype) +
      geom_text_repel(aes(logFC, -log10(adj.P.Val), label=ifelse(adj.P.Val<0.00001,V1, '')), 
                      colour='red', size=3)
    ggsave('NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_volcanoplot.pdf', height=3, width=4)
    
    # combine dataframes
    if (exists('full_results')){
      full_results <- rbind(full_results, results)
    } else {full_results <- results}
  }
}
# volcano plot of all them together
ggplot(full_results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
  theme_bw() + ylab('-log10(adjusted p-value)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
  geom_hline(yintercept=1.30103, color='red')
ggsave('NI_IVAxRV_limma_facetgrid_volcanoplot.pdf', height=5, width=8)

# prepare dataframe to compare logFC between conditions
rv_results <- full_results %>% filter(condition=='RV') %>% select(V1, celltype, logFC)
colnames(rv_results)[3] <- c('logFC.RV')
iva_results <- full_results %>% filter(condition=='IVA') %>% select(V1, celltype, logFC)
colnames(iva_results)[3] <- c('logFC.IVA')
long_results <- full_join(rv_results, iva_results, by=c('V1', 'celltype'))

# scatter plot with point density 
ggplot(long_results) + geom_pointdensity(aes(logFC.RV, logFC.IVA), show.legend=F) +
  stat_smooth(aes(logFC.RV, logFC.IVA), method='lm', geom='smooth', formula=y~x) +
  geom_abline(slope=1, color='red') + theme_bw() + facet_wrap(~celltype, scales='free') + scale_color_viridis()
ggsave('IVAxRV_limma_logFCcorrelation_scatterplot.pdf', height=5, width=6)
rm(results, long_results, iva_results, rv_results)

colnames(full_results)[1] <- c('Gene')
### BARPLOT OF DE GENES BASED ON FDR AND LOGFC CUTOFFS
# filter results 
filtered_results <- full_results %>% filter(adj.P.Val<0.05)

# get summary
summary_results <- filtered_results %>% group_by(celltype, condition, direction) %>%
  summarise(n_genes=n())

# modify df
summary_results$n_genes <- ifelse(summary_results$direction=='UP', summary_results$n_genes, -summary_results$n_genes)

# make bar plot
ggplot(summary_results) + geom_col(aes(x=celltype, y=n_genes, fill=direction), position='stack') +
  geom_text(aes(x=celltype, y=n_genes, label=abs(n_genes), group=direction), size=4) + 
  theme_bw() + facet_wrap(~condition)
ggsave('NI_IVAxRV_NumOfDEgenes_barplot.pdf', height=4, width=6)

# upset plot
for (ctype in cells_seurat){
  IVA <- filtered_results %>% filter(celltype==ctype, condition=='IVA') %>% pull(Gene)
  RV <- filtered_results %>% filter(celltype==ctype, condition=='RV') %>% pull(Gene)
  
  if (length(IVA)==0 || length(RV)==0){
    next
  }
  
  l <- list(IVA=IVA, RV=RV)
  pdf('NI_IVAxRV_SharedDEgenes_'%&%ctype%&%'_upsetplot.pdf', height=3, width=4, onefile=F)
  plot.new() 
  upset(fromList(l), point.size=3.5, line.size=2, text.scale=c(1.3, 1.3, 1, 1, 2, 1.5))
  grid.text(ctype,x=0.2, y=0.95, gp=gpar(fontsize=20))
  dev.off()
}