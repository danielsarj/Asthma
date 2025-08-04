library(Seurat)
library(Hmisc)
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

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt') %>% pull('hgnc_symbol')

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
  geom_abline(slope=1, color='red') + theme_bw() + facet_wrap(~celltype) + scale_color_viridis()
ggsave('IVAxRV_limma_logFCcorrelation_scatterplot.pdf', height=5, width=6)
rm(results, long_results, iva_results, rv_results)

colnames(full_results)[1] <- c('Gene')
### COMPUTE AVG LOGCPM FOR EACH PSEUDOBULK, AND PLOT THEIR BINS/LOGFCS
bulk_objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')
for (i in 1:length(conditions)){
  boxplots <- list()
  for (c in 1:length(cells_seurat)){
    
    # get specific celltype pseudobulk count matrix
    tmp_bulk <- bulk_objs %>% subset(celltype==cells_seurat[c] & condition==conditions[i])
    tmp_count <- tmp_bulk@assays$RNA$counts
    tmp_count <- tmp_count[rownames(tmp_count) %in% annotations,]
    
    # compute average logCPM per gene
    count_data <- DGEList(tmp_count)
    count_data <- calcNormFactors(count_data)
    logCPM <- cpm(count_data, log=TRUE)
    average_logCPM <- rowMeans(logCPM) %>% as.data.frame() %>% rownames_to_column()
    colnames(average_logCPM) <- c('GENE', 'AVG_logCPM')
    
    # define 10 bins based on average logCPM
    tmp_full_results <- full_results %>% filter(celltype==cells_seurat[c], condition==conditions[i]) %>% 
      inner_join(average_logCPM, by=c('Gene'='GENE'))
    tmp_full_results$bin <- cut2(tmp_full_results$AVG_logCPM, g=10)
    
    # combine dataframes
    if (exists('full_results_w.avglogCPM')){
      full_results_w.avglogCPM <- rbind(full_results_w.avglogCPM, tmp_full_results)
    } else {full_results_w.avglogCPM <- tmp_full_results}
    
    # create boxplots
    boxplots[[c]] <- ggplot(tmp_full_results) + geom_boxplot(aes(x=bin, y=logFC)) + theme_bw() +
      geom_hline(yintercept=0, color='red') + ggtitle(conditions[i]%&%' - '%&%cells_seurat[c]) +
        xlab('logCPM bins') + coord_flip()
    ggsave('NI_'%&%conditions[i]%&%'_'%&%cells_seurat[c]%&%'_logCPMbins.logFC_boxoplot.pdf', 
           boxplots[[c]], height=4, width=5)
  }
  # plot all boxplots (per condition)
  pdf('NI_'%&%conditions[i]%&%'_logCPMbins.logFC_boxoplot.pdf', height=6, width=12)
  do.call(grid.arrange, c(boxplots[1:length(cells_seurat)], ncol=3))
  dev.off()
}
avglogCPMdf <- full_results_w.avglogCPM %>% select(Gene, celltype, condition, AVG_logCPM)
fwrite(avglogCPMdf, 'genes_avglogCPM.txt', sep=' ')
rm(bulk_objs, tmp_bulk, tmp_count, tmp_full_results, boxplots, full_results, avglogCPMdf)

# define minimum logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','T-CD4','T-CD8','Mono','NK',
                                            'B','T-CD4','T-CD8','Mono','NK'),
                                 threshold=c(6.0,1.9,0.9,3.7,5.2,
                                             5.1,0.1,1.6,3.7,5.8),
                                 condition=c(rep('IVA',5),rep('RV',5)))
for (i in 1:length(conditions)){
  for (ctype in cells_seurat){
    
    # filter results by the logCPM threshold
    tmp_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
      pull(threshold)
    tmp <- full_results_w.avglogCPM %>% filter(celltype==ctype, condition==conditions[i], 
                                               AVG_logCPM>=tmp_threshold)
    
    # compute new adjusted FDR
    tmp$adj.P.Val <- p.adjust(tmp$P.Value, method='BH', n=nrow(tmp))
    
    # volcano plot
    ggplot(tmp) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
      theme_bw() + ylab('-log10(FDR)') + ggtitle(conditions[i]%&%' - '%&%ctype) +
      geom_text_repel(aes(logFC, -log10(adj.P.Val), label=ifelse(adj.P.Val<0.00001,
                                                                Gene, '')), colour='red', size=3)
    ggsave('NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_avglogCPM.filtered_volcanoplot.pdf', height=6, width=8)

    # combine dataframes
    if (exists('full_results_avglogCPM.filtered')){
      full_results_avglogCPM.filtered <- rbind(full_results_avglogCPM.filtered, tmp)
    } else {full_results_avglogCPM.filtered <- tmp}
  }
}
rm(full_results_w.avglogCPM, tmp)

# volcano plot of all them together
ggplot(full_results_avglogCPM.filtered) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
  theme_bw() + ylab('-log10(FDR)') + facet_grid(cols=vars(celltype), rows=vars(condition)) +
  geom_hline(yintercept=1.30103, color='red')
ggsave('NI_IVAxRV_limma_facetgrid_avglogCPM.filtered_volcanoplot.pdf', height=5, width=8)
fwrite(full_results_avglogCPM.filtered, 'NI_IVAxRV_limma_results_avglogCPM.filtered.txt', sep=' ')

# histogram of pvalues
ggplot(full_results_avglogCPM.filtered, aes(x=P.Value)) + geom_histogram() + theme_bw() +
  facet_grid(rows=vars(condition), cols=vars(celltype), scales='free')
ggsave('NI_IVAxRV_pval_histogram.pdf', height=4, width=10)

### BARPLOT OF DE GENES BASED ON FDR AND LOGFC CUTOFFS
# filter results 
filtered_results <- full_results_avglogCPM.filtered %>% filter(abs(logFC)>=1, adj.P.Val<0.05)

# get summary
summary_results <- filtered_results %>% group_by(celltype, condition, direction) %>%
  summarise(n_genes=n())

# make bar plot
ggplot(summary_results) + geom_col(aes(x=celltype, y=n_genes, fill=direction), position='dodge') +
  geom_text(aes(x=celltype, y=n_genes, label=n_genes, group=direction), position=position_dodge(width=0.9),
            vjust=-0.5, size=4) + theme_bw() + facet_wrap(~condition)
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