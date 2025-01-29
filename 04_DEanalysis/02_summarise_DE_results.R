library(Seurat)
library(Hmisc)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','CD4-T','CD8-T','DC','Mono','NK')

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt') %>% pull('hgnc_symbol')

### INTEGRATE LIMMA RESULTS AND MAKE VOLCANO PLOT
for (i in 1:length(conditions)){
  for (ctype in c('B','CD4-T','CD8-T','DC','Mono','NK')){
    
    # read results per celltype/condition, add metadata 
    results <- fread('NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_results.txt') %>% 
      mutate(celltype=ctype, condition=conditions[i], direction=ifelse(logFC>0, 'UP', 'DOWN'))
  
    # volcano plot
    ggplot(results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
      theme_bw() + ylab('-log10(FDR)') + ggtitle(conditions[i]%&%' - '%&%ctype) +
      geom_text_repel(aes(logFC, -log10(adj.P.Val), label=ifelse(adj.P.Val<0.00001,V1, '')), 
                    colour='red', size=3)
    ggsave('NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_volcanoplot.pdf', height=6, width=8)
  
    # combine dataframes
    if (exists('full_results')){
      full_results <- rbind(full_results, results)
    } else {full_results <- results}
  }
}
# volcano plot of all them together
ggplot(full_results) + geom_point(aes(logFC, -log10(adj.P.Val)), size=0.5, alpha=0.5) +
  theme_bw() + ylab('-log10(adjusted p-value)') + facet_grid(cols=vars(celltype), rows=vars(condition))
ggsave('NI_IVAxRV_limma_facetgrid_volcanoplot.pdf', height=6, width=10)
colnames(full_results)[1] <- c('Gene')
rm(results)

### COMPUTE AVG LOGCPM FOR EACH PSEUDOBULK, AND PLOT THEIR BINS/LOGFCS
for (i in 1:length(conditions)){
  bulk_objs <- readRDS('NI_'%&%conditions[i]%&%'_pseudobulks.rds')
  boxplots <- list()
  
  for (c in 1:length(cells_seurat)){
    
    # get specific celltype pseudobulk count matrix
    tmp_bulk <- bulk_objs %>% subset(predicted.celltype.l1==cells_seurat[c])
    tmp_count <- tmp_bulk@assays$RNA$counts
    tmp_count <- tmp_count[rownames(tmp_count) %in% annotations,]
    
    # get average logcpm
    library_sizes <- Matrix::colSums(tmp_count)
    cpm <- t(t(tmp_count) / library_sizes) * 1e6
    log_cpm <- log2(cpm + 1)
    average_log_cpm <- Matrix::rowMeans(log_cpm) %>% as.data.frame() %>% rownames_to_column()
    colnames(average_log_cpm) <- c('GENE', 'AVG_logCPM')
    
    # define 10 bins based on average logCPM
    tmp_full_results <- full_results %>% filter(celltype==cells_seurat[c], condition==conditions[i]) %>% 
      inner_join(average_log_cpm, by=c('Gene'='GENE'))
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
           boxplots[[c]], height=4, width=8)
  }
  # plot all boxplots (per condition)
  pdf('NI_'%&%conditions[i]%&%'_logCPMbins.logFC_boxoplot.pdf', height=6, width=12)
  do.call(grid.arrange, c(boxplots[1:length(cells_seurat)], ncol=3))
  dev.off()
}
rm(bulk_objs, tmp_bulk, tmp_count, library_sizes, cpm, log_cpm, average_log_cpm, 
   tmp_full_results, boxplots, full_results)

# define minimum logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','CD4-T','CD8-T','DC','Mono','NK',
                                            'B','CD4-T','CD8-T','DC','Mono','NK'),
                                 threshold=c(1.5,0.2,3.2,2.7,0.2,1.5,
                                             1.6,0.6,1.9,0.3,0.2,3.2),
                                 condition=c(rep('IVA',6),rep('RV',6)))
for (i in 1:length(conditions)){
  for (ctype in c('B','CD4-T','CD8-T','DC','Mono','NK')){
    
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
  theme_bw() + ylab('-log10(FDR)') + facet_grid(cols=vars(celltype), rows=vars(condition))
ggsave('NI_IVAxRV_limma_facetgrid_avglogCPM.filtered_volcanoplot.pdf', height=6, width=10)

### BARPLOT OF DE GENES BASED ON FDR AND LOGFC CUTOFFS
# filter results 
filtered_results <- full_results_avglogCPM.filtered %>% filter(abs(logFC)>1, adj.P.Val<0.05)

# get summary
summary_results <- filtered_results %>% group_by(celltype, condition, direction) %>%
  summarise(n_genes=n())

# make bar plot
ggplot(summary_results) + geom_col(aes(x=celltype, y=n_genes, fill=direction), position='dodge') +
  geom_text(aes(x=celltype, y=n_genes, label=n_genes, group=direction), position=position_dodge(width=0.9),
            vjust=-0.5, size=4) + theme_bw() + facet_wrap(~condition)
ggsave('NI_IVAxRV_NumOfDEgenes_barplot.pdf', height=4, width=6)
