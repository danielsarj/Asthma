library(Seurat)
library(Hmisc)
library(tidyverse)
library(data.table)
library(edgeR)
library(matrixStats)
library(gridExtra)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')

bulk_objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# integrate raw limma results
for (i in 1:length(conditions)){
  for (ctype in cells_seurat){
    
    # read results per celltype/condition, add metadata 
    results <- fread('NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_results.txt') %>% 
      mutate(celltype=ctype, condition=conditions[i])
    
    # combine dataframes
    if (exists('full_results')){
      full_results <- rbind(full_results, results)
    } else {full_results <- results}
  }
}
colnames(full_results)[1] <- c('Gene')

for (i in 1:length(conditions)){
  boxplots_avg <- list()
  boxplots_median <- list()
  for (c in 1:length(cells_seurat)){
    
    # get specific celltype pseudobulk count matrix
    tmp_bulk <- bulk_objs %>% subset(celltype==cells_seurat[c] & condition==conditions[i])
    tmp_count <- tmp_bulk@assays$RNA$counts
    tmp_count <- tmp_count[rownames(tmp_count) %in% annotations,]
    
    # transform count matrix into dge object
    count_data <- DGEList(tmp_count)

    # compute mean and median logCPM per gene
    logCPM <- cpm(count_data, log=TRUE)
    average_logCPM <- rowMeans(logCPM) %>% as.data.frame() %>% rownames_to_column()
    colnames(average_logCPM) <- c('Gene', 'avg_logCPM')
    median_logCPM <- rowMedians(logCPM) %>% as.data.frame() %>% rownames_to_column()
    colnames(median_logCPM) <- c('Gene', 'median_logCPM')
    
    # define 10 bins based on average/median logCPM
    tmp_full_results <- full_results %>% filter(celltype==cells_seurat[c], condition==conditions[i]) %>% 
      inner_join(average_logCPM, by=c('Gene')) %>% inner_join(median_logCPM, by=c('Gene'))
    
    tmp_full_results$avg_bin <- cut2(tmp_full_results$avg_logCPM, g=10)
    tmp_full_results$median_bin <- cut2(tmp_full_results$median_logCPM, g=10)
    
    # combine dataframes
    if (exists('full_results_w.avglogCPM')){
      full_results_w.avglogCPM <- rbind(full_results_w.avglogCPM, tmp_full_results)
    } else {full_results_w.avglogCPM <- tmp_full_results}
    
    # create boxplots
    boxplots_avg[[c]] <- ggplot(tmp_full_results) + geom_boxplot(aes(x=avg_bin, y=logFC)) + theme_bw() +
      geom_hline(yintercept=0, color='red') + ggtitle(conditions[i]%&%' - '%&%cells_seurat[c]) +
      xlab('Average logCPM bins') + coord_flip()
    ggsave('NI_'%&%conditions[i]%&%'_'%&%cells_seurat[c]%&%'_avg_logCPMbins.logFC_boxoplot.pdf', 
           boxplots_avg[[c]], height=4, width=5)
    
    boxplots_median[[c]] <- ggplot(tmp_full_results) + geom_boxplot(aes(x=median_bin, y=logFC)) + theme_bw() +
      geom_hline(yintercept=0, color='red') + ggtitle(conditions[i]%&%' - '%&%cells_seurat[c]) +
      xlab('Median logCPM bins') + coord_flip()
    ggsave('NI_'%&%conditions[i]%&%'_'%&%cells_seurat[c]%&%'_median_logCPMbins.logFC_boxoplot.pdf', 
           boxplots_median[[c]], height=4, width=5)
  }
  # plot all boxplots (per condition)
  pdf('NI_'%&%conditions[i]%&%'_avg_logCPMbins.logFC_boxoplot.pdf', height=6, width=12)
  do.call(grid.arrange, c(boxplots_avg[1:length(cells_seurat)], ncol=3))
  dev.off()
  
  pdf('NI_'%&%conditions[i]%&%'_median_logCPMbins.logFC_boxoplot.pdf', height=6, width=12)
  do.call(grid.arrange, c(boxplots_median[1:length(cells_seurat)], ncol=3))
  dev.off()
}

# define minimum average logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','T-CD4','T-CD8','Mono','NK',
                                            'B','T-CD4','T-CD8','Mono','NK'),
                                 threshold=c(4.9,1.9,1,3.4,5.6,
                                             3.5,3.6,3.1,3.4,5.6),
                                 condition=c(rep('IVA',5),rep('RV',5)))
