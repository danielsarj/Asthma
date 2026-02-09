library(Seurat)
library(Hmisc)
library(tidyverse)
library(data.table)
library(limma)
library(edgeR)
library(matrixStats)
library(gridExtra)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

# read seurat object
bulk_objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks_new.rds')

# remove batch 4
bulk_objs <- subset(bulk_objs, subset= batch!='B4')

# load sample metadata
sample_m <- fread('../sample_metadata.txt')
# merge metadata
mdata <- bulk_objs@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
mdata$batch <- factor(mdata$batch, levels=c('B1','B2','B3'))
bulk_objs@meta.data <- mdata

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# run limma once per condition/celltype
for (i in 1:length(conditions)){
  for (ctype in celltypes){
    
    # subset pseudobulk object
    tmp <- subset(bulk_objs, celltype==ctype & (condition==conditions[i] | condition=='NI'))
    
    # extract metadata and count matrices
    mdata <- tmp@meta.data
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    
    # remove non protein coding genes from count matrix and genes with variance == 0
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% annotations,]
    zero_var_genes <- apply(count, 1, var) == 0
    count <- count[!zero_var_genes, ]
    count <- DGEList(counts=count)
    count <- calcNormFactors(count)
    
    # define design matrix
    design <- model.matrix(~batch+age+gender+n+prop+avg_mt+condition, data=mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # get results
    results <- topTable(fit, coef=ncol(fit), number=Inf) %>% 
      rownames_to_column('Gene') %>% mutate(condition=conditions[i], celltype=ctype)
    
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
  for (c in 1:length(celltypes)){
    
    # get specific celltype pseudobulk count matrix
    tmp_bulk <- bulk_objs %>% subset(celltype==celltypes[c] & condition==conditions[i])
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
    median_logCPM$Gene <- average_logCPM$Gene
    
    # define 10 bins based on average/median logCPM
    tmp_full_results <- full_results %>% filter(celltype==celltypes[c], condition==conditions[i]) %>% 
      inner_join(average_logCPM, by=c('Gene')) %>% inner_join(median_logCPM, by=c('Gene'))
    tmp_full_results$avg_bin <- cut2(tmp_full_results$avg_logCPM, g=10)
    tmp_full_results$median_bin <- cut2(tmp_full_results$median_logCPM, g=10)
    
    # combine dataframes
    if (exists('full_results_w.avglogCPM')){
      full_results_w.avglogCPM <- rbind(full_results_w.avglogCPM, tmp_full_results)
    } else {full_results_w.avglogCPM <- tmp_full_results}
    
    # create boxplots
    boxplots_avg[[c]] <- ggplot(tmp_full_results) + geom_boxplot(aes(x=avg_bin, y=logFC)) + theme_bw() +
      geom_hline(yintercept=0, color='red') + ggtitle(conditions[i]%&%' - '%&%celltypes[c]) +
      xlab('Average logCPM bins') + coord_flip()
    ggsave('NI_'%&%conditions[i]%&%'_'%&%celltypes[c]%&%'_avg_logCPMbins.logFC_boxoplot.pdf', 
           boxplots_avg[[c]], height=4, width=5)
    
    boxplots_median[[c]] <- ggplot(tmp_full_results) + geom_boxplot(aes(x=median_bin, y=logFC)) + theme_bw() +
      geom_hline(yintercept=0, color='red') + ggtitle(conditions[i]%&%' - '%&%celltypes[c]) +
      xlab('Median logCPM bins') + coord_flip()
    ggsave('NI_'%&%conditions[i]%&%'_'%&%celltypes[c]%&%'_median_logCPMbins.logFC_boxoplot.pdf', 
           boxplots_median[[c]], height=4, width=5)
  }
  # plot all boxplots (per condition)
  pdf('NI_'%&%conditions[i]%&%'_avg_logCPMbins.logFC_boxoplot.pdf', height=6, width=12)
  do.call(grid.arrange, c(boxplots_avg[1:length(celltypes)], ncol=3))
  dev.off()
  
  pdf('NI_'%&%conditions[i]%&%'_median_logCPMbins.logFC_boxoplot.pdf', height=6, width=12)
  do.call(grid.arrange, c(boxplots_median[1:length(celltypes)], ncol=3))
  dev.off()
}

# define minimum average logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','CD4-T','CD8-T','Mono','NK',
                                            'B','CD4-T','CD8-T','Mono','NK'),
                                 threshold=c(2.7,-0.5,1.0,3.9,2.8,
                                             4.0,-0.7,1.9,4.0,2.8),
                                 condition=c(rep('IVA',5),rep('RV',5)))
