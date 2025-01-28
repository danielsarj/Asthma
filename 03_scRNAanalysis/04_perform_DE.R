library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('RV', 'IVA')

for (i in 1:length(conditions)){
  print(c(conditions[i]))
  
  objs <- readRDS('NI_'%&%conditions[i]%&%'.integrated.w_celltype.rds')
  
  # clean and organize metadata
  meta.data <- objs@meta.data
  meta.data$condition <- factor(meta.data$condition, levels=c('NI', conditions[i]))
  meta.data$predicted.celltype.l1 <- gsub(' ', '_', meta.data$predicted.celltype.l1)
  meta.data$predicted.celltype.l2 <- gsub(' ', '_', meta.data$predicted.celltype.l2)
  objs@meta.data <- meta.data
  
  # aggregate cell types per individual per condition (pseudobulk)
  bulk_objs <- AggregateExpression(objs, group.by=c('IDs','predicted.celltype.l1','condition'), 
                                   slot='counts', assays='RNA', return.seurat=T)
  Idents(bulk_objs) <- 'condition'

  # celltype specific DESeq2
  for (ctype in unique(bulk_objs$predicted.celltype.l1)){
    # subset pseudobulk object
    tmp <- subset(bulk_objs, predicted.celltype.l1==ctype)
    
    # DESeq2
    de_markers <- FindMarkers(tmp, ident.1=conditions[i], ident.2='NI', slot='counts', 
                              test.use='DESeq2')
    de_markers$gene <- rownames(de_markers)
    
    # volcano plot
    ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size=0.5, alpha=0.5) +
      theme_bw() +ylab('-log10(unadjusted p-value)') + ggtitle(ctype) + 
      geom_text_repel(aes(label=ifelse(p_val_adj<0.01, gene, '')), colour='red', size=3)
    ggsave('NI_'%&%conditions[i]%&%'_DESeq2_'%&%ctype%&%'_volcanoplot.pdf', height=6, width=8)
  }
}
