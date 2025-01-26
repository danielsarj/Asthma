library(Seurat)
library(SeuratDisk)
library(patchwork)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('NI_RV', 'NI_IVA')

# load reference
pbmcref <- LoadH5Seurat('azimuthrefs/pbmcref/multi.h5seurat', assays='SCT')

for (i in 1:length(conditions)){
  print(c(conditions[i]))
  
  objs <- readRDS(conditions[i] %&%'_integrated.rds')

  # annotate celltypes using a PBMC reference
  anchors <- FindTransferAnchors(reference=pbmcref, query=objs, normalization.method='SCT', 
                                 reference.reduction='spca', dims=1:30)
  
  objs <- MapQuery(anchorset=anchors, query=objs, reference=pbmcref, 
                   refdata=list(celltype.l1='celltype.l1', celltype.l2='celltype.l2'),
                     reference.reduction='spca', reduction.model='wnn.umap')
  
  saveRDS(objs, file=conditions[i] %&%'.integrated.w_celltype.rds')
  
  # save UMAP viz
  DimPlot(objs, reduction='umap', group.by='predicted.celltype.l1', split.by='condition',
          label=TRUE, label.size=5, repel=TRUE)
  ggsave(filename='UMAP_'%&%conditions[i]%&%'_celltypes.pdf', height=6, width=18)
  
  # analyze metadata
  metadata <- objs@meta.data %>% 
    filter(predicted.celltype.l1 %in% c('other', 'other T')==FALSE)
  
  ## proportion of celltype per condition
  summ_condition <- metadata %>% select(condition, predicted.celltype.l1) %>%
    group_by(condition, predicted.celltype.l1) %>% summarise(n=n())
  ggplot(summ_condition) + geom_col(aes(x=condition, y=n, fill=predicted.celltype.l1), position='fill') +
    theme_bw()
  ggsave(filename='BarPlot_'%&%conditions[i]%&%'_proportion.celltypes.perCond.pdf', height=6, width=8)
  
  ## proportion of celltype per indv. and condition
  summ_indv_condition <- metadata %>% select(condition, IDs, predicted.celltype.l1) %>%
    group_by(condition, IDs, predicted.celltype.l1) %>% summarise(n=n())
  ggplot(summ_indv_condition) + geom_col(aes(x=IDs, y=n, fill=predicted.celltype.l1), position='fill') +
    facet_wrap(~condition) + theme_bw() + coord_flip()
  ggsave(filename='BarPlot_'%&%conditions[i]%&%'_proportion.celltypes.perIDandCond.pdf', height=10, width=14)
}