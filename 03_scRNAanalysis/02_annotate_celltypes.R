library(Seurat)
library(SeuratDisk)
library(patchwork)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('NI','RV','IVA')
umaps_unfilt <- list()
umaps_filt <- list()
hists <- list()

# load reference
pbmcref <- LoadH5Seurat('azimuthrefs/pbmcref/multi.h5seurat', assays='SCT')

for (i in 1:length(conditions)){
  print(c(conditions[i]))
  
  objs <- readRDS(conditions[i] %&%'_QCintegrated.rds')

  # annotate celltypes using a PBMC reference
  print(objs)
  anchors <- FindTransferAnchors(reference=pbmcref, query=objs, normalization.method='SCT', 
                                 reference.reduction='spca', dims=1:30)
  objs <- MapQuery(anchorset=anchors, query=objs, reference=pbmcref, 
                   refdata=list(celltype.l1='celltype.l1', celltype.l2='celltype.l2'),
                     reference.reduction='spca', reduction.model='wnn.umap')
  saveRDS(objs, file=conditions[i] %&%'.integrated.w_celltype.rds')
  
  # save initial UMAP viz
  umaps_unfilt[[i]] <- DimPlot(objs, reduction='umap', group.by='predicted.celltype.l1', 
                               label=TRUE, label.size=5, repel=TRUE) + NoLegend() + 
    ggtitle('Unfiltered ' %&% conditions[i])
  
  # save histogram of cell type mapping scores
  hists[[i]] <- ggplot() + geom_histogram(aes(x=objs$predicted.celltype.l1.score)) + 
    geom_vline(xintercept=median(objs$predicted.celltype.l1.score), color='red') +
    theme_bw() + ggtitle('Cell type mapping score for ' %&% conditions[i])
  
  # filter cells based on mapping score 
  objs <- objs %>% subset(subset=predicted.celltype.l1.score>0.5)
  print(objs)

  # save filtered UMAP viz
  umaps_filt[[i]] <- DimPlot(objs, reduction='umap', group.by='predicted.celltype.l1', 
                               label=TRUE, label.size=5, repel=TRUE) + NoLegend() + 
    ggtitle('Filtered (>0.5) ' %&% conditions[i])
}
  
umaps_unfilt[[1]] | umaps_unfilt[[2]] | umaps_unfilt[[3]]
ggsave(filename='UMAP_unfilt.celltypes.allconditions.pdf', height=6, width=18)

umaps_filt[[1]] | umaps_filt[[2]] | umaps_filt[[3]]
ggsave(filename='UMAP_filt.celltypes.allconditions.pdf', height=6, width=18)

hists[[1]] | hists[[2]] | hists[[3]]
ggsave(filename='Histograms_celltypesmappingscores.allconditions.pdf', height=4, width=12)