library(Seurat)
library(SeuratDisk)
library(patchwork)
library(tidyverse)
library(celldex)
library(SingleR)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

# load seurat object
obj <- readRDS('NI_RV_IVA_integrated.rds') %>% JoinLayers()
obj@meta.data$IDs <- gsub('SEA3', 'SEA-3', obj@meta.data$IDs)

# find cluster-specific markers
Idents(obj) <- 'seurat_clusters' 
markers_seurat <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.10)

# find top markers for seurat clusters
top_seurat <- markers_seurat %>% group_by(cluster) %>% slice_min(p_val_adj, n=100) %>% slice_max(avg_log2FC, n=20)
DimPlot(obj, reduction='rna.umap')

# annotate using different references using SingleR
ref_list <- list(celldex::BlueprintEncodeData(), celldex::DatabaseImmuneCellExpressionData(), 
                 celldex::HumanPrimaryCellAtlasData(), celldex::MonacoImmuneData(), celldex::NovershternHematopoieticData())
for (i in seq(length(ref_list))){
  print(i)
  prediction <- SingleR(test=GetAssayData(obj, layer='data'), ref=ref_list[[i]], labels=ref_list[[i]]$label.fine, clusters=obj$seurat_clusters)
  
  # retrieve labels and their respective scores
  labels <- data.frame(seurat_clusters=as.factor(prediction@rownames), labs=prediction$pruned.labels)
  scores <- prediction$scores %>% as.data.frame() %>% rownames_to_column(var='seurat_clusters') %>% 
    mutate(seurat_clusters=as.factor(as.numeric(seurat_clusters)-1)) %>% pivot_longer(cols=-seurat_clusters, names_to='labs', values_to='score')
  full_anno <- left_join(labels, scores, by=c('seurat_clusters', 'labs'))
  full_anno$labs <- as.factor(full_anno$labs)
    
  if (i == 1){
    colnames(full_anno) <- c('seurat_clusters', 'blueprint_ct', 'blueprint_ct_score')
    obj@meta.data <- left_join(obj@meta.data, full_anno, by=c('seurat_clusters'))
    DimPlot(obj, reduction='rna.umap', group.by='blueprint_ct')
    } else if (i == 2){
      colnames(full_anno) <- c('seurat_clusters', 'dice_ct', 'dice_ct_score')
      obj@meta.data <- left_join(obj@meta.data, full_anno, by=c('seurat_clusters'))
    } else if (i == 3){
      colnames(full_anno) <- c('seurat_clusters', 'hpca_ct', 'hpca_ct_score')
      obj@meta.data <- left_join(obj@meta.data, full_anno, by=c('seurat_clusters'))
    } else if (i == 4){
      colnames(full_anno) <- c('seurat_clusters', 'monaco_ct', 'monaco_ct_score')
      obj@meta.data <- left_join(obj@meta.data, full_anno, by=c('seurat_clusters'))
    } else {
      colnames(full_anno) <- c('seurat_clusters', 'novershtern_ct', 'novershtern_ct_score')
      obj@meta.data <- left_join(obj@meta.data, full_anno, by=c('seurat_clusters'))
    }
}
rm(prediction, labels, scores, full_anno)

# fix metadata missing cell IDs as row names
rownames(obj@meta.data) <- rownames(obj@reductions$rna.umap@cell.embeddings)
Idents(obj) <- 'seurat_clusters' 

# make UMAPs per annotation source
for (i in seq(length(ref_list))){
  if (i == 1){
    DimPlot(obj, reduction='rna.umap', group.by='blueprint_ct', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_blueprint_celltypes.pdf', height=6, width=8) 
  } else if (i == 2){
    DimPlot(obj, reduction='rna.umap', group.by='dice_ct', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_dice_celltypes.pdf', height=6, width=8)
  } else if (i == 3){
    DimPlot(obj, reduction='rna.umap', group.by='hpca_ct', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_hpca_celltypes.pdf', height=6, width=8)
  } else if (i == 4){
    DimPlot(obj, reduction='rna.umap', group.by='monaco_ct', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_monaco_celltypes.pdf', height=6, width=8)
  } else {
    DimPlot(obj, reduction='rna.umap', group.by='novershtern_ct', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_novershtern_celltypes.pdf', height=6, width=8)
  }
}

# remove non-annotated clusters
meta_df <- obj@meta.data
filtered_meta <- meta_df %>% filter(!is.na(celltype))
matching_cells <- rownames(filtered_meta)
obj <- subset(obj, cells=matching_cells)
obj$condition <- factor(obj$condition, levels=c('NI', 'IVA', 'RV'))

# visualize new UMAP
DimPlot(obj, reduction='rna.umap', group.by='celltype', split.by='condition', 
        label=TRUE, label.size=5, repel=TRUE)
ggsave(filename='UMAP_NI_IVA_RV_celltypes.pdf', height=4, width=9)
ggsave(filename='UMAP_NI_IVA_RV_celltypes.png', height=4, width=9)

# analyze meta.data
## proportion of celltype per condition
summ_condition <- filtered_meta %>% select(condition, celltype) %>%
  group_by(condition, celltype) %>% summarise(n=n())
(ggplot(summ_condition) + geom_col(aes(x=condition, y=n, fill=celltype), position='dodge') +
  theme_bw() + guides(fill='none')) +
(ggplot(summ_condition) + geom_col(aes(x=condition, y=n, fill=celltype), position='fill') +
  theme_bw() + ylab('p'))
ggsave(filename='BarPlot_NI_IVA_RV_absolute.proportion.celltypes.perCond.pdf', height=4, width=9)
  
## proportion of celltype per indv. and condition
summ_indv_condition <- filtered_meta %>% select(batch, condition, IDs, celltype, percent.mt) %>%
  group_by(batch, condition, IDs, celltype) %>% summarise(n=n(), avg_mt=mean(percent.mt))
summ_indv_condition <- summ_indv_condition %>%
  mutate(batch_ID = paste0(batch, '_', IDs))
summ_indv_condition$condition <- factor(summ_indv_condition$condition, levels=c('NI', 'IVA', 'RV'))

ggplot(summ_indv_condition) + geom_col(aes(x=batch_ID, y=n, fill=celltype, group=batch), position='fill') +
  facet_wrap(~condition) + theme_bw() + coord_flip() + xlab(NULL)
ggsave(filename='BarPlot_NI_IVA_RV_proportion.celltypes.perIDandCond.pdf', height=5, width=7)

ggplot(summ_indv_condition) + geom_col(aes(x=batch_ID, y=n, fill=celltype), position='stack') +
  facet_wrap(~condition) + theme_bw() + coord_flip() + xlab(NULL)
ggsave(filename='BarPlot_NI_IVA_RV_absolute.celltypes.perIDandCond.pdf', height=5, width=7)

(ggplot(summ_indv_condition) + geom_col(aes(x=batch_ID, y=n, fill=celltype), position='stack') +
    facet_wrap(~condition) + theme_bw() + coord_flip() + xlab(NULL) + guides(fill='none')) +
(ggplot(summ_indv_condition) + geom_col(aes(x=batch_ID, y=n, fill=celltype, group=batch), position='fill') +
    facet_wrap(~condition) + theme_bw() + coord_flip() + ylab('p') + xlab(NULL) +
   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()))
ggsave(filename='BarPlot_NI_IVA_RV_absolute.proportion.celltypes.perCond.pdf', height=5, width=10)

# save object
saveRDS(obj, file='NI_IVA_RV.integrated.w_celltype.rds')

# get pseudobulk counts per individual per condition 
bulk_obj <- AggregateExpression(obj, group.by=c('IDs','condition','celltype'), 
                                slot='counts', assays='RNA', return.seurat=T)

# add number of cells to pseudobulk metadata
summ_indv_condition <- inner_join(bulk_obj@meta.data, summ_indv_condition, 
                                  by=c('condition', 'IDs', 'celltype')) %>% select(-c(batch_ID))
rownames(summ_indv_condition) <- summ_indv_condition$orig.ident
bulk_obj@meta.data <- summ_indv_condition

# save object
saveRDS(bulk_obj, file='NI_IVA_RV.integrated.pseudobulks.rds')
