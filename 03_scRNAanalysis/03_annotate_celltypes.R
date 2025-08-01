library(Seurat)
library(patchwork)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

# load seurat object
obj <- readRDS('NI_RV_IVA_integrated.rds') %>% JoinLayers()

# find cluster-specific markers at different resolutions
Idents(obj) <- 'seurat_clusters' 
markers_seurat <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

Idents(obj) <- 'RNA_snn_res.0.8' 
markers_res08 <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

Idents(obj) <- 'RNA_snn_res.0.4' 
markers_res04 <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

Idents(obj) <- 'RNA_snn_res.0.1' 
markers_res01 <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)

# define markers manually
cd4t <- c('CD3D', 'CD3E', 'CD4')
cd8t <- c('CD3D', 'CD8A', 'CD8B')
nk <- c('NKG7', 'GNLY', 'KLRF1')
b <- c('MS4A1', 'CD79A')
mono <- c('CD14', 'LYZ')

# make feature plot
Idents(obj) <- 'RNA_snn_res.0.4' 
FeaturePlot(obj, features=c(cd4t, cd8t, nk, b, mono), label = TRUE)
ggsave(filename='UMAP_celltypeannotation_FeaturePlot.pdf', height=15, width=15)

# make violin plot
VlnPlot(obj, features=c(cd4t, cd8t, nk, b, mono), group.by='RNA_snn_res.0.4', pt.size=0.01)
ggsave(filename='UMAP_celltypeannotation_ViolinPlot.pdf', height=15, width=20)

# add manual annotation labels to metadata
manual_anno <- data.frame(c('Mono','B','NK','T-CD4','T-CD4','T-CD4','T-CD4',
                            'T-CD8','T-CD8','T-CD8','T-CD8','T-CD8','T-CD8'),
                          c(13,11,18,0,1,3,15,2,8,10,12,14,16))
colnames(manual_anno) <- c('manual_anno', 'RNA_snn_res.0.4')
obj[['celltype']]<- manual_anno$manual_anno[match(as.vector(obj[['RNA_snn_res.0.4']])[[1]], manual_anno$RNA_snn_res.0.4)]

# remove non-annotated clusters
obj@meta.data$IDs <- gsub('SEA3', 'SEA-3', obj@meta.data$IDs)
meta_df <- obj@meta.data
filtered_meta <- meta_df %>% filter(!is.na(celltype))
matching_cells <- rownames(filtered_meta)
obj <- subset(obj, cells=matching_cells)

# visualize new UMAP
DimPlot(obj, reduction='rna.umap', group.by='celltype', split.by='condition', 
        label=TRUE, label.size=5, repel=TRUE)
ggsave(filename='UMAP_NI_IVA_RV_celltypes.pdf', height=6, width=12)

# analyze meta.data
## proportion of celltype per condition
summ_condition <- filtered_meta %>% select(condition, celltype) %>%
  group_by(condition, celltype) %>% summarise(n=n())
ggplot(summ_condition) + geom_col(aes(x=condition, y=n, fill=celltype), position='fill') +
  theme_bw()
ggsave(filename='BarPlot_NI_IVA_RV_proportion.celltypes.perCond.pdf', height=5, width=6)
  
## proportion of celltype per indv. and condition
summ_indv_condition <- filtered_meta %>% select(batch, condition, IDs, celltype, percent.mt) %>%
  group_by(batch, condition, IDs, celltype) %>% summarise(n=n(), avg_mt=mean(percent.mt))
summ_indv_condition <- summ_indv_condition %>%
  mutate(batch_ID = paste0(batch, '_', IDs))
ggplot(summ_indv_condition) + geom_col(aes(x=batch_ID, y=n, fill=celltype, group=batch), position='fill') +
  facet_wrap(~condition) + theme_bw() + coord_flip()
ggsave(filename='BarPlot_NI_IVA_RV_proportion.celltypes.perIDandCond.pdf', height=8, width=10)
ggplot(summ_indv_condition) + geom_col(aes(x=batch_ID, y=n, fill=celltype), position='stack') +
  facet_wrap(~condition) + theme_bw() + coord_flip()
ggsave(filename='BarPlot_NI_IVA_RV_absolute.celltypes.perIDandCond.pdf', height=8, width=10)

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
