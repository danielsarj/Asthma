library(Seurat)
library(SeuratDisk)
library(patchwork)
library(tidyverse)
library(celldex)
library(SingleR)
library(viridis)
'%&%' <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

# load seurat object
obj <- readRDS('NI_RV_IVA_integrated.rds') %>% JoinLayers()
obj@meta.data$IDs <- gsub('SEA3', 'SEA-3', obj@meta.data$IDs)

# find cluster-specific markers
Idents(obj) <- 'seurat_clusters' 
markers_seurat <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.10)

# find top markers for seurat clusters
top_seurat <- markers_seurat %>% group_by(cluster) %>% slice_min(p_val_adj, n=100) %>% slice_max(avg_log2FC, n=20)
DimPlot(obj, reduction='rna.umap', label=T, repel=T, shuffle=T)
ggsave(filename='UMAP_NI_IVA_RV_seurat_clusters.pdf', height=6, width=8)

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

# define gene markers 
HannahMarkerGenes <- c(
  'CD34', # Hematopoietic stem cells (HSPCs)
  'GATA1','HBB', 'PPBP', # Erythrocytes or Megakaryocytes
  'PAX5', 'RAG1', 'MS4A1', 'EBF1', 'MME', 'CD19', 'CD38', # B cells
  'CD14', 'FCGR3A', 'CIITA', 'LYZ', # Monocytes
  'NCAM1', 'KLRB1', 'NKG7', 'GZMA', 'GZMB', 'PRF1', 'EOMES', 'TBX21', 'KIR3DL1', 'KIR2DL1', # NK cells
  'CD3D', 'IL7R', 'CD8A', 'CCR7', 'S100A4', 'CTLA4', # without CD8, CD4+ T cell
  'TRGV4', # γδ T
  'FOXP3', # Treg
  'MR1', 'TRBV2', # MAIT (TRBV13 may not exist in human, check below) #MR1 also expressed in Neutrophils
  'TRAC', 'TRBC1', 'TRBC2', # TCR (instead of TRA/TRB shorthand)
  'LILRA4','IRF8', # Plasmacytoid dendritic cells (pDCs)
  'CD1C', # Myeloid dendritic cells (mDCs),
  'CEACAM8' # Neutrophils
)
#VlnPlot(obj, features=HannahMarkerGenes[1])
labels <- data.frame(seurat_clusters=factor(seq(from=0, to=35)), 
                     hannah=c('CD4+ T cell','CD4+ T cell','CD4+ T cell','CD4+ T cell','CD8+ T cell',
                              'CD8+ T cell','NK cell','B cell','CD4+ T cell','CD8+ T cell','CD8+ T cell',
                              'CD4+ T cell','NK-like CD8+ T cell','B cell','CD4+ T cell','NK-like CD8+ T cell',
                              'CD4+ T cell','Treg','CD4+ T cell','Monocyte','Monocyte','CD8+ T cell','CD4+ T cell',
                              'CD4+ T cell','NK cell','Monocyte','B cell','NK cell','NK cell','CD8+ T cell','CD4+ T cell',
                              'CD8+ T cell','B cell','NK-like CD8+ T cell','CD4+ T cell','NK cell'))
obj@meta.data <- left_join(obj@meta.data, labels, by=c('seurat_clusters'))

# fix metadata missing cell IDs as row names
rownames(obj@meta.data) <- rownames(obj@reductions$rna.umap@cell.embeddings)
OneK1KMarkerGenes <- c(
  'CST3', 'FCER1A', 'SERPINF1', # Dendritic cell 
  'CD16', # Non-Classical Monocyte
  'CD14', 'LYZ', # Classical Monocyte
  'TNFRSF17', 'IGJ', # Plasma cell 
  'MS4A1', # Shared B cell
  'CD27', 'TNFRSF13B', # Memory B cell 
  'TCL1A', 'FCER2', 'IL4R', # Immature and Naïve B cell 
  'NCAM1', # Shared NK cell 
  'GZMA', 'GZMB', # Natural killer cell 
  'GZMK', 'XCL1', 'XCL2', # Natural killer cell Recruiting 
  'CD3D', # Shared T cells
  'CD8A', # Shared CD8+ T cells
  'S100B', 'KLRB1', 'LTB', 'IL7R', 'GZMK', # CD8+ S100B T cell
  'GNLY', 'NKG7', 'KLRB1', # CD8+ Effector memory T cell 
  'LTB', 'CCR7', 'PASK', # CD8+ Naïve and Central memory T cell 
  'CD4', # Shared CD4+ T cell
  'SOX4', 'ID2', 'SELL', # SOX4 CD4+ T cells
  'KLRB1', 'GZMK', 'TNFSF13B', 'IL7R', # CD4+ Effector memory and central memory T cell 
  'CCR7', 'SELL', 'LRRN3' # CD4+ Naïve and Central Memory T cell 
)
#VlnPlot(obj, features=OneK1KMarkerGenes[1])
labels <- data.frame(seurat_clusters=factor(seq(from=0, to=35)), 
                     onek1k=c('CD4+ Naïve and Central Memory T cells','CD4+ Naïve and Central Memory T cells',
                              'CD4+ Naïve and Central Memory T cells','Double-negative T cells (gamma delta?)',
                              'CD8+ Naïve and Central memory T cells','CD8+ Naïve and Central memory T cells',
                              'NK-like CD8+ Effector memory T cells','Immature and naïve B cells','Double-negative T cells (gamma delta?)',
                              'CD8+ Naïve and Central memory T cells','CD8+ Naïve and Central memory T cells','Double-negative T cells (gamma delta?)',
                              'CD8+ Effector memory T cells','Immature and naïve B cells','CD4+ Naïve and Central Memory T cells',
                              'CD8+ Effector memory T cells','Double-negative T cells (gamma delta?)','CD4+ Naïve and Central Memory T cells',
                              'Pro-thymocytes? ','Classical monocytes','Memory B cells','CD8+ Naïve and Central memory T cells',
                              'Double-negative T cells (gamma delta?)','Classical monocytes','Recruiting NK cells','Classical monocytes',
                              'Immature and naïve B cells','NK cells','NK cells','CD8+ Naïve and Central memory T cells',
                              'Double-negative T cells (gamma delta?)','CD8+ Effector memory or Naïve and Central memory T cells','B cells',
                              'CD8+ Effector memory T cells','CD4+ Naïve and Central Memory T cells','Gamma delta?'))
obj@meta.data <- left_join(obj@meta.data, labels, by=c('seurat_clusters'))

# fix metadata missing cell IDs as row names
rownames(obj@meta.data) <- rownames(obj@reductions$rna.umap@cell.embeddings)
# make UMAPs per annotation source
for (i in seq(length(ref_list)+2)){
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
  } else if (i == 5){
    DimPlot(obj, reduction='rna.umap', group.by='novershtern_ct', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_novershtern_celltypes.pdf', height=6, width=8)
  } else if (i == 6){
    DimPlot(obj, reduction='rna.umap', group.by='hannah', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_hannah_celltypes.pdf', height=6, width=8)
  } else {
    DimPlot(obj, reduction='rna.umap', group.by='onek1k', label=T, repel=T, shuffle=T) + NoLegend()
    ggsave(filename='UMAP_NI_IVA_RV_onek1k_celltypes.pdf', height=6, width=8)
  }
}

# FINAL CLUSTER ANNOTATION
labels <- data.frame(
  seurat_clusters = factor(0:35, levels = 0:35),
  celltype = c(
    'CD4-T', #0 - consensus btwn all refs
    'CD4-T', #1 - consensus btwn all refs
    'CD4-T', #2 - consensus btwn all refs
    'CD4-T', #3 - low CD4, but clusters close to CD4+ and expresses Tcell markers
    'CD8-T', #4 - CD8A+/CD3D+/CCR7+/TRAC+
    'CD8-T', #5 - CD8A+/CD3D+/CCR7+/TRAC+
    'NK', #6 - CD3D-/TRAC-/KLRB1+/NKG7+
    'B', #7 - consensus btwn all refs
    'CD4-T', #8 - low CD4, but clusters close to CD4+ and expresses Tcell markers
    'CD8-T', #9 - CD8A+/CD3D+/CCR7+/TRAC+
    NA, #10 - TRDV1+/TRAC+, likely a cluster of mixed CD8+ and gamma-delta T cells
    'CD4-T', #11 - low CD4, but clusters close to CD4+ and expresses Tcell markers
    'CD8-T', #12 CD3D+/CD8A+/NCAM1-/FCGR3A- 
    'B', #13 - consensus btwn all refs
    'CD4-T', #14 CD4+/FOXP3- 
    'MAIT', #15 CD3D+/CD8A+/NCAM1-/TRAC+/CXCR6+
    'DN-T', #16 CD3D+/TRAC+/CD8A-/CD4-/FOXP3- 
    'CD4-Treg', #17 CD3D+/CD4+/FOXP3+
    'Lymphoid_progenitors', #18 TRAC-/TRBC1-/CD3D-/FOXP3-/CCR7+/LTB+/IL7R+
    'Mono', #19 - consensus btwn all refs
    'B', #20 - consensus btwn all refs
    NA, #21 CD3D+/TRAC+/CD8A+/CD4-/FOXP3-/IRF4+/lowFOS+; looks like a cluster of CD8+ and cells
    # that look like DN-T and progenitors (based on UMAP )
    NA, #22 although CD8A- and low CD4+, it clusters with CD8+ and CD4+ cells (including Tregs)
    # likely a functional / transcriptional program cluster, not a strict lineage cluster
    'Immature_neutrophils', #23 CD3D-/S100A8+/S100A9+/MPO-/CD47–/CXCR4–/IL1B+/CXCL8+
    'NK', #24 CD3D-/CCR7-/CD8A-/KLRB1+/NKG7+
    'Mono', #25 - consensus btwn all refs
    'B', #26 - consensus btwn all refs
    'NK', #27 CD3D-/CCR7-/CD8A+/KLRB1+/NKG7+/GNLY+/GZMB+/FCGR3A+
    'NK', #28 - consensus btwn all refs
    'CD8-T', #29 CD3D+/CD4-/CD8A+/FOXP3-/TRAC+/TRDV1-
    NA, #30 - TRDV1+/TRAC+, likely a cluster of mixed CD8+ and gamma-delta T cells
    NA, #31 - TRDV1+/TRAC+, likely a cluster of mixed CD8+ and gamma-delta T cells
    'B', #32 - consensus btwn all refs
    'GD-T', #33 CD3D+/CD4-/low CD8A+/TRAC+/TRDV1-/TRDV2+/TRGV9+/ZBTB16+/KLRB1+
    'CD4-T', #34 CD3D+/CD4+/CD8A-/FOXP3-
    NA  #35 CD3D-/CD4-/CD8A-/FOXP3-/PTPRC+/MS4A1-/NKG7-/TGFBR3+.. 
    # rare stromal/endothelial/circulating non-hematopoietic contaminants?
  ),
  stringsAsFactors = FALSE
)
obj@meta.data <- left_join(obj@meta.data, labels, by=c('seurat_clusters'))

# fix metadata missing cell IDs as row names
rownames(obj@meta.data) <- rownames(obj@reductions$rna.umap@cell.embeddings)
DimPlot(obj, reduction='rna.umap', group.by='celltype', label=T, repel=T, shuffle=T) + NoLegend()
ggsave(filename='UMAP_NI_IVA_RV_manual_celltypes.pdf', height=6, width=8) 

# remove non-annotated clusters
meta_df <- obj@meta.data
filtered_meta <- meta_df %>% filter(celltype %in% c('B','CD4-T','CD8-T','Mono','NK'))
matching_cells <- rownames(filtered_meta)
obj <- subset(obj, cells=matching_cells)
obj$condition <- factor(obj$condition, levels=c('NI', 'IVA', 'RV'))

# visualize new UMAP
DimPlot(obj, reduction='rna.umap', group.by='celltype', label=T, repel=T, shuffle=T) + NoLegend()
ggsave(filename='UMAP_NI_IVA_RV_manual_filtered_celltypes.pdf', height=6, width=8) 

DimPlot(obj, reduction='rna.umap', group.by='celltype', split.by='condition', 
        label=TRUE, label.size=5, repel=TRUE)
ggsave(filename='UMAP_NI_IVA_RV_celltypes.pdf', height=4, width=9)
ggsave(filename='UMAP_NI_IVA_RV_celltypes.png', height=4, width=9)

# analyze meta.data
## proportion of celltype per condition
summ_condition <- filtered_meta %>% select(condition, celltype) %>%
  group_by(condition, celltype) %>% summarise(n=n())
summ_condition$condition <- factor(summ_condition$condition, levels=c('NI', 'IVA', 'RV'))
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

# remove unnecessary columns from metadata
meta_df <- obj@meta.data %>% select(-contains('_ct'), -contains('RNA_snn'), -hannah, -onek1k)
obj@meta.data <- meta_df
obj$condition <- factor(obj$condition, levels=c('NI', 'IVA', 'RV'))

# save object
saveRDS(obj, file='NI_IVA_RV.integrated.w_celltype_new.rds')

# get pseudobulk counts per individual per condition 
bulk_obj <- AggregateExpression(obj, group.by=c('IDs','condition','celltype'), 
                                slot='counts', assays='RNA', return.seurat=T)

# add number of cells, proportion, and avg_mt to pseudobulk metadata
summ_indv_condition <- inner_join(bulk_obj@meta.data, summ_indv_condition, 
                                  by=c('condition', 'IDs', 'celltype')) %>% select(-c(batch_ID))
ct_prop <- summ_indv_condition %>% select(IDs, condition, celltype, n) %>%
  group_by(IDs, condition) %>% mutate(total_cells=sum(n), prop=n/total_cells) %>%
  ungroup() %>% select(-total_cells)
summ_indv_condition <- left_join(summ_indv_condition, ct_prop)
rownames(summ_indv_condition) <- summ_indv_condition$orig.ident

bulk_obj@meta.data <- summ_indv_condition

# save object
saveRDS(bulk_obj, file='NI_IVA_RV.integrated.pseudobulks_new.rds')