library(tidyverse)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(glmGamPoi)
library(future)
"%&%" <- function(a,b) paste(a,b, sep = "")
options(future.globals.maxSize=200*1024^3)
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('NI','RV','IVA')
seurat_objs <- list()

# read each file
for (c in 1:length(conditions)){
  seurat_objs <- append(seurat_objs, readRDS(conditions[c]%&%'_raw.objects.batches.rds'))
}

# merge seurat objects
merged_seurat <- merge(x=seurat_objs[[1]], y=seurat_objs[-1], 
                       merge.data=FALSE,
                       merge.dr=FALSE)
rm(seurat_objs)

# normalize each condition and regress out batch and percent.mt
merged_seurat <- merged_seurat %>% NormalizeData() %>% 
  FindVariableFeatures(selection.method='vst', nfeatures=2000) %>% 
  ScaleData(vars.to.regress=c('batch', 'percent.mt')) %>% 
  RunPCA(features=VariableFeatures(object=merged_seurat), npcs=60)

# integrate assay 
merged_seurat <- IntegrateLayers(object=merged_seurat, assay='RNA', 
                                 method=RPCAIntegration, orig.reduction='pca',
                                 new.reduction='integrated.rpca', 
                                 normalization.method='LogNormalize')

# UMAP
merged_seurat <- RunUMAP(merged_seurat, reduction='integrated.rpca', dims=1:30,
                         assay='RNA', reduction.name='rna.umap', reduction.key='rnaUMAP_')

# find clusters
merged_seurat <- merged_seurat %>% FindNeighbors(dims=1:30) %>% 
  FindClusters(resolution=c(0.1, 0.4, 0.8, 1.2))
merged_seurat$condition <- factor(merged_seurat$condition, levels=c('NI', 'IVA', 'RV'))
Idents(merged_seurat) <- 'condition' 

# plots
DimPlot(merged_seurat, reduction='rna.umap', group.by='condition')
ggsave(filename='UMAP_mergedNI_RV_IVA.pdf', height=5, width=6)

DimPlot(merged_seurat, reduction='rna.umap', group.by='RNA_snn_res.0.4')
ggsave(filename='UMAP_mergedNI_RV_IVA_byclusters.pdf', height=5, width=6)

DimPlot(merged_seurat, reduction='rna.umap', group.by=c('condition','batch'))
ggsave(filename='UMAP_mergedNI_RV_IVA_byconditionandbatch.pdf', height=4, width=10)

DimPlot(merged_seurat, reduction='rna.umap', group.by='batch', split.by='condition')
ggsave(filename='UMAP_mergedNI_RV_IVA_byconditionandbatch.pdf', height=4, width=10)

VlnPlot(merged_seurat, features=c('nFeature_RNA','nCount_RNA'), ncol=2, alpha=0.01)
ggsave(filename='VlnPlot_mergedNI_RV_IVA_nFeatureRNA.nCountRNA.pdf', height=4, width=6)
ggsave(filename='VlnPlot_mergedNI_RV_IVA_nFeatureRNA.nCountRNA.png', height=4, width=6)

# save integrated Seurat object
saveRDS(merged_seurat, file='NI_RV_IVA_integrated.rds')