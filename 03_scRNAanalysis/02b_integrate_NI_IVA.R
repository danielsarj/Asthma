library(tidyverse)
library(Seurat)
library(patchwork)
library(glmGamPoi)
library(future)
plan('multicore', workers=4)
"%&%" <- function(a,b) paste(a,b, sep = "")
options(future.globals.maxSize=200*1024^3)
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('NI','IVA')
seurat_objs <- list()

# read each file
for (c in 1:length(conditions)){
  seurat_objs <- append(seurat_objs, readRDS(conditions[c]%&%'_raw.objects.batches.rds'))
}

# SCTransform each batch/condition
for (i in 1:length(seurat_objs)){
  seurat_objs[[i]] <- seurat_objs[[i]] %>% SCTransform(vars.to.regress=c('percent.mt'), 
                                                       conserve.memory=TRUE, method='glmGamPoi')
}

# select features for integration
features <- SelectIntegrationFeatures(object.list=seurat_objs, nfeatures=3000)

# prepare for SCT integration
seurat_objs <- PrepSCTIntegration(object.list=seurat_objs, anchor.features=features)

# find anchors
anchors <- FindIntegrationAnchors(object.list=seurat_objs, normalization.method='SCT', 
                                  anchor.features=features)

# integrate data
integrated_obj <- IntegrateData(anchorset=anchors, normalization.method='SCT')

# PCA, clustering, viz
integrated_obj <- integrated_obj %>% RunPCA(npcs=30) %>% FindNeighbors(dims=1:30) %>% 
  FindClusters() %>% RunUMAP(dims=1:30)
DimPlot(integrated_obj, reduction='umap', group.by=c('batch', 'condition'))
ggsave(filename='UMAP_QCintegratedNI_IVA.pdf', height=6, width=13)

DimPlot(integrated_obj, reduction='umap', group.by='batch', split.by='condition')
ggsave(filename='UMAP_QCintegratedNI_IVA_splitcondition.pdf', height=6, width=18)

# save list of integrated batches/ditions
saveRDS(integrated_obj, file='NI_IVA_integrated.rds')