library(tidyverse)
library(Seurat)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('NI','RV','IVA')
files <- list.files(pattern='\\_integrated.rds$')
objs <- list()

for (i in 1:length(files)){
  a <- readRDS('../scRNAanalysis/'%&% conditions[i] %&%'_integrated.rds')
  a[['RNA']] <- JoinLayers(a[['RNA']])
  a <- ScaleData(a)
  a <- RunPCA(a)
  a <- FindNeighbors(a, dims = 1:30) %>%
    FindClusters()
  a <- RunUMAP(a, dims = 1:30)
  objs[[i]] <- a
}

DimPlot(objs[[1]], reduction='umap', group.by='orig.ident') | 
  DimPlot(objs[[2]], reduction='umap', group.by='orig.ident') | 
  DimPlot(objs[[3]], reduction='umap', group.by='orig.ident')

ggsave(filename='UMAP_integratedbatches.pdf', height=5, width=14)
