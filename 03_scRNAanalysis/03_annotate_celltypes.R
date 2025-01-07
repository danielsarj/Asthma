library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

conditions <- c('NI','RV','IVA')
files <- list.files(pattern='\\_allbatches.scaled.clustered.rds$')
umaps <- list()

for (i in 1:length(files)){
  print(c(conditions[i]))
  
  objs <- readRDS(conditions[i] %&%'_allbatches.scaled.clustered.rds')
  DefaultAssay(objs) <- 'RNA'
  
  # annotate celltypes using a PBMC reference and filter cells by mapping.score > 50%
  print(objs)
  objs <- RunAzimuth(objs, reference='pbmcref', verbose=FALSE)
  objs <- objs %>% subset(subset=mapping.score>0.5)
  print(objs)
  saveRDS(objs, file=conditions[i] %&%'_allbatches.scaled.clustered.celltypeannotation.rds')
  
  # save UMAP viz
  umaps[[i]] <- DimPlot(objs, group.by='predicted.celltype.l1', label=TRUE, label.size=5) + 
    ggtitle(conditions[i])
}

umaps[[1]] | umaps[[2]] | umaps[[3]]
ggsave(filename='UMAP_celltypes.allconditions.pdf', height=6, width=18)