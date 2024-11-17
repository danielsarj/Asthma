library(tidyverse)
library(Seurat)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/alignment')

batches <- c('B1','B2','B3','B4')
conditions <- c('NI','RV','IVA')

for (c in 1:length(conditions)){
  for (b in 1:length(batches)){
    print(c(batches[b], conditions[c]))
    if (b==1){
      seurat_objs <- list()
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0))
      cells_to_keep <- cells_to_keep %>% filter(!grepl('\\+', IDs)) %>% 
        select(BARCODE) %>% pull()
    
      # load cellranger alignment results for batch 1
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix') %>% 
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep) %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c])
      
      # perform QC and save object in a list
      counts[['percent.mt']] <- PercentageFeatureSet(counts, pattern='^MT-')
      counts <- counts %>% subset(subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5) %>%
        NormalizeData() %>% FindVariableFeatures()
      seurat_objs[[b]] <- counts
      
      } else if (b==4){
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0))
      cells_to_keep <- cells_to_keep %>% filter(!grepl('\\+', IDs)) %>% 
        select(BARCODE) %>% pull()
      
      # load cellranger alignment results for batch 4
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix') %>% 
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep) %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c])
      
      # perform QC and save object in a list
      counts[['percent.mt']] <- PercentageFeatureSet(counts, pattern='^MT-')
      counts <- counts %>% subset(subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5) %>%
        NormalizeData() %>% FindVariableFeatures()
      seurat_objs[[b]] <- counts
      
      # save list of QC'd batches
      saveRDS(seurat_objs, file='../scRNAanalysis/'%&%conditions[c]%&%'_QCd.batches.rds')
      
      # find integration anchors
      anchors <- FindIntegrationAnchors(object.list=seurat_objs, dims=1:30)
      
      # integrate the datasets & save
      combined <- IntegrateData(anchorset=anchors, dims=1:30)
      
      # save list of QC'd batches
      saveRDS(combined, file='../scRNAanalysis/'%&%conditions[c]%&%'_integrated.rds')
      
    } else {
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0))
      cells_to_keep <- cells_to_keep %>% filter(!grepl('\\+', IDs)) %>% 
        select(BARCODE) %>% pull()
      
      # load cellranger alignment results other batches
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix') %>% 
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep) %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c])
      
      # perform QC and save object in a list
      counts[['percent.mt']] <- PercentageFeatureSet(counts, pattern='^MT-')
      counts <- counts %>% subset(subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5) %>%
        NormalizeData() %>% FindVariableFeatures()
      seurat_objs[[b]] <- counts
    }
  }
}