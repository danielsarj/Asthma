library(tidyverse)
library(Seurat)
library(data.table)
library(patchwork)
"%&%" <- function(a,b) paste(a,b, sep = "")
options(future.globals.maxSize=120*1024^3)
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
        rename(IDs=as.character(0)) %>% filter(!grepl('\\+', IDs))
      cells_to_keep$IDs <- gsub('\\_.*', '', cells_to_keep$IDs)
    
      # load cellranger alignment results for batch 1
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix') %>% 
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep$BARCODE)
      
      # add metadata to seurat object
      counts[['percent.mt']] <- PercentageFeatureSet(counts, pattern='^MT-')
      counts[['batch']] <- as.factor(batches[b])
      counts[['condition']] <- as.factor(conditions[c])
      mdata <- counts@meta.data %>% rownames_to_column('BARCODE') %>% 
        inner_join(cells_to_keep, by='BARCODE') %>% column_to_rownames('BARCODE')
      counts@meta.data <- mdata
      rm(mdata)

      # perform initial QC and save object in a list
      counts <- counts %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<10) 
      seurat_objs[[b]] <- counts
      
      } else if (b==length(batches)){
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0)) %>% filter(!grepl('\\+', IDs))
      cells_to_keep$IDs <- gsub('\\_.*', '', cells_to_keep$IDs)
      
      # load cellranger alignment results for the last batch
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix') %>% 
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep$BARCODE)
      
      # add metadata to seurat object
      counts[['percent.mt']] <- PercentageFeatureSet(counts, pattern='^MT-')
      counts[['batch']] <- as.factor(batches[b])
      counts[['condition']] <- as.factor(conditions[c])
      mdata <- counts@meta.data %>% rownames_to_column('BARCODE') %>% 
        inner_join(cells_to_keep, by='BARCODE') %>% column_to_rownames('BARCODE')
      counts@meta.data <- mdata
      rm(mdata)
      
      # perform initial QC and save object in a list
      counts <- counts %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<10) 
      seurat_objs[[b]] <- counts
      
      # save list of raw batches
      saveRDS(seurat_objs, file='../scRNAanalysis/'%&%conditions[c]%&%'_raw.objects.batches.rds')

      } else {
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0)) %>% filter(!grepl('\\+', IDs))
      cells_to_keep$IDs <- gsub('\\_.*', '', cells_to_keep$IDs)
      
      # load cellranger alignment results for corresponding batch
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix') %>% 
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep$BARCODE)
      
      # add metadata to seurat object
      counts[['percent.mt']] <- PercentageFeatureSet(counts, pattern='^MT-')
      counts[['batch']] <- as.factor(batches[b])
      counts[['condition']] <- as.factor(conditions[c])
      mdata <- counts@meta.data %>% rownames_to_column('BARCODE') %>% 
        inner_join(cells_to_keep, by='BARCODE') %>% column_to_rownames('BARCODE')
      counts@meta.data <- mdata
      rm(mdata)
      
      # perform initial QC and save object in a list
      counts <- counts %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<10) 
      seurat_objs[[b]] <- counts
    }
  }
}

umaps[[1]] | umaps[[2]] | umaps[[3]]
ggsave(filename='../scRNAanalysis/UMAP_QCintegratedbatches.pdf', height=6, width=18)