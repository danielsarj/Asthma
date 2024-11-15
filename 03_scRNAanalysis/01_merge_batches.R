library(tidyverse)
library(Seurat)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/alignment')

batches <- c('B1','B2','B3','B4')
conditions <- c('NI','RV','IVA')
seurat_objs <- list()

for (c in 1:length(conditions)){
  for (b in 1:length(batches)){
    print(c(batches[b], conditions[c]))
    if (b==1){
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0))
      cells_to_keep <- cells_to_keep %>% filter(!grepl('\\+', IDs)) %>% 
        select(BARCODE) %>% pull()
    
      # load cellranger alignment results for batch 1
      counts <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix')  %>%
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep) %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c])
    } else if (b==4){
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0))
      cells_to_keep <- cells_to_keep %>% filter(!grepl('\\+', IDs)) %>% 
        select(BARCODE) %>% pull()
      
      # load cellranger alignment results for batch 4 and merge with previous batches
      tmp <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix')  %>%
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep) %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c])
      counts <- merge(counts, tmp)
      
      # save merged seurat object in a list
      seurat_objs[[c]] <- counts
    } else {
      # get barcodes with donors assigned
      cells_to_keep <- fread(batches[b]%&%'-'%&%conditions[c]%&%'_GRCh38/demuxalot/assignments_refined.tsv.gz', header=T) %>%
        rename(IDs=as.character(0))
      cells_to_keep <- cells_to_keep %>% filter(!grepl('\\+', IDs)) %>% 
        select(BARCODE) %>% pull()
      
      # load cellranger alignment results and merge with previous batch(es)
      tmp <- Read10X(batches[b]%&%'-'%&%conditions[c]%&%'/outs/filtered_feature_bc_matrix')  %>%
        CreateSeuratObject(project=batches[b]%&%'_'%&%conditions[c]) %>% 
        subset(cells=cells_to_keep) %>% RenameCells(add.cell.id=batches[b]%&%'_'%&%conditions[c])
      counts <- merge(counts, tmp)
    }
  }
}

# save list of merged batches 
saveRDS(seurat_objs, file='../scRNAanalysis/seurat_objs_raw.rds')