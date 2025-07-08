library(Seurat)
library(tidyverse)
library(data.table)
library(PCAForQTL)
library(irlba)
library(Matrix)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

# read files, update metadata
mdata <- readRDS('HALEYs/pseudobulks/all_metadata.rds')
obj <- readRDS('HALEYs/pseudobulks/all_raw_counts.rds') %>% CreateSeuratObject()
obj@meta.data <- mdata

# keep only NI 
obj <- obj %>% subset(subset = SOC_infection_status=='NI')

# ID metadata
mdata <- fread('HALEYs/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt') %>% 
  select(indiv_ID, age_Scale, YRI_Scale) %>% drop_na() %>% unique() 

# join ID metadata to Seurat's metadat
s_mdata <- obj@meta.data %>% full_join(mdata, by=c('SOC_indiv_ID'='indiv_ID'), relationship='many-to-many')
rownames(s_mdata) <- rownames(obj@meta.data)
obj@meta.data <- s_mdata
rm(s_mdata, mdata)

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# work with one cell type at a time
for (ctype in c('B','CD4_T','CD8_T','monocytes','NK')){
  print(ctype)
  
  # subset seurat object
  subset_obj <- obj %>% subset(subset = celltype == ctype)
  
  # extract count 
  count_mat <- subset_obj@assays$RNA$counts
  
  # keep only protein-coding genes with variance > 0
  keep_genes <- rownames(count_mat) %in% annotations
  count_mat <- count_mat[keep_genes, ]
  gene_vars <- apply(count_mat, 1, function(x) var(as.numeric(x)))
  count_mat <- count_mat[gene_vars > 0, ]
  rm(gene_vars, keep_genes)
  
  # convert to dense matrix for PC computing
  count_mat <- count_mat %>% as.matrix() %>% t()
  
  # compute expression PCs
  exp_pcs <- prcomp_irlba(count_mat, n=20, scale.=TRUE, center=TRUE)
  
  # find best K 
  class(exp_pcs) <- 'prcomp'
  K_elbow <- runElbow(prcompResult=exp_pcs)
  exp_pcs <- exp_pcs$x[,c(1:K_elbow)] %>% as.data.frame() %>% mutate(cell_ID=rownames(count_mat))

  # adjust count_mat to append metadata
  count_mat <- count_mat %>% as.data.frame() %>% rownames_to_column('cell_ID')
  
  # subset metadata
  mdata <- subset_obj@meta.data %>% as.data.frame() %>% rownames_to_column('cell_ID') %>%
    select(cell_ID, SOC_indiv_ID, age_Scale, YRI_Scale)
  
  # join metadata, exp_pcs, and count_mat
  full_df <- inner_join(mdata, exp_pcs, by=c('cell_ID')) %>% inner_join(count_mat, by=c('cell_ID')) %>% 
    select(-cell_ID) %>% arrange(SOC_indiv_ID)
  rm(count_mat, exp_pcs, mdata)
  
  fwrite(full_df, 'Saige/step1/inputs/'%&%ctype%&%'_NI_counts.w.covs.txt', sep='\t', col.names=TRUE)
}

length(full_df$NOC2L)
length(full_df$age_Scale)
length(full_df$YRI_Scale)
length(full_df$PC1)

