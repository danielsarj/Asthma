library(Seurat)
library(SeuratData)
library(limma)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/Volumes/daniel/asthma_project/QTLmapping')
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

# define rank-inverse normalization function
rin_transform <- function(x) {
  r <- rank(x, ties.method='average')
  n <- length(x)
  qnorm((r-0.5)/n) 
}

# load sample metadata
sample_m <- fread('../DEanalysis/sample_metadata.txt')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# load genotype PCs
geno_pcs <- fread('PCAIR.eigenvec') %>% select(sample_id, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)
geno_pcs$sample_id <- gsub('SEA3', 'SEA-3', geno_pcs$sample_id)

for (i in 1:length(conditions)){
  # load pseudobulk object
  print(conditions[i])
  
  if (conditions[i]=='IVA'){
    obj <- readRDS('../DEanalysis/NI_IVA_pseudobulks.rds')
  } else {
    obj <- readRDS('../DEanalysis/NI_RV_pseudobulks.rds')
  }

  for (j in 1:length(celltypes)){
    print(celltypes[j])
    
    # extract metadata
    meta_df <- obj@meta.data
    filtered_meta <- meta_df[meta_df$predicted.celltype.l1 == celltypes[j] & meta_df$condition == conditions[i], ]
    
    # subset bulk object
    matching_cells <- rownames(filtered_meta)
    tmp <- subset(obj, cells=matching_cells)
    
    # edit metadata
    filtered_meta$IDs <- gsub('SEA3', 'SEA-3', filtered_meta$IDs)
    filtered_meta <- left_join(filtered_meta, sample_m, by=c('IDs'='ID')) %>% 
      left_join(geno_pcs, by=c('IDs'='sample_id'))
    filtered_meta$gender <- as.factor(filtered_meta$gender)
    
    # extract count 
    count_df <- tmp@assays$RNA$counts %>% as.data.frame()
    count_df <- count_df[rownames(count_df) %in% annotations,]
    zero_var_genes <- apply(count_df, 1, var) == 0
    low_exp_genes <- rowMeans(count_df) > 0.4
    count_df <- count_df[!zero_var_genes & low_exp_genes, ]
    exp_pcs <- prcomp(count_df, scale.=TRUE)$rotation %>% as.data.frame() %>% rownames_to_column('IDs') %>%
      select(IDs, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
    exp_pcs$IDs <- gsub('SEA3', 'SEA-3', exp_pcs$IDs)
    exp_pcs$IDs <- gsub('_'%&%celltypes[j]%&%'_'%&%conditions[i], '', exp_pcs$IDs)
    
    # finalize metadata dataframe
    filtered_meta <- left_join(filtered_meta, exp_pcs, by=c('IDs'))
    
    # get linear regression residuals after adjusting for age, sex, 
    # first 10 genotype PCs, first 10 expression PCs, and obtain rin residuals
    log_counts <- log2(count_df+1)
    design <- model.matrix(~age+gender+V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=filtered_meta)
    fit <- lmFit(log_counts, design)
    residual_matrix <- residuals.MArrayLM(fit, log_counts)
    rin_residuals <- apply(residual_matrix, 1, rin_transform)
    
    # edit matrix and save results
    rin_residuals <- rin_residuals %>% t() %>% 
      as.data.frame() %>% rownames_to_column(var='GENES')
    tmp_colnames <- colnames(rin_residuals)
    tmp_colnames <- gsub('SEA3', 'SEA-3', tmp_colnames)
    tmp_colnames <- gsub('_'%&%celltypes[j]%&%'_'%&%conditions[i], '', tmp_colnames)
    colnames(rin_residuals) <- tmp_colnames
    fwrite(rin_residuals, conditions[i]%&%'_'%&%celltypes[j]%&%'_rinResiduals.txt', col.names=T, sep='\t')
    
  }
}