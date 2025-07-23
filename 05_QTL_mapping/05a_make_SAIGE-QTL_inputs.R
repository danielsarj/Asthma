library(Seurat)
library(tidyverse)
library(data.table)
library(PCAForQTL)
library(irlba)
library(Matrix)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

# load Seurat object
obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.w_celltype.rds') 

# get ID metadata
sample_m <- fread('../sample_metadata.txt') %>% select(ID, age, gender)
sample_m$gender <- as.factor(sample_m$gender) %>% as.numeric()

# get genetic PCs
genetic_pcs <- fread('PCAIR.eigenvec') %>% select(sample_id, V1, V2, V3, V4)
genetic_pcs$sample_id <- gsub('SEA3', 'SEA-3', genetic_pcs$sample_id)

# merge all metadatas 
seurat_md <- obj@meta.data %>% rownames_to_column('cell_ID')
mdata <- inner_join(seurat_md, genetic_pcs, by=c('IDs'='sample_id')) %>% 
  inner_join(sample_m, by=c('IDs'='ID')) %>% 
  select(-c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, batch, RNA_snn_res.0.1, 
            RNA_snn_res.0.4, RNA_snn_res.0.8, RNA_snn_res.1.2, seurat_clusters))
rownames(mdata) <- mdata$cell_ID
obj@meta.data <- mdata
rm(genetic_pcs, mdata, sample_m, seurat_md)

# remove individuals for which there are NAs in covariates
complete_cells <- complete.cases(obj@meta.data)
obj <- obj %>% subset(cells = colnames(obj)[complete_cells])
rm(complete_cells)

# get gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt') %>% filter(gene_biotype=='protein_coding',
                hgnc_symbol!='', !grepl('^MT-', hgnc_symbol))
annotations$chromosome_name <- as.numeric(annotations$chromosome_name)
annotations <- annotations %>% drop_na()

# work with one condition-cell type at a time
for (cond in c('NI', 'IVA', 'RV')){
  print(cond)
  for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
    print(ctype)
  
    # subset seurat object
    subset_obj <- subset(obj, subset = celltype == ctype & condition == cond)
    
    # extract count matrix
    count_mat <- subset_obj@assays$RNA$counts
  
    # keep only protein-coding genes with variance > 0
    keep_genes <- rownames(count_mat) %in% annotations$hgnc_symbol
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
  
    # join metadata, exp_pcs, and count_mat
    full_df <- inner_join(subset_obj@meta.data, exp_pcs, by=c('cell_ID')) %>% 
      inner_join(count_mat, by=c('cell_ID')) %>% arrange(IDs) %>% select(-c(condition, celltype))
  
    # make chr-gene df 
    sub_anno <- annotations %>% dplyr::select(hgnc_symbol, chromosome_name) %>% 
      filter(hgnc_symbol %in% colnames(full_df))
    
    # save files
    fwrite(full_df, 'SAIGE_results/'%&%cond%&%'_'%&%ctype%&%'_counts.w.covs.txt', sep='\t', col.names=TRUE)
    fwrite(sub_anno, 'SAIGE_results/'%&%cond%&%'_'%&%ctype%&%'_gene_list.txt', sep='\t', col.names=FALSE)
  
    rm(count_mat, exp_pcs, sub_anno, full_df)
  }
}