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
obj$log_total_counts <- log(obj$nCount_RNA)

# keep only NI 
obj <- obj %>% subset(subset = SOC_infection_status=='NI')

# ID metadata
mdata <- fread('HALEYs/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt') %>% 
  select(indiv_ID, age_Scale, YRI_Scale) %>% unique() 

# join ID metadata to Seurat's metadata
s_mdata <- obj@meta.data %>% rownames_to_column('cell_ID') %>% 
  inner_join(mdata, by=c('SOC_indiv_ID'='indiv_ID'), relationship='many-to-many')
rownames(s_mdata) <- rownames(obj@meta.data)
obj@meta.data <- s_mdata
rm(s_mdata, mdata)

# remove individuals for which there are NAs in covariates
complete_cells <- complete.cases(obj@meta.data)
obj <- obj %>% subset(cells = colnames(obj)[complete_cells])
rm(complete_cells)

# load PLINK .fam file to create 10 permutations of genotypes
plink_fam <- fread('Saige/step2/inputs/Haley_filtered_genotypes.fam')
for (permutation in seq(1:10)){
  new_fam <- plink_fam
  ids <- new_fam[, 1:2] 
  shuffled_ids <- ids[sample(nrow(ids)), ]
  new_fam[, 1:2] <- shuffled_ids
  
  fwrite(new_fam, 'Saige/step2/inputs/Haley_filtered_genotypes_'%&%permutation%&%'.fam')
  file.copy('Saige/step2/inputs/Haley_filtered_genotypes.bed', 'Saige/step2/inputs/Haley_filtered_genotypes_'%&%permutation%&%'.bed', 
            overwrite=TRUE)
  file.copy('Saige/step2/inputs/Haley_filtered_genotypes.bim', 'Saige/step2/inputs/Haley_filtered_genotypes_'%&%permutation%&%'.bim', 
            overwrite=TRUE)
}

# get gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt') %>% filter(gene_biotype=='protein_coding',
                     hgnc_symbol!='', !grepl('^MT-', hgnc_symbol))
annotations$chromosome_name <- as.numeric(annotations$chromosome_name)
annotations <- annotations %>% drop_na()

# work with one cell type at a time
for (ctype in c('B','CD4_T','CD8_T','monocytes','NK')){
  print(ctype)
  
  # subset seurat object
  subset_obj <- obj %>% subset(subset = celltype == ctype)
  
  # extract count 
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
  exp_pcs <- exp_pcs$x[, 1:10] %>% as.matrix()
  colnames(exp_pcs) <- paste0("PC", 1:10)
  exp_pcs <- as.data.frame(exp_pcs) %>% mutate(cell_ID = rownames(count_mat))
    
  # adjust count_mat to append metadata
  count_mat <- count_mat %>% as.data.frame() %>% rownames_to_column('cell_ID')
    
  # join metadata, exp_pcs_df, and count_mat
  full_df <- inner_join(subset_obj@meta.data, exp_pcs, by=c('cell_ID')) %>% 
    inner_join(count_mat, by=c('cell_ID')) %>% arrange(SOC_indiv_ID) %>% 
    select(-c(orig.ident, nCount_RNA, nFeature_RNA, batchID, SOC_status,
                SOC_infection_status, SOC_genetic_ancestry, CEU, YRI, nCount_SCT, nFeature_SCT,
                integrated_snn_res.0.5, cluster_IDs, celltype, sample_condition))
    
  # save count file
  fwrite(full_df, 'Saige/step1/inputs/'%&%ctype%&%'_NI_counts.w.covs_upto10PCs.txt', sep='\t', col.names=TRUE)

  # make chr-gene df 
  sub_anno <- annotations %>% dplyr::select(hgnc_symbol, chromosome_name) %>% 
    filter(hgnc_symbol %in% colnames(full_df))
  fwrite(sub_anno, 'Saige/step1/inputs/'%&%ctype%&%'_NI_gene_list.txt', sep='\t', col.names=FALSE)
  
  rm(full_df, sub_anno, exp_pcs, count_mat, subset_obj)
}
