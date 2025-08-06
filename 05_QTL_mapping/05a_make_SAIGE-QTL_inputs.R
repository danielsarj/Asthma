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

# define chromosome sizes
chrom_sizes <- data.frame(
  chromosome = seq(1:22),
  size = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468
  )
)

# create regions file
cis_regions <- annotations %>% select(hgnc_symbol, chromosome_name, start_position, end_position) %>%
  mutate(start_position=start_position-(5e5), end_position=end_position+(5e5))

# fix out-of-bonds cis regions 
cis_regions <- cis_regions %>% left_join(chrom_sizes, by=c('chromosome_name'='chromosome')) %>%
  mutate(start_position=pmax(start_position, 1), end_position=pmin(end_position, size)) %>% select(-size)
fwrite(cis_regions, 'SAIGE_results/cis_regions.txt', sep=' ', col.names=F)
rm(cis_regions, chrom_sizes)

# load PLINK .fam file to create 10 permutations of genotypes
plink_fam <- fread('SAIGE_results/filtered_genotypes.fam')
for (permutation in seq(1:10)){
  new_fam <- plink_fam
  ids <- new_fam[, 1:2] 
  shuffled_ids <- ids[sample(nrow(ids)), ]
  new_fam[, 1:2] <- shuffled_ids
  
  fwrite(new_fam, 'SAIGE_results/filtered_genotypes_'%&%permutation%&%'.fam')
  file.copy('SAIGE_results/filtered_genotypes.bed', 'SAIGE_results/filtered_genotypes_'%&%permutation%&%'.bed', 
            overwrite=TRUE)
  file.copy('SAIGE_results/filtered_genotypes.bim', 'SAIGE_results/filtered_genotypes_'%&%permutation%&%'.bim', 
            overwrite=TRUE)
}

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
    exp_pcs <- prcomp_irlba(count_mat, n=10, scale.=TRUE, center=TRUE)
    class(exp_pcs) <- 'prcomp'
    
    # try different PCs
    for (pcs in seq(1:10)){
      print(pcs)
      exp_pcs_df <- exp_pcs$x[,c(1:pcs)] %>% as.data.frame() %>% mutate(cell_ID=rownames(count_mat))

      # adjust count_mat to append metadata
      count_mat_df <- count_mat %>% as.data.frame() %>% rownames_to_column('cell_ID')
  
      # join metadata, exp_pcs_df, and count_mat
      full_df <- inner_join(subset_obj@meta.data, exp_pcs_df, by=c('cell_ID')) %>% 
        inner_join(count_mat_df, by=c('cell_ID')) %>% arrange(IDs) %>% select(-c(condition, celltype))
  
      # save files
      fwrite(full_df, 'SAIGE_results/'%&%cond%&%'_'%&%ctype%&%'_counts.w.covs_'%&%pcs%&%'.txt', sep='\t', col.names=TRUE)
    }
    rm(exp_pcs, exp_pcs_df, count_mat, count_mat_df, subset_obj)
    # make chr-gene df 
    sub_anno <- annotations %>% dplyr::select(hgnc_symbol, chromosome_name) %>% 
      filter(hgnc_symbol %in% colnames(full_df))
    fwrite(sub_anno, 'SAIGE_results/'%&%cond%&%'_'%&%ctype%&%'_gene_list.txt', sep='\t', col.names=FALSE)
    rm(full_df, sub_anno)
  }
}
