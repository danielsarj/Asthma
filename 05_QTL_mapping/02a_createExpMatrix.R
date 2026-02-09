library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(tidyverse)
library(data.table)
library(PCAForQTL)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

# function to regress out k number of pcs from expression data frame
pca_rm <- function(input_data, pc_set) {
  pca <- prcomp(t(input_data))
  new <- input_data
  new <- apply(new, 1, FUN = function(x){return(lm(as.numeric(x) ~ -1 + pca$x[, as.numeric(pc_set)])$resid)})
  new <- t(new)
  colnames(new) <- colnames(input_data)
  rownames(new) <- rownames(input_data)
  return(new)
}

# function for quantile normalization
quantile_norm <- function(df){
  quantile_expression <- matrix(, nrow = nrow(df), ncol = ncol(df))
  for (j in 1:nrow(quantile_expression)){
    exp <- df[j,]
    exp_QN <- qqnorm(exp, plot = FALSE)
    exp_QN <- exp_QN$x
    quantile_expression[j,] <- exp_QN
  }
  rownames(quantile_expression) <- rownames(df)
  colnames(quantile_expression) <- colnames(df)
  return(quantile_expression)
}

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# load pseudobulk object
obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks_new.rds')

# remove batch 4
obj <- subset(obj, subset= batch!='B4')

for (i in 1:length(conditions)){
  print(conditions[i])
  for (j in 1:length(celltypes)){
    print(celltypes[j])
    
    # extract metadata
    meta_df <- obj@meta.data
    filtered_meta <- meta_df[meta_df$celltype == celltypes[j] & meta_df$condition == conditions[i], ]
    
    # subset bulk object
    matching_cells <- rownames(filtered_meta)
    tmp <- subset(obj, cells=matching_cells)
    
    # edit metadata
    filtered_meta <- left_join(filtered_meta, sample_m, by=c('IDs'='ID'))
    filtered_meta$gender <- as.factor(filtered_meta$gender)
    
    # extract count 
    count_df <- tmp@assays$RNA$counts %>% as.data.frame()
    count_df <- count_df[rownames(count_df) %in% annotations,]
    zero_var_genes <- apply(count_df, 1, var) == 0
    count_df <- count_df[!zero_var_genes, ]
    
    # normalize and remove genes with low expression (min. of 2 CPM across 5 samples)
    dge <- DGEList(counts=count_df)
    dge <- calcNormFactors(dge)
    cpm_values <- cpm(dge)
    threshold <- 2
    min_samples <- 5
    keep <- rowSums(cpm_values > threshold) >= min_samples
    dge <- dge[keep, ]

    # adjust for batch, age, gender, number of cells, average MT content, and k expression PCs
    design <- model.matrix(~batch+age+gender+n+avg_mt, data=filtered_meta)
    expression <- voom(dge, design, plot=FALSE)
    fit <- lmFit(expression, design)
    expression <- residuals.MArrayLM(fit, expression)
      
    # find best number of PCs according to the elbow method
    pc_set <- prcomp(expression, center=TRUE, scale.=TRUE)
    resultRunElbow <- runElbow(prcompResult=pc_set)
    print('number of PCs per elbow method: '%&% resultRunElbow)
    expression <- quantile_norm(expression)
    
    for (k in seq(from=0, to=20)){
      if (k>0){
        adj_expression <- pca_rm(expression, c(1:k))
      } else {
        adj_expression <- expression
      }

      # edit matrix and save results
      adj_expression <- adj_expression %>% as.data.frame() %>% rownames_to_column(var='GENES')
      tmp_colnames <- colnames(adj_expression)
      tmp_colnames <- gsub('_'%&%conditions[i]%&%'_'%&%celltypes[j], '', tmp_colnames)
      colnames(adj_expression) <- tmp_colnames
      fwrite(adj_expression, conditions[i]%&%'_'%&%celltypes[j]%&%'_'%&%k%&%'PCs_new.txt', col.names=T, sep='\t')
    }
  }
}
