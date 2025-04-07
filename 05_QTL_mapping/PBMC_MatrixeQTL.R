library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(tidyverse)
library(data.table)
library(PCAForQTL)
library(janitor)
library(MatrixEQTL)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')
conditions <- c('NI', 'RV', 'IVA')

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

# load sample metadata
sample_m <- fread('../DEanalysis/sample_metadata.txt')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# #for (i in 1:length(conditions)){
#   # load pseudobulk object
#   print(conditions[i])
#   
#   if (conditions[i]=='IVA'){
#     obj <- readRDS('../DEanalysis/NI_IVA_pseudobulks.rds')
#   } else {
#     obj <- readRDS('../DEanalysis/NI_RV_pseudobulks.rds')
#   }
#   
#   # get metadata before summing cell types
#   pre_metadata <- obj@meta.data %>% filter(condition==conditions[i]) %>%
#     select(IDs, predicted.celltype.l1, n_cells)
#   
#   # get pseudobulk for all cell types
#   bulk_objs <- AggregateExpression(obj, group.by=c('IDs','condition'), 
#                                    slot='counts', assays='RNA', return.seurat=T)
#   
#   # extract metadata
#   meta_df <- bulk_objs@meta.data
#   filtered_meta <- meta_df[meta_df$condition == conditions[i], ]
#   
#   # subset bulk object
#   matching_cells <- rownames(filtered_meta)
#   tmp <- subset(bulk_objs, cells=matching_cells)
#   
#   # edit metadata
#   filtered_meta$IDs <- gsub('SEA3', 'SEA-3', filtered_meta$IDs)
#   filtered_meta <- left_join(filtered_meta, sample_m, by=c('IDs'='ID'))
#   filtered_meta$gender <- as.factor(filtered_meta$gender)
#   
#   # add quantity of cell types to PBMC metadata
#   pre_metadata$predicted.celltype.l1 <- gsub('-', '_', pre_metadata$predicted.celltype.l1)
#   pre_metadata$IDs <- gsub('SEA3', 'SEA-3', pre_metadata$IDs)
#   pre_metadata <- pre_metadata %>% pivot_wider(names_from=predicted.celltype.l1, 
#                                                values_from=n_cells, values_fill=0) %>%
#     select(-c('DC', 'other', 'other_T'))
#   filtered_meta <- full_join(filtered_meta, pre_metadata, by=c('IDs'))
#   
#   # extract count 
#   count_df <- tmp@assays$RNA$counts %>% as.data.frame()
#   count_df <- count_df[rownames(count_df) %in% annotations,]
#   zero_var_genes <- apply(count_df, 1, var) == 0
#   count_df <- count_df[!zero_var_genes, ]
#   exp_pcs <- prcomp(count_df, scale.=TRUE, center=TRUE)
# 
#   # find best K 
#   K_elbow <- runElbow(prcompResult=exp_pcs)
#   pc_set <- c(1:K_elbow)
#   
#   # normalize and remove genes with low expression (min. of 2 CPM across 3 samples)
#   dge <- DGEList(counts=count_df)
#   dge <- calcNormFactors(dge)
#   cpm_values <- cpm(dge)
#   threshold <- 2
#   min_samples <- 3
#   keep <- rowSums(cpm_values > threshold) >= min_samples
#   dge <- dge[keep, ]
#   dge <- DGEList(counts=count_df)
#   dge <- calcNormFactors(dge)
#   
#   # adjust for age, gender, cell type composition, and expression PCs
#   design <- model.matrix(~age+gender+B+CD4_T+CD8_T+Mono+NK, data=filtered_meta)
#   expression <- voom(dge, design, plot=FALSE)$E
#   expression <- pca_rm(expression, pc_set)
#   
#   # edit matrix and save results
#   expression <- expression %>% as.data.frame() %>% rownames_to_column(var='GENES')
#   tmp_colnames <- colnames(expression)
#   tmp_colnames <- gsub('SEA3', 'SEA-3', tmp_colnames)
#   tmp_colnames <- gsub('_'%&%conditions[i], '', tmp_colnames)
#   colnames(expression) <- tmp_colnames
#   fwrite(expression, conditions[i]%&%'_PBMC_elbowPCs.txt', col.names=T, sep='\t')  
# }
# #rm(bulk_objs, count_df, design, dge, exp_pcs, expression, filtered_meta, 
# #   meta_df, obj, pre_metadata, sample_m, tmp, annotations, i, K_elbow,
# #   low_exp_genes, matching_cells, pc_set, tmp_colnames, zero_var_genes, pca_rm)
  
### MATRIXEQTL PART

# load snp and gene location files
snp_local <- fread('../genotypes/imputed_vcfs/snp_location.txt')
gene_local <- fread('gene_location.txt')

n_permutations <- 10
for (i in 1:length(conditions)){
  # load pseudobulk object
  print(conditions[i])
  
  # load expression matrix and dosage file
  # make sure the columns are in the same order
  exp_matrix <- fread(conditions[i]%&%'_PBMC_elbowPCs.txt')
  dos_matrix <- fread('../genotypes/imputed_vcfs/imputed_dosage.txt')
  common_cols <- intersect(names(exp_matrix), names(dos_matrix))  
  dos_matrix <- cbind(dos_matrix[,1], dos_matrix[, ..common_cols])
  setcolorder(dos_matrix, c(names(dos_matrix)[1], setdiff(common_cols, names(dos_matrix)[1])))
  
  # load genotype PCs and make sure columns are in the correct order
  geno_pcs <- fread('PCAIR.eigenvec') %>% select(sample_id, V1, V2, V3, V4)
  geno_pcs$sample_id <- gsub('SEA3', 'SEA-3', geno_pcs$sample_id)
  geno_pcs <- geno_pcs %>% t() %>% as.data.frame() %>% row_to_names(row_number=1) %>%
    rownames_to_column() %>% setDT()
  common_cols <- intersect(names(exp_matrix), names(geno_pcs))  
  geno_pcs <- cbind(geno_pcs[,1], geno_pcs[, ..common_cols])
  setcolorder(geno_pcs, c(names(geno_pcs)[1], setdiff(common_cols, names(geno_pcs)[1])))
  geno_pcs[, (2:ncol(geno_pcs)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(geno_pcs)]
  
  # turn into matrices
  dos_matrix_mat <- as.matrix(dos_matrix[,-1])
  rownames(dos_matrix_mat) <- dos_matrix[[1]]
  exp_matrix_mat <- as.matrix(exp_matrix[,-1])
  rownames(exp_matrix_mat) <- exp_matrix[[1]]
  geno_pcs_mat <- as.matrix(geno_pcs[,-1])
  rownames(geno_pcs_mat) <- geno_pcs[[1]]
  
  # create SlicedData objects
  ## genotype data
  snp_d <- SlicedData$new()    
  snp_d$CreateFromMatrix(dos_matrix_mat)
  rm(dos_matrix, dos_matrix_mat)
  ## expression data
  exp_d <- SlicedData$new()
  exp_d$CreateFromMatrix(exp_matrix_mat)
  rm(exp_matrix, exp_matrix_mat)
  ## covariates data
  cov_d <- SlicedData$new()
  cov_d$CreateFromMatrix(geno_pcs_mat)
  rm(geno_pcs, geno_pcs_mat)
  
  # run main MatrixeQTL function
  for (j in 1:(n_permutations+1)){
    print(j)
    if (j==1){
      # QTL mapping without permutation of gene expression labels
      me <- Matrix_eQTL_main(
        snps = snp_d,
        gene = exp_d,
        cvrt = cov_d,
        pvOutputThreshold = 0,
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = FALSE,
        pvOutputThreshold.cis = 1,
        snpspos = snp_local,
        genepos = gene_local,
        cisDist = 1e6,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
      
      # retrieve and save results
      cis_qtls <- me$cis$eqtls %>% mutate(condition=conditions[i], celltype='PBMC', SE=abs(beta/qnorm(pvalue/2)))
      cis_qtls <- inner_join(cis_qtls, snp_local, by=c('snps'='snpid')) %>% 
        select(snps, chr, pos, gene, statistic, pvalue, FDR, beta, SE, condition, celltype) %>% arrange(chr, pos)
      fwrite(cis_qtls, 'matrixEQTL_results/'%&%conditions[i]%&%'_PBMC_elbowPCs_cisQTL_sumstats.txt', quote=F, sep='\t', na='NA')
      
    } else {
      
      # reorder expression data frame
      permuted_gene <- exp_d$Clone()
      permuted_gene$ColumnSubsample( sample(1:permuted_gene$nCols(), permuted_gene$nCols()) )
      
      # QTL mapping with permuted expression labels
      me <- Matrix_eQTL_main(
        snps = snp_d,
        gene = permuted_gene,
        cvrt = cov_d,
        pvOutputThreshold = 0,
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = FALSE,
        pvOutputThreshold.cis = 1,
        snpspos = snp_local,
        genepos = gene_local,
        cisDist = 1e6,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
      
      # retrieve and save results
      cis_qtls <- me$cis$eqtls %>% mutate(condition=conditions[i], celltype='PBMC', SE=abs(beta/qnorm(pvalue/2)))
      cis_qtls <- inner_join(cis_qtls, snp_local, by=c('snps'='snpid')) %>% 
        select(snps, chr, pos, gene, statistic, pvalue, FDR, beta, SE, condition, celltype) %>% arrange(chr, pos)
      fwrite(cis_qtls, 'matrixEQTL_results/'%&%conditions[i]%&%'_PBMC_Perm'%&%as.character(as.numeric(j)-1)%&%'_elbowPCs_cisQTL_sumstats.txt', quote=F, sep='\t', na='NA')
    }
  }
}
rm(cis_qtls, me, permuted_gene, common_cols, cov_d, exp_d, snp_d)
  
### QVALUE SECTION
permutations <- c(0:10)

for (cond in conditions){
  print(cond)
    for (perm in permutations){
      print(perm)
      
      # select the top SNP for each gene in true and permuted files
      if (perm==0){
        best_true <- fread('matrixEQTL_results/'%&%cond%&%'_PBMC_elbowPCs_cisQTL_sumstats.txt') 
        
        #histogram of unadjusted pvalues
        pdf('matrixEQTL_results/plots/'%&%cond%&%'_PBMC_cisQTL_pvalhistogram.pdf', width=4, height=4)
        hist(best_true$pvalue, main = cond%&%' PBMC cisQTL pvalues', breaks=100)
        dev.off()
        
        # select top SNP per gene based on pvalue
        best_true <- best_true %>% group_by(gene) %>% slice_min(pvalue, with_ties=FALSE)
        best_true <- best_true %>% arrange(gene)
        
      } else {
        tmp <- fread('matrixEQTL_results/'%&%cond%&%'_PBMC_Perm'%&%perm%&%'_elbowPCs_cisQTL_sumstats.txt') 
        
        #histogram of unadjusted pvalues
        pdf('matrixEQTL_results/plots/'%&%cond%&%'_PBMC_Perm'%&%perm%&%'_cisQTL_pvalhistogram.pdf', width=4, height=4)
        hist(tmp$pvalue, main = cond%&%' PBMC Perm'%&%perm%&%' cisQTL pvalues', breaks=100)
        dev.off()
        
        # select top SNP per gene based on pvalue
        tmp <- tmp %>% group_by(gene) %>% slice_min(FDR, with_ties=FALSE)
        tmp <- tmp %>% arrange(gene) %>% ungroup() %>% select(pvalue)
        
        if (exists('compiled.perm')){
          compiled.perm <- cbind(compiled.perm, tmp)
        } else {compiled.perm <- tmp}
      }
    }
    
  # compute qvalues with the top SNPs
  empP <- empPvals(stat=-log10(best_true$pvalue), stat0=-log10(as.matrix(compiled.perm)), pool=TRUE)
  best_true$qvals <- qvalue(empP)$qvalue
  fwrite(best_true, cond%&%'_PBMC_best_cisQTL_sumstats.txt', sep=' ')
    
  #histogram of qvalues
  pdf('matrixEQTL_results/plots/'%&%cond%&%'_PBMC_cisQTL_qvalhistogram.pdf', width=4, height=4)
  hist(best_true$qvals, main = cond%&%' PBMC cisQTL qvalues')
  dev.off()
    
  # qqplot with best true SNPs and best perm SNPs
  best_perm <- data.frame(min_value = apply(compiled.perm, 1, min, na.rm = TRUE))
  pdf('matrixEQTL_results/plots/'%&%cond%&%'_PBMC_best_SNPs_qqplot.pdf', width=4, height=4)
  qqplot(x=-log10(best_perm$min_value), y=-log10(best_true$pvalue), main = cond%&%' PBMC Best SNPs qqplot', 
         xlab='-log10(best permuted p-values)', ylab = '-log10(best true p-values)')
  abline(c(0,1), col='red')
  dev.off()
    
  rm(compiled.perm)
}
