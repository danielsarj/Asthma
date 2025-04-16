library(MatrixEQTL)
library(tidyverse)
library(data.table)
library(argparse)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

parser <- ArgumentParser()
parser$add_argument('--cond')
parser$add_argument('--ctype')
args <- parser$parse_args()

# load expression matrix and dosage file
# make sure the columns are in the same order
exp_matrix <- fread(args$cond%&%'_'%&%args$ctype%&%'_elbowPCs.txt')
dos_matrix <- fread('../genotypes/imputed_vcfs/imputed_dosage.txt')
tmp_names <- colnames(dos_matrix)
tmp_names <- gsub('SEA3', 'SEA-3', tmp_names)
colnames(dos_matrix) <- tmp_names
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

# load snp and gene location files
snp_local <- fread('../genotypes/imputed_vcfs/snp_location.txt')
gene_local <- fread('gene_location.txt')

# save temporary files
fwrite(dos_matrix_mat, cond%&%'_'%&%ctype%&%'_dosage_TMP.txt', quote=F, sep='\t', row.names=TRUE)
fwrite(exp_matrix_mat, cond%&%'_'%&%ctype%&%'_expression_TMP.txt', quote=F, sep='\t', row.names=TRUE)
fwrite(geno_pcs_mat, cond%&%'_'%&%ctype%&%'_covariates_TMP.txt', quote=F, sep='\t', row.names=TRUE)
rm(dos_matrix, dos_matrix_mat, exp_matrix, exp_matrix_mat, geno_pcs, geno_pcs_mat)

# create SlicedData objects
## genotype data
snp_d <- SlicedData$new()    
snp_d$fileDelimiter='\t'
snp_d$fileOmitCharacters='NA'
snp_d$fileSkipRows=1
snp_d$fileSkipColumns=1
snp_d$fileSliceSize=2000
snp_d$LoadFile(cond%&%'_'%&%ctype%&%'_dosage_TMP.txt')
## expression data
exp_d <- SlicedData$new()
exp_d$fileDelimiter='\t'
exp_d$fileOmitCharacters='NA'
exp_d$fileSkipRows=1
exp_d$fileSkipColumns=1
exp_d$fileSliceSize=2000
exp_d$LoadFile(cond%&%'_'%&%ctype%&%'_expression_TMP.txt')
## covariates data
cov_d <- SlicedData$new()
cov_d$fileDelimiter='\t'
cov_d$fileOmitCharacters='NA'
cov_d$fileSkipRows=1
cov_d$fileSkipColumns=1
cov_d$fileSliceSize=2000
cov_d$LoadFile(cond%&%'_'%&%ctype%&%'_covariates_TMP.txt')

# run main MatrixeQTL function
n_permutations <- 10
for (i in 1:(n_permutations+1)){
  print(i)
  if (i==1){
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
      cisDist = 1e5,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)
    
    # retrieve and save results
    cis_qtls <- me$cis$eqtls %>% mutate(condition=args$cond, celltype=args$ctype, SE=abs(beta/qnorm(pvalue/2)))
    cis_qtls <- inner_join(cis_qtls, snp_local, by=c('snps'='snpid')) %>% 
      select(snps, chr, pos, gene, statistic, pvalue, FDR, beta, SE, condition, celltype) %>% arrange(chr, pos)
    fwrite(cis_qtls, 'matrixEQTL_results/'%&%args$cond%&%'_'%&%args$ctype%&%'_elbowPCs_cisQTL_sumstats.txt', quote=F, sep='\t', na='NA')
    
  } else {
    
    # reorder expression data frame
    permuted_gene <- exp_d$Clone()
    while (sum(colnames(permuted_gene) == colnames(exp_d))>0){
      permuted_gene$ColumnSubsample( sample(1:permuted_gene$nCols(), permuted_gene$nCols()) )
      print(sum(colnames(permuted_gene) == colnames(exp_d)))
    }
    
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
      cisDist = 1e5,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)
    
    # retrieve and save results
    cis_qtls <- me$cis$eqtls %>% mutate(condition=args$cond, celltype=args$ctype, SE=abs(beta/qnorm(pvalue/2)))
    cis_qtls <- inner_join(cis_qtls, snp_local, by=c('snps'='snpid')) %>% 
      select(snps, chr, pos, gene, statistic, pvalue, FDR, beta, SE, condition, celltype) %>% arrange(chr, pos)
    fwrite(cis_qtls, 'matrixEQTL_results/'%&%args$cond%&%'_'%&%args$ctype%&%'_Perm'%&%as.character(as.numeric(i)-1)%&%'_elbowPCs_cisQTL_sumstats.txt', quote=F, sep='\t', na='NA')
  }
}