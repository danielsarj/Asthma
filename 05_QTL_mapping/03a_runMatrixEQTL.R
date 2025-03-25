library(MatrixEQTL)
library(tidyverse)
library(data.table)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

parser <- ArgumentParser()
parser$add_argument('--cond')
parser$add_argument('--ctype')
args <- parser$parse_args()

# load expression matrix and dosage file
# make sure the columns are in the same order
exp_matrix <- fread(args$cond%&%'_'%&%args$ctype%&%'_elbowPCs_rinResiduals.txt')
dos_matrix <- fread('../genotypes/imputed_vcfs/imputed_dosage.txt')
common_cols <- intersect(names(exp_matrix), names(dos_matrix))  
dos_matrix <- cbind(dos_matrix[,1], dos_matrix[, ..common_cols])
setcolorder(dos_matrix, c(names(dos_matrix)[1], setdiff(common_cols, names(dos_matrix)[1])))

# turn into matrices
dos_matrix_mat <- as.matrix(dos_matrix[,-1])
rownames(dos_matrix_mat) <- dos_matrix[[1]]
exp_matrix_mat <- as.matrix(exp_matrix[,-1])
rownames(exp_matrix_mat) <- exp_matrix[[1]]

# load snp and gene location files
snp_local <- fread('../genotypes/imputed_vcfs/snp_location.txt')
gene_local <- fread('gene_location.txt')

# create SlicedData objects
## genotype data
snp_d <- SlicedData$new()    
snp_d$CreateFromMatrix(dos_matrix_mat)
rm(dos_matrix, dos_matrix_mat)
## expression data
exp_d <- SlicedData$new()
exp_d$CreateFromMatrix(exp_matrix_mat)
rm(exp_matrix, exp_matrix_mat)

# run main MatrixeQTL function
n_permutations <- 10
for (i in 1:n_permutations+1){
  print(i)
  if (i==1){
    # QTL mapping without permutation of gene expression labels
    me <- Matrix_eQTL_main(
      snps = snp_d,
      gene = exp_d,
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
    cis_qtls <- me$cis$eqtls %>% mutate(condition=args$cond, celltype=args$ctype, SE=abs(beta/qnorm(pvalue/2)))
    cis_qtls <- inner_join(cis_qtls, snp_local, by=c('snps'='snpid')) %>% 
      select(snps, chr, pos, gene, statistic, pvalue, FDR, beta, SE, condition, celltype) %>% arrange(chr, pos)
    fwrite(cis_qtls, 'matrixEQTL_results/'%&%args$cond%&%'_'%&%args$ctype%&%'_elbowPCs_cisQTL_sumstats.txt', quote=F, sep='\t', na='NA')
    
  } else {
    
    # reorder expression data frame
    permuted_gene <- exp_d$Clone()
    permuted_gene$ColumnSubsample( sample(1:permuted_gene$nCols(), permuted_gene$nCols()) )
    
    # QTL mapping with permuted expression labels
    me <- Matrix_eQTL_main(
      snps = snp_d,
      gene = permuted_gene,
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
    cis_qtls <- me$cis$eqtls %>% mutate(condition=args$cond, celltype=args$ctype, SE=abs(beta/qnorm(pvalue/2)))
    cis_qtls <- inner_join(cis_qtls, snp_local, by=c('snps'='snpid')) %>% 
      select(snps, chr, pos, gene, statistic, pvalue, FDR, beta, SE, condition, celltype) %>% arrange(chr, pos)
    fwrite(cis_qtls, 'matrixEQTL_results/'%&%args$cond%&%'_'%&%args$ctype%&%'_Perm'%&%i-1%&%'_elbowPCs_cisQTL_sumstats.txt', quote=F, sep='\t', na='NA')
  }
}