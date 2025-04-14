library(tidyverse)
library(data.table)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')
permutations <- c(0:10)

for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    for (perm in permutations){
      print(perm)
      
      # select the top SNP for each gene in true and permuted files
      if (perm==0){
        best_true <- fread(cond%&%'_'%&%ctype%&%'_elbowPCs_cisQTL_sumstats.txt') 
        
        #histogram of unadjusted pvalues
        pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_cisQTL_pvalhistogram.pdf', width=4, height=4)
        hist(best_true$pvalue, main = cond%&%' '%&%ctype%&%' cisQTL pvalues', breaks=100)
        dev.off()
        
        # select top SNP per gene based on pvalue
        best_true <- best_true %>% group_by(gene) %>% slice_min(pvalue, with_ties=FALSE)
        best_true <- best_true %>% arrange(gene)
        
      } else {
        tmp <- fread(cond%&%'_'%&%ctype%&%'_Perm'%&%perm%&%'_elbowPCs_cisQTL_sumstats.txt') 
        
        #histogram of unadjusted pvalues
        pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_Perm'%&%perm%&%'_cisQTL_pvalhistogram.pdf', width=4, height=4)
        hist(tmp$pvalue, main = cond%&%' '%&%ctype%&%' Perm'%&%perm%&%' cisQTL pvalues', breaks=100)
        dev.off()
        
        # select top SNP per gene based on pvalue
        tmp <- tmp %>% group_by(gene) %>% slice_min(pvalue, with_ties=FALSE)
        tmp <- tmp %>% arrange(gene) %>% ungroup() %>% select(pvalue)
        
        if (exists('compiled.perm')){
          compiled.perm <- cbind(compiled.perm, tmp)
        } else {compiled.perm <- tmp}
      }
    }
    
    # compute qvalues with the top SNPs
    empP <- empPvals(stat=-log10(best_true$pvalue), stat0=-log10(as.matrix(compiled.perm)), pool=TRUE)
    best_true$qvals <- qvalue(empP)$qvalue
    fwrite(best_true, cond%&%'_'%&%ctype%&%'_best_cisQTL_sumstats.txt', sep=' ')
    
    #histogram of qvalues
    pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_cisQTL_qvalhistogram.pdf', width=4, height=4)
    hist(best_true$qvals, main = cond%&%' '%&%ctype%&%' cisQTL qvalues')
    dev.off()
    
    # qqplot with best true SNPs and best perm SNPs
    best_perm <- data.frame(min_value = apply(compiled.perm, 1, min, na.rm = TRUE))
    pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_best_SNPs_qqplot.pdf', width=4, height=4)
    qqplot(x=-log10(best_perm$min_value), y=-log10(best_true$pvalue), main = cond%&%' '%&%ctype%&%' Best SNPs qqplot', 
           xlab='-log10(best permuted p-values)', ylab = '-log10(best true p-values)')
    abline(c(0,1), col='red')
    dev.off()
    
    rm(compiled.perm)
  }
}
