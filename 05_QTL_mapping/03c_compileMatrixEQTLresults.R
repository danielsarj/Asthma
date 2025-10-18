library(tidyverse)
library(data.table)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')
pcs <- c(1:20)
permutations <- c(0:10)

for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    for (pc in pcs){
      print(pc)
      for (perm in permutations){
        print(perm)
      
        # select the top SNP for each gene in true and permuted files
        if (perm==0){
          best_true <- fread(cond%&%'_'%&%ctype%&%'_'%&%pc%&%'PCs_cisQTL_sumstats.txt') 
        
          #histogram of unadjusted pvalues
          pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_'%&%pc%&%'cisQTL_pvalhistogram.pdf', width=4, height=4)
          hist(best_true$pvalue, main = cond%&%' '%&%ctype%&%' '%&%pc%&%' cisQTL pvalues', breaks=100)
          dev.off()
        
          # select top SNP per gene based on pvalue
          best_true <- best_true %>% group_by(gene) %>% slice_min(pvalue, with_ties=FALSE)
          best_true <- best_true %>% arrange(gene)
        
        } else {
          tmp <- fread(cond%&%'_'%&%ctype%&%'_Perm'%&%perm%&%'_'%&%pc%&%'PCs_cisQTL_sumstats.txt') 
        
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
      fwrite(best_true, cond%&%'_'%&%ctype%&%'_'%&%pc%&%'_best_cisQTL_sumstats.txt', sep=' ')
    
      #histogram of qvalues
      pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_'%&%pc%&%'_cisQTL_qvalhistogram.pdf', width=4, height=4)
      hist(best_true$qvals, main = cond%&%' '%&%ctype%&%' '%&%pc%&%' cisQTL qvalues')
      dev.off()
    
      # qqplot with best true SNPs and best perm1 SNPs
      pdf('plots/'%&%cond%&%'_'%&%ctype%&%'_'%&%pc%&%'_best_SNPs_qqplot.pdf', width=4, height=4)
      qqplot(x=-log10(compiled.perm[,1]), y=-log10(best_true$pvalue), main = cond%&%' '%&%ctype%&%' '%&%pc%&%' Best SNPs qqplot', 
            xlab='-log10(best permuted p-values)', ylab = '-log10(best true p-values)')
      abline(c(0,1), col='red')
      dev.off()
    
      rm(compiled.perm)
    }
  }
}

# assess results after computing qvalues
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    for (pc in pcs){
      print(pc)
      
      result <- fread(cond%&%'_'%&%ctype%&%'_'%&%pc%&%'_best_cisQTL_sumstats.txt', sep=' ') %>% 
        summarise(qval_5e02=sum(qvals<0.05), qval_1e01=sum(qvals<0.1)) %>% 
        mutate(celltype=ctype, condition=cond, n_pcs=pc)

      if (exists('compiled.results')){
        compiled.results <- rbind(compiled.results, result)
      } else {compiled.results <- result}
      
    }
  }
}

# transform dataframe into longer format 
compiled.results_long <- compiled.results %>% pivot_longer(cols=c(qval_5e02, qval_1e01))

# plot number of sig eGenes by nPCs
ggplot(compiled.results_long, aes(x=n_pcs, y=value, color=name)) + geom_line() + 
  geom_point() + theme_bw() + facet_grid(rows=vars(condition), cols=vars(celltype))
ggsave('eGenes_by_nPCs.pdf', height=4, width=10)

# find the smallest PC (in case of ties) that has the largest number of eGenes for qvalue<0.1
compiled.results %>% group_by(celltype, condition) %>% slice_max(qval_1e01, with_ties=TRUE) %>%
  slice_min(n_pcs, with_ties=FALSE)
