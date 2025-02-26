library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')
chromosomes <- seq(1:22)

for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    for (chrom in chromosomes){
      print(chrom)
      
      # read file
      tmp <- fread(cond%&%'_'%&%ctype%&%'_'%&%chrom%&%'_cisQTL_sumstats.txt')
      
      # save into compiled dataframe
      if (exists('full_results')){
        full_results <- rbind(full_results, tmp)
      } else {full_results <- tmp}
    }
    # save compiled results
    fwrite(full_results, cond%&%'_'%&%ctype%&%'_cisQTL_sumstats.txt', sep=' ')
    rm(full_results)
  }
}
