library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/')
celltypes <- c('NK_NI')
pcs <- c('PCs1', 'PCs2', 'PCs3', 'PCs4', 'PCs5')

for (ct in celltypes){
  print(ct)
  for (pc in pcs){
    print(pc)
    # find all SAIGE-QTL step3 outputs
    saige_files <- list.files(ct%&%'/exp_'%&%pc, full.names=T)
  
    for (f in saige_files){
      # read file
      tmp <- fread(f) %>% mutate(celltype=ct, exp_pcs=pc)
    
      # combine into a sigle dataframe
      if (exists('final_df')){
      final_df <- rbind(final_df, tmp)
      } else {final_df <- tmp}
    }
  }
}
fwrite(final_df, 'compiled_step3_outputs.txt', sep=' ', col.names=T)

summary_df <- final_df %>% group_by(celltype, exp_pcs) %>% summarise(n_eGenes=sum(top_pval<1e-5))
ggplot(summary_df, aes(x=exp_pcs, y=n_eGenes, group=celltype)) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(~celltype)