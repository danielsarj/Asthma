library(tidyverse)
library(data.table)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/mashr')

# define condition vector
conditions <- c('NI', 'RV', 'IVA')

# read sig mashr genes
sig_mashr <- fread('sig_mashr_results.lfsr0.05.txt') %>% 
  select(gene, snps, celltype) %>% unique()

# load dosage file
dos_matrix <- fread('../../genotypes/imputed_vcfs/imputed_dosage.txt')

for (i in 1:nrow(sig_mashr)){
  print(i/nrow(sig_mashr)*100)
  
  # subset dosage file for the specific SNP
  subset_dosage <- dos_matrix %>% filter(snpid==sig_mashr$snps[i]) %>% t() %>% 
    as.data.frame() %>% rownames_to_column() %>% row_to_names(row_number=1) %>%
    rename(ID=snpid)
  
  # get expression levels for the specific gene across all conditions
  for (cond in conditions){
    expression <- fread('../'%&%cond%&%'_'%&%sig_mashr$celltype[i]%&%'_elbowPCs.txt') %>%
      filter(GENES==sig_mashr$gene[i]) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% 
      row_to_names(row_number=1) %>% rename(ID=GENES) %>% mutate(condition=cond)
    
    if (ncol(expression)==3){
      if (exists('compiled.exp')){
        compiled.exp <- rbind(compiled.exp, expression)
      } else {compiled.exp <- expression}
    }
  }
  
  # join dosage and expression tbl
  full_tbl <- full_join(subset_dosage, compiled.exp, by=c('ID')) %>% drop_na()
  full_tbl[,2] <- as.factor(full_tbl[,2])
  full_tbl[,3] <- as.numeric(full_tbl[,3])
  full_tbl$condition <- factor(full_tbl$condition, levels=c('NI','IVA','RV'))
  colnames(full_tbl)[2:3] <- c('SNP', 'Gene')
  
  # make boxplot
  ggplot(full_tbl, aes(x=SNP, y=Gene)) + geom_boxplot() + 
    xlab(sig_mashr$snps[i]) + ylab(sig_mashr$gene[i]) +
    theme_bw() + facet_wrap(~condition) 
  
  # save plots
  ggsave('QTL_boxplots/'%&%sig_mashr$celltype[i]%&%'_'%&%sig_mashr$gene[i]%&%'_QTL_boxplot.pdf', 
         height=3, width=6)
  
  rm(compiled.exp)
}
