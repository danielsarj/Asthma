library(tidyverse)
library(data.table)
library(reshape2)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

ctype <- list('B', 'NK')
Gene <- 'CCL2'

for (i in 1:length(ctype)){
  # load expression matrix
  iva_exp_matrix <- fread('IVA_'%&%ctype[i]%&%'_rinResiduals.txt') %>% filter(GENES==Gene) %>%
    melt() %>% mutate(condition='IVA', celltype=as.character(ctype[i]))
  
  rv_exp_matrix <- fread('RV_'%&%ctype[i]%&%'_rinResiduals.txt') %>% filter(GENES==Gene) %>%
    melt() %>% mutate(condition='RV', celltype=as.character(ctype[i]))
  
  ni_exp_matrix <- fread('NI_'%&%ctype[i]%&%'_rinResiduals.txt') %>% filter(GENES==Gene) %>%
    melt() %>% mutate(condition='NI', celltype=as.character(ctype[i]))

  if (exists('exp_matrix')){
    exp_matrix <- rbind(exp_matrix, iva_exp_matrix, rv_exp_matrix, ni_exp_matrix)
  } else {
  exp_matrix <- rbind(iva_exp_matrix, rv_exp_matrix, ni_exp_matrix)
  }
}
  
# retrieve gene-associated snp
snp <- fread('mashr_out_beta_df.txt') %>% filter(gene==Gene) %>% pull(snps)
chrom <-  str_extract(snp, '^[^:]+')
  
# get dosage
dos_matrix <- fread('../genotypes/imputed_vcfs/imputed_chr'%&%chrom%&%'_dosage.txt') %>% 
  filter(snpid==snp) %>% select(-snpid) %>% melt()
dos_matrix$value <- factor(dos_matrix$value, levels=c('0','1','2'))

# combine tables
combined_tbl <- inner_join(exp_matrix, dos_matrix, by=c('variable'))

ggplot(combined_tbl, aes(x=value.y, y=value.x)) + geom_boxplot() + 
  facet_grid(rows=vars(celltype), cols=vars(condition)) + theme_bw() +
  geom_jitter(color='black', size=0.1, alpha=0.4) +
  xlab('Dosage ('%&%snp%&%')') + ylab(Gene)
ggsave(Gene%&%'_boxplot.pdf', width=5, height=4)
rm(exp_matrix)
