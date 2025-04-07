library(tidyverse)
library(data.table)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/matrixEQTL_results')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)

    tmp <- fread(cond%&%'_'%&%ctype%&%'_best_cisQTL_sumstats.txt') %>%
      filter(qvals<0.1)
    
    if (nrow(tmp)>0){
      if (exists('compiled.QTL')){
        compiled.QTL <- rbind(compiled.QTL, tmp)
      } else {compiled.QTL <- tmp}
    }
  }
}

# load dosage file
dos_matrix <- fread('../../genotypes/imputed_vcfs/imputed_dosage.txt')
short.qtl <- compiled.QTL %>% select(snps, gene, celltype) %>% unique()

summ_shortqtl <- short.qtl %>% group_by(celltype) %>% summarise(n=n())
ggplot(summ_shortqtl, aes(x=celltype, y=n)) + geom_col() + 
  geom_text(aes(label=n), position=position_dodge(width=0.9),
            vjust=-0.5, size=4) + theme_bw()
ggsave('plots/QTL_boxplots/eGenes_per_celltype_barplot.pdf', height=3, width=4)

for (i in 1:nrow(short.qtl)){
  print(i/nrow(short.qtl)*100)
  
  # subset dosage file for the specific SNP
  subset_dosage <- dos_matrix %>% filter(snpid==short.qtl$snps[i]) %>% t() %>% 
    as.data.frame() %>% rownames_to_column() %>% row_to_names(row_number=1) %>%
    rename(ID=snpid)
  
  # get expression levels for the specific gene across all conditions
  for (cond in conditions){
    expression <- fread('../'%&%cond%&%'_'%&%short.qtl$celltype[i]%&%'_elbowPCs.txt') %>%
      filter(GENES==short.qtl$gene[i]) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% 
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
  colnames(full_tbl)[2:3] <- c('SNP', 'Gene')

  # make boxplot
  ggplot(full_tbl, aes(x=SNP, y=Gene)) + geom_boxplot() + 
    xlab(short.qtl$snps[i]) + ylab(short.qtl$gene[i]) +
    theme_bw() + facet_wrap(~condition) 
  
  # save plots
  ggsave('plots/QTL_boxplots/'%&%short.qtl$celltype[i]%&%'_'%&%short.qtl$gene[i]%&%'_QTL_boxplot.pdf', 
         height=4, width=12)
  
  rm(compiled.exp)
}
