library(tidyverse)
library(data.table)
library(patchwork)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/mashr')
conditions <- c('NI', 'RV', 'IVA')
input_prefix <- c('IVA_B_2',
                  'NI_B_0',
                  'RV_B_0',
                  'IVA_CD4-T_0',
                  'NI_CD4-T_0',
                  'RV_CD4-T_1',
                  'IVA_CD8-T_2',
                  'NI_CD8-T_0',
                  'RV_CD8-T_1',
                  'IVA_Mono_19',
                  'NI_Mono_20',
                  'RV_Mono_0',
                  'IVA_NK_2',
                  'NI_NK_9',
                  'RV_NK_5')
# load dosage file
dos_matrix <- fread('../../genotypes/imputed_vcfs/imputed_dosage.txt')

# read mashr dfs
mash_df <- fread('mashr_out_allstats_df_new.txt')

# remove snps in which all lfsr are >0.05 for a given gene
mash_df <- mash_df %>% group_by(gene) %>% filter(sum(lfsr<0.05)<15) %>% ungroup()

# total number of eGenes per condition and celltype
mash_total <- mash_df %>% group_by(condition, celltype) %>% filter(lfsr<0.05) %>% 
  summarise(n_eGenes=n()) %>% ungroup()
mash_total$condition <- factor(mash_total$condition, levels=c('NI','IVA','RV'))

ggplot(mash_total, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
  theme_bw() + ggtitle('Sig. eGenes')
ggsave('sig_eGenes_per_cond.ctype_new.pdf', height=4, width=5)

# unique eQTLs per infection 
mash_unique_inf <- mash_df %>% group_by(gene, celltype) %>% filter(sum(lfsr<0.05)==1) %>% 
  ungroup() %>% filter(lfsr<0.05) %>% group_by(condition, celltype) %>% summarise(n_eGenes=n()) %>% ungroup()
mash_unique_inf$condition <- factor(mash_unique_inf$condition, levels=c('NI','IVA','RV'))

ggplot(mash_unique_inf, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
  theme_bw() + ggtitle('Sig. eGenes that are unique per infection status')
ggsave('sig_eGenes_per_ctype.unique.condition_new.pdf', height=4, width=5)

# unique eQTLs per celltype
mash_unique_ct <- mash_df %>% group_by(gene, condition) %>% filter(sum(lfsr<0.05)==1) %>% 
  ungroup() %>% filter(lfsr<0.05) %>% group_by(condition, celltype) %>% summarise(n_eGenes=n()) %>% ungroup()
mash_unique_ct$condition <- factor(mash_unique_ct$condition, levels=c('NI','IVA','RV'))

ggplot(mash_unique_ct, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
  theme_bw() + ggtitle('Sig. eGenes that are unique per celltype')
ggsave('sig_eGenes_per_cond.unique.ctype_new.pdf', height=4, width=5)

(ggplot(mash_total, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
    theme_bw() + ggtitle('eGenes') + guides(fill='none')) +
  (ggplot(mash_unique_ct, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
     theme_bw() + ggtitle('Unique eGenes per celltype') + ylab(NULL) + guides(fill='none')) +
  (ggplot(mash_unique_inf, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
     theme_bw() + ggtitle('Unique eGenes per infection status') + ylab(NULL))
ggsave(filename='sig_eGenes_allbarplots_new.pdf', height=3, width=10)
ggsave(filename='sig_eGenes_allbarplots_new.png', height=3, width=10)
rm(mash_total, mash_unique_ct, mash_unique_inf)

# box plots of significant eGenes per celltype (significant in the celltype in at least one infection status)
mash_reduced <- mash_df %>% group_by(gene, condition, celltype) %>% filter(lfsr<0.05) %>% ungroup() %>%
  select(gene, snps, celltype) %>% unique()
for (i in 1:nrow(mash_reduced)){
  # subset dosage file for the specific SNP
  subset_dosage <- dos_matrix %>% filter(snpid==mash_reduced$snps[i]) %>% t() %>% 
    as.data.frame() %>% rownames_to_column() %>% row_to_names(row_number=1) %>%
    rename(ID=snpid)
  
  # get expression levels for the specific gene across all conditions
  for (cond in c('NI','IVA','RV')){
    pc <- str_extract(input_prefix[str_detect(input_prefix, paste0('^', cond, '_', mash_reduced$celltype[i], '_'))],'(?<=_)\\d+$')
    
    expression <- fread('../'%&%cond%&%'_'%&%mash_reduced$celltype[i]%&%'_'%&%pc%&%'PCs_new.txt') %>%
      filter(GENES==mash_reduced$gene[i]) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% 
      row_to_names(row_number=1) %>% rename(ID=GENES) %>% mutate(condition=cond)
    
    if (ncol(expression)==3){
      if (exists('compiled.exp')){
        compiled.exp <- rbind(compiled.exp, expression)
      } else {compiled.exp <- expression}
    }
  }
  
  # join dosage and expression tbl
  full_tbl <- full_join(subset_dosage, compiled.exp, by=c('ID')) %>% drop_na()
  full_tbl[,2] <- as.numeric(full_tbl[,2])
  full_tbl[,3] <- as.numeric(full_tbl[,3])
  full_tbl$condition <- factor(full_tbl$condition, levels=c('NI','IVA','RV'))
  colnames(full_tbl)[2:3] <- c('SNP', 'Gene')
  
  # make boxplot
  ggplot(full_tbl, aes(x=as.factor(SNP), y=Gene)) + geom_violin(fill='lightblue', alpha=0.5) + geom_boxplot(fill='white', width=0.3) +
    xlab(mash_reduced$snps[i]) + ylab(mash_reduced$gene[i]) + theme_bw() + facet_wrap(~condition)
  
  # save plots
  ggsave('QTL_boxplots/'%&%mash_reduced$celltype[i]%&%'_'%&%mash_reduced$gene[i]%&%'_QTL_boxplot_new.pdf', 
         height=3, width=6)
  
  rm(compiled.exp)
}
