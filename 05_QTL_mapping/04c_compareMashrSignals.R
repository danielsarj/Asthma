library(tidyverse)
library(data.table)
library(patchwork)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/mashr')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')
conditions <- c('NI', 'RV', 'IVA')
input_prefix <- c('IVA_B_7', 'NI_B_7', 'RV_B_16', 'IVA_Mono_1', 'NI_Mono_2', 'RV_Mono_1', 
                  'IVA_NK_2', 'NI_NK_3', 'RV_NK_4', 'IVA_T-CD4_3', 'NI_T-CD4_2', 'RV_T-CD4_10',
                  'IVA_T-CD8_1', 'NI_T-CD8_6', 'RV_T-CD8_7')

# load dosage file
dos_matrix <- fread('../../genotypes/imputed_vcfs/imputed_dosage.txt')

# read mashr dfs
beta_df <- fread('mashr_out_beta_df.txt')
sd_df <- fread('mashr_out_sd_df.txt')
lfsr_df <- fread('mashr_out_lfsr_df.txt')

# remove snps in which all lfsr are >0.05
lfsr_df <- lfsr_df %>% filter(!apply(.[, 3:15], 1, function(row) all(row > 0.05)))
beta_df <- left_join(lfsr_df[,1:2], beta_df, by=c('gene', 'snps'))
sd_df <- left_join(lfsr_df[,1:2], sd_df, by=c('gene', 'snps'))

# analyze each cell type separately 
for (i in 1:length(celltypes)){
  
  # get gene-snp pairs in which lfsr < 0.05 for only one condition
  sigs <- lfsr_df %>% select(gene, snps, contains(celltypes[i])) %>%
    filter(rowSums(.[, 3:5] < 0.05) == 1) %>% select(gene, snps)
  
  tmp <- inner_join(beta_df, sd_df, by=c('gene', 'snps')) %>% right_join(sigs, by=c('gene', 'snps')) %>%
    select(gene, snps, contains('_'%&%celltypes[i]%&%'_'), -contains('lfsr'))
  
  tmp <- tmp %>% pivot_longer(cols=ends_with(c('_beta', '_SD')), 
                              names_to=c('condition', '.value'), 
                              names_pattern='(.*)_(beta|SD)') %>% 
    mutate(lower=beta-1.96*SD, upper=beta+1.96*SD, celltype=celltypes[i])
  tmp$condition <- gsub('_'%&%celltypes[i], '', tmp$condition)
  tmp$condition <- factor(tmp$condition, levels=c('NI','IVA','RV'))
  
  # identify genes where at least one condition has a CI that does not include zero
  significant_genes <- tmp %>% filter(lower>0 | upper<0) %>% distinct(gene) %>% pull()

  if (length(significant_genes)>0){
    tmp <- tmp %>% filter(gene %in% significant_genes)
    
    if (length(significant_genes) <= 8){
      # forest plot
      ggplot(tmp, aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
        geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
        geom_vline(xintercept=0, linetype='dashed', color='black') + 
        theme_bw() + labs(x='Effect Size', y='Gene', color='Condition')
      ggsave(celltypes[i]%&%'_eGenes_forestplot.pdf', height=6, width=6)
    
    } else {
      # split df into two for viz purposes
      last_idx <- tmp %>% mutate(row=row_number()) %>% 
        filter(gene == significant_genes[floor(length(significant_genes)/2)]) %>% 
        summarise(last_row = max(row)) %>% pull(last_row)
      tmp_before <- tmp %>% slice(1:last_idx)
      tmp_after  <- tmp %>% slice((last_idx + 1):n())
      
      p1 <- ggplot(tmp_before, aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
          geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
          geom_vline(xintercept=0, linetype='dashed', color='black') + 
          theme_bw() + labs(x='Effect Size', y='Gene', color='Condition') +
          theme(legend.position='none')
      p2 <- ggplot(tmp_after, aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
          geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
          geom_vline(xintercept=0, linetype='dashed', color='black') + 
          theme_bw() + labs(x='Effect Size', y='Gene', color='Condition')

      p1 | p2
      ggsave(celltypes[i]%&%'_eGenes_forestplot.pdf', height=6, width=8)
      }
    
    if(exists('sig_mashr')){
      sig_mashr <- rbind(sig_mashr, tmp)
    } else {sig_mashr <- tmp}
  }
}
sig_mashr <- sig_mashr %>% drop_na() %>% select(gene, snps, celltype) %>% unique()

# box plot
for (i in 1:nrow(sig_mashr)){
  print(i/nrow(sig_mashr)*100)
  
  # subset dosage file for the specific SNP
  subset_dosage <- dos_matrix %>% filter(snpid==sig_mashr$snps[i]) %>% t() %>% 
    as.data.frame() %>% rownames_to_column() %>% row_to_names(row_number=1) %>%
    rename(ID=snpid)
  
  # get expression levels for the specific gene across all conditions
  for (cond in conditions){
    pc <- str_extract(input_prefix[str_detect(input_prefix, paste0('^', cond, '_', sig_mashr$celltype[i], '_'))],'(?<=_)\\d+$')
    
    expression <- fread('../'%&%cond%&%'_'%&%sig_mashr$celltype[i]%&%'_'%&%pc%&%'PCs.txt') %>%
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

# reformat dataframe 
mash_df <- inner_join(beta_df, sd_df) %>% inner_join(lfsr_df)
for (cond in c('NI','IVA','RV')){
  for (ctype in celltypes){
    sub <- mash_df %>% select(gene, snps, contains(cond%&%'_'%&%ctype)) %>% 
      mutate(condition=cond, celltype=ctype)
    colnames(sub) <- gsub(cond%&%'_'%&%ctype%&%'_', '', colnames(sub))
    
    if(exists('mash.long.df')){
      mash.long.df <- rbind(mash.long.df, sub)
    } else {mash.long.df <- sub}
  }
}

# summary
mash_summary <- mash.long.df %>% group_by(condition, celltype) %>% filter(lfsr < 0.05) %>% 
  summarise(n_eGenes=n())
mash_summary$condition <- factor(mash_summary$condition, levels=c('NI','IVA','RV'))
ggplot(mash_summary, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
  theme_bw() + ggtitle('Sig. eGenes')
ggsave('sig_eGenes_per_cond.ctype.pdf', height=4, width=5)

# unique eQTLs per infection
mash_unique <- mash.long.df %>% group_by(gene, snps, celltype) %>%
  filter(sum(lfsr< 0.05)==1) %>% ungroup() %>% filter(lfsr<0.05) %>% 
    group_by(condition, celltype) %>% summarise(n_eGenes=n())
mash_unique$condition <- factor(mash_unique$condition, levels=c('NI','IVA','RV'))
ggplot(mash_unique, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
  theme_bw() + ggtitle('Sig. eGenes that are unique per infection status')
ggsave('sig_eGenes_per_ctype.unique.condition.pdf', height=4, width=5)

# unique eQTLs per celltype
mash_unique <- mash.long.df %>% group_by(gene, snps, condition) %>%
  filter(sum(lfsr< 0.05)==1) %>% ungroup() %>% filter(lfsr<0.05) %>% 
  group_by(condition, celltype) %>% summarise(n_eGenes=n())
mash_unique$condition <- factor(mash_unique$condition, levels=c('NI','IVA','RV'))
ggplot(mash_unique, aes(x=celltype, y=n_eGenes, fill=condition)) + geom_col(position='dodge') +
  theme_bw() + ggtitle('Sig. eGenes that are unique per celltype')
ggsave('sig_eGenes_per_cond.unique.ctype.pdf', height=4, width=5)

# forest plot
mash_filt <- mash.long.df %>% group_by(gene, snps, celltype) %>%
  filter(any(lfsr<0.05) & any(lfsr>0.3)) %>% ungroup()
mash_filt <- mash_filt %>% mutate(lower=beta-1.96*SD, upper=beta+1.96*SD)

ggplot(mash_filt, aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
  geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
  geom_vline(xintercept=0, linetype='dashed', color='black') + 
  theme_bw() + labs(x='Effect Size', y='Gene', color='Condition') +
  facet_wrap(~celltype)

# box plots
mash_filt <- mash_filt %>% select(gene, snps, celltype) %>% unique()
for (i in 1:nrow(mash_filt)){
  print(i/nrow(mash_filt)*100)
  
  # subset dosage file for the specific SNP
  subset_dosage <- dos_matrix %>% filter(snpid==mash_filt$snps[i]) %>% t() %>% 
    as.data.frame() %>% rownames_to_column() %>% row_to_names(row_number=1) %>%
    rename(ID=snpid)
  
  # get expression levels for the specific gene across all conditions
  for (cond in c('NI','IVA','RV')){
    pc <- str_extract(input_prefix[str_detect(input_prefix, paste0('^', cond, '_', sig_mashr$celltype[i], '_'))],'(?<=_)\\d+$')
    
    expression <- fread('../'%&%cond%&%'_'%&%sig_mashr$celltype[i]%&%'_'%&%pc%&%'PCs.txt') %>%
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
    xlab(mash_filt$snps[i]) + ylab(mash_filt$gene[i]) +
    theme_bw() + facet_wrap(~condition) 
  
  # save plots
  ggsave('QTL_boxplots/'%&%mash_filt$celltype[i]%&%'_'%&%mash_filt$gene[i]%&%'_QTL_boxplot.pdf', 
         height=3, width=6)
  
  rm(compiled.exp)
}
