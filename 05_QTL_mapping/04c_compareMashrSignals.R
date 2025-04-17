library(tidyverse)
library(data.table)
library(patchwork)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/mashr')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')

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
  
  # identify genes where at least one condition has a CI that does not include zero
  significant_genes <- tmp %>% filter(lower>0 | upper<0) %>% distinct(gene) %>% pull()

  if (length(significant_genes)>0){
    tmp <- tmp %>% filter(gene %in% significant_genes)
    
    if (nrow(tmp) <= 140){
      # forest plot
      ggplot(tmp, aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
        geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
        geom_vline(xintercept=0, linetype='dashed', color='black') + 
        theme_bw() + labs(x='Effect Size', y='Gene', color='Condition')
      ggsave(celltypes[i]%&%'_eGenes_forestplot.pdf', height=8, width=6)
    
    } else {
      p1 <- ggplot(tmp[1:floor((nrow(tmp)/2)),], aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
          geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
          geom_vline(xintercept=0, linetype='dashed', color='black') + 
          theme_bw() + labs(x='Effect Size', y='Gene', color='Condition') +
          theme(legend.position='none')
      p2 <- ggplot(tmp[(floor((nrow(tmp)/2))+1):nrow(tmp),], aes(y=gene, x=beta, xmin=lower, xmax=upper, color=condition)) +
          geom_pointrange(position=position_dodge(width=0.5), size=0.5) +
          geom_vline(xintercept=0, linetype='dashed', color='black') + 
          theme_bw() + labs(x='Effect Size', y='Gene', color='Condition')

    p1 | p2
    ggsave(celltypes[i]%&%'_eGenes_forestplot.pdf', height=8, width=8)
    }
    

    if(exists('final.df')){
      final.df <- rbind(final.df, tmp)
    } else {final.df <- tmp}
  }
}
final.df <- final.df %>% drop_na()
fwrite(final.df, 'sig_mashr_results.lfsr0.05.txt', col.names=T, sep='\t')
