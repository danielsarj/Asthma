library(tidyverse)
library(data.table)
library(patchwork)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/')

# define all vectors
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

#################
# who are they? #
#################
removed <- c('SEA-3-438', 'SEA-3-427', 'SEA-3-406', 'SEA-3-375', 'SEA-3-369',
             'SEA-3-279', 'SEA-3-223', 'SEA-3-139', 'SEA-3-12')
mdata <- fread('sample_metadata.txt') %>% filter(ID %in% removed)

###############
# DE analysis #
###############
DEresults_withB4 <- fread('DEanalysis/NI_IVAxRV_integrated_limma_results.txt')
DEresults_noB4 <- fread('DEanalysis/NI_IVAxRV_integrated_limma_results_noB4.txt')

# first, compare number of sig. hits per interaction/cell type
sighits_DEresults_withB4 <- DEresults_withB4 %>% filter(sig==TRUE) %>%
  group_by(condition, interaction, celltype) %>% summarise(n_hits=n()) %>% 
  ungroup() %>% mutate(dataset='with_B4')
sighits_DEresults_noB4 <- DEresults_noB4 %>% filter(sig==TRUE) %>%
  group_by(condition, interaction, celltype) %>% summarise(n_hits=n()) %>% 
  ungroup() %>% mutate(dataset='no_B4')

sighits_joint <- rbind(sighits_DEresults_withB4, sighits_DEresults_noB4)
sighits_joint$interaction <- factor(sighits_joint$interaction, levels=c('none', 'income'))
ggplot(sighits_joint, aes(x=celltype, y=n_hits, fill=dataset)) + geom_col(position='dodge') +
  theme_bw() + facet_wrap(~condition+interaction)
ggsave('DEanalysis/withB4_noB4_number.sig.hits.png', height=3, width=8)

# compare logFCs
logFCs_withB4 <- DEresults_withB4 %>% select(Gene, logFC, condition, interaction, celltype, sig) %>%
  rename(logFC_withB4=logFC, sig_withB4=sig)
logFCs_noB4 <- DEresults_noB4 %>% select(Gene, logFC, condition, interaction, celltype, sig) %>%
  rename(logFC_noB4=logFC, sig_noB4=sig)
logFCs_joint <- full_join(logFCs_noB4, logFCs_withB4) %>% 
  mutate(sig = case_when(
    sig_withB4==TRUE & sig_noB4==TRUE ~ 'shared',
    sig_withB4==TRUE & 
      (sig_noB4==FALSE | is.na(sig_noB4)) ~ 'withB4_only',
    (sig_withB4==FALSE | is.na(sig_withB4)) & 
      sig_noB4==TRUE ~ 'noB4_only',
    TRUE ~ 'neither'
  ))
logFCs_joint$interaction <- factor(logFCs_joint$interaction, levels=c('none', 'asthma','income'))
logFCs_joint$condition <- factor(logFCs_joint$condition, levels=c('IVA', 'RV'))
logFCs_joint$sig <- factor(logFCs_joint$sig, levels=c('shared','noB4_only','withB4_only', 'neither'))

summary_logFCs_joint <- logFCs_joint %>%  
  group_by(condition, interaction, celltype, sig) %>% summarise(n_hits=n()) %>% ungroup()
summary_logFCs_joint$sig <- factor(summary_logFCs_joint$sig, levels=c('shared','noB4_only','withB4_only'))
summary_logFCs_joint %>% filter(sig!='neither') %>%
  ggplot(., aes(x=celltype, y=n_hits, fill=sig)) + geom_col(position='dodge') +
  theme_bw() + facet_wrap(~condition+interaction)
ggsave('DEanalysis/summary_withB4_noB4_number.sig.hits.png', height=3, width=8)

(logFCs_joint %>% filter(condition=='IVA') %>%
  ggplot(., aes(x=logFC_withB4, y=logFC_noB4, color=sig)) + geom_point(alpha=0.6) +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
  geom_abline(slope=1) + ggtitle('IVA') + theme(legend.position='none')) +
(logFCs_joint %>% filter(condition=='RV') %>%
  ggplot(., aes(x=logFC_withB4, y=logFC_noB4, color=sig)) + geom_point(alpha=0.6) +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
  geom_abline(slope=1) + ggtitle('RV'))
ggsave('DEanalysis/withB4_noB4_logFCs.png', height=5, width=10)

(logFCs_joint %>% filter(condition=='IVA', sig!='neither') %>%
  ggplot(., aes(x=logFC_withB4, y=logFC_noB4, color=sig)) + geom_point(alpha=0.6) +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
  geom_abline(slope=1) + ggtitle('IVA') + theme(legend.position='none')) +
(logFCs_joint %>% filter(condition=='RV', sig!='neither') %>%
  ggplot(., aes(x=logFC_withB4, y=logFC_noB4, color=sig)) + geom_point(alpha=0.6) +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
  geom_abline(slope=1) + ggtitle('RV'))
ggsave('DEanalysis/withB4_noB4_sig.logFCs.png', height=5, width=10)

###############
#     GSEA    #
###############
GSEAresults_withB4 <- fread('DEanalysis/NI_IVAxRV_integrated_descGSEAresults.txt') 
CAMERAresults_withB4 <- fread('DEanalysis/NI_IVAxRV_integrated_descCAMERAresults.txt') %>% filter(FDR<0.05)
GSEAresults_noB4 <- fread('DEanalysis/NI_IVAxRV_integrated_descGSEAresults_noB4.txt') 
CAMERAresults_nohB4 <- fread('DEanalysis/NI_IVAxRV_integrated_descCAMERAresults_noB4.txt') %>% filter(FDR<0.05)

# first, compare number of sig. hits per interaction/cell type
sigGSEAresults_withB4 <- inner_join(GSEAresults_withB4, CAMERAresults_withB4, 
                                    by=c('pathway', 'condition', 'interaction', 'celltype')) %>%
  group_by(celltype, condition, interaction) %>% filter(padj<0.05, (Direction=='Up' & NES>0) | (Direction=='Down' & NES<0)) %>% 
  summarise(n_hits=n()) %>% mutate(dataset='with_B4')
sigGSEAresults_noB4 <- inner_join(GSEAresults_noB4, CAMERAresults_nohB4, 
                                    by=c('pathway', 'condition', 'interaction', 'celltype')) %>%
  group_by(celltype, condition, interaction) %>% filter(padj<0.05, (Direction=='Up' & NES>0) | (Direction=='Down' & NES<0)) %>% 
  summarise(n_hits=n()) %>% mutate(dataset='no_B4')

sighits_joint <- rbind(sigGSEAresults_withB4, sigGSEAresults_noB4)
sighits_joint$interaction <- factor(sighits_joint$interaction, levels=c('none','asthma','income'))
ggplot(sighits_joint, aes(x=celltype, y=n_hits, fill=dataset)) + geom_col(position='dodge') +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(condition))
ggsave('DEanalysis/withB4_noB4_GSEA_n.sig.hits.png', height=4, width=8)

# compare NES
NES_withB4 <- GSEAresults_withB4 %>% select(pathway, NES, condition, interaction, celltype, padj) %>%
  rename(NES_withB4=NES, padj_withB4=padj)
NES_noB4 <- GSEAresults_noB4 %>% select(pathway, NES, condition, interaction, celltype, padj) %>%
  rename(NES_noB4=NES, padj_noB4=padj)

NES_joint <- full_join(NES_withB4, NES_noB4) %>% 
  mutate(sig = case_when(
    padj_withB4 < 0.05 & padj_noB4 < 0.05 ~ 'shared',
    padj_withB4 < 0.05 & 
      (padj_noB4 >= 0.05 | is.na(padj_noB4)) ~ 'withB4_only',
    (padj_withB4 >= 0.05 | is.na(padj_withB4)) & 
      padj_noB4 < 0.05 ~ 'noB4_only',
    TRUE ~ 'neither'
  ))
NES_joint$interaction <- factor(NES_joint$interaction, levels=c('none', 'asthma','income'))
NES_joint$condition <- factor(NES_joint$condition, levels=c('IVA', 'RV'))
NES_joint$sig <- factor(NES_joint$sig, levels=c('shared','noB4_only','withB4_only', 'neither'))
summary_NES_joint <- NES_joint %>%  
  group_by(condition, interaction, celltype, sig) %>% summarise(n_hits=n()) %>% ungroup()
summary_NES_joint %>% filter(sig!='neither') %>%
  ggplot(., aes(x=celltype, y=n_hits, fill=sig)) + geom_col(position='dodge') +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(condition))
ggsave('DEanalysis/withB4_noB4_GSEA_n.hits.png', height=4, width=9)

(NES_joint %>% filter(condition=='IVA') %>%
  ggplot(., aes(x=NES_withB4, y=NES_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
  geom_abline(slope=1) + ggtitle('IVA') + theme(legend.position='none')) +
(NES_joint %>% filter(condition=='RV') %>%
  ggplot(., aes(x=NES_withB4, y=NES_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
  geom_abline(slope=1) + ggtitle('RV'))
ggsave('DEanalysis/withB4_noB4_GSEA_NES.png', height=5, width=10)

(NES_joint %>% filter(condition=='IVA', sig!='neither') %>%
    ggplot(., aes(x=NES_withB4, y=NES_noB4, color=sig)) + geom_point() +
    theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
    geom_abline(slope=1) + ggtitle('IVA') + theme(legend.position='none')) +
  (NES_joint %>% filter(condition=='RV', sig!='neither') %>%
     ggplot(., aes(x=NES_withB4, y=NES_noB4, color=sig)) + geom_point() +
     theme_bw() + facet_grid(cols=vars(interaction), rows=vars(celltype)) +
     geom_abline(slope=1) + ggtitle('RV'))
ggsave('DEanalysis/withB4_noB4_GSEA_sig.NES.png', height=5, width=10)

###############
# QTL Mapping #
###############
QTLs_withB4 <- fread('QTLmapping/mashr/mashr_out_allstats_df.txt')
QTLs_noB4 <- fread('QTLmapping/mashr/mashr_out_allstats_df_noB4.txt')

# first, compare number of sig. hits per interaction/cell type
sigQTLs_withB4 <- QTLs_withB4 %>% filter(lfsr<0.05) %>% 
  group_by(condition, celltype) %>% mutate(n_hits=n(), dataset='with_B4') %>%
  select(condition, celltype, dataset, n_hits) %>% unique()
sigQTLs_noB4 <- QTLs_noB4 %>% filter(lfsr<0.05) %>% select(condition, celltype) %>%
  group_by(condition, celltype) %>% mutate(n_hits=n(), dataset='no_B4') %>%
  select(condition, celltype, dataset, n_hits) %>% unique()

sighQTLs_joint <- rbind(sigQTLs_withB4, sigQTLs_noB4)
sighQTLs_joint$condition <- factor(sighQTLs_joint$condition, levels=c('NI','IVA','RV'))
ggplot(sighQTLs_joint, aes(x=celltype, y=n_hits, fill=dataset)) + geom_col(position='dodge') +
  theme_bw() + facet_wrap(~condition)
ggsave('QTLmapping/withB4_noB4_n.sig.eGenes.png', height=3, width=9)

# compare effect sizes between the same top gene=snp pair
betaQTLs_withB4 <- QTLs_withB4 %>% select(gene, snps, condition, celltype, beta, lfsr) %>%
  rename(beta_withB4=beta, lfsr_withB4=lfsr)
betaQTLs_noB4 <- QTLs_noB4 %>% select(gene, snps, condition, celltype, beta, lfsr) %>%
  rename(beta_noB4=beta, lfsr_noB4=lfsr)

betaQTLs_joint <- inner_join(betaQTLs_withB4, betaQTLs_noB4) %>% 
  mutate(sig = case_when(
    lfsr_withB4 < 0.05 & lfsr_noB4 < 0.05 ~ 'shared',
    lfsr_withB4 < 0.05 & 
      (lfsr_noB4 >= 0.05 | is.na(lfsr_noB4)) ~ 'withB4_only',
    (lfsr_withB4 >= 0.05 | is.na(lfsr_withB4)) & 
      lfsr_noB4 < 0.05 ~ 'noB4_only',
    TRUE ~ 'neither'
))

summary_betaQTLs_joint <- betaQTLs_joint %>%  
  group_by(condition, celltype, sig) %>% summarise(n_hits=n())
summary_betaQTLs_joint$sig <- factor(summary_betaQTLs_joint$sig, levels=c('shared', 'withB4_only', 'noB4_only', 'neither'))
summary_betaQTLs_joint$condition <- factor(summary_betaQTLs_joint$condition, levels=c('NI','IVA','RV'))
summary_betaQTLs_joint %>% filter(sig!='neither') %>%
  ggplot(., aes(x=celltype, y=n_hits, fill=sig)) + geom_col(position='dodge') +
  theme_bw() + facet_wrap(~condition)
ggsave('QTLmapping/withB4_noB4_n.sig.eGenes_split.png', height=3, width=9)

betaQTLs_joint$condition <- factor(betaQTLs_joint$condition, levels=c('NI','IVA','RV'))
betaQTLs_joint$sig <- factor(betaQTLs_joint$sig, levels=c('shared', 'withB4_only', 'noB4_only', 'neither'))

betaQTLs_joint %>% ggplot(., aes(x=beta_withB4, y=beta_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(slope=1) 
betaQTLs_joint %>% filter(sig!='neither') %>%
  ggplot(., aes(x=beta_withB4, y=beta_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(slope=1) 
ggsave('QTLmapping/withB4_noB4_eGenes_betas.png', height=5, width=10)

# let's check MatrixeQTL outputs
noB4_prefix <- c('IVA_B_14','NI_B_11','RV_B_2','IVA_CD4-T_0','NI_CD4-T_1','RV_CD4-T_1',
                  'IVA_CD8-T_20','NI_CD8-T_1','RV_CD8-T_4','IVA_Mono_3',
                  'NI_Mono_0','RV_Mono_0','IVA_NK_1','NI_NK_15','RV_NK_0')
withB4_prefix <- c('IVA_B_17','NI_B_13','RV_B_17','IVA_CD4-T_8','NI_CD4-T_4','RV_CD4-T_5',
                   'IVA_CD8-T_1','NI_CD8-T_6','RV_CD8-T_2','IVA_Mono_1','NI_Mono_0',
                   'RV_Mono_20','IVA_NK_5','NI_NK_2','RV_NK_2')

# compile results
for (p in noB4_prefix){
  tmp <- fread('QTLmapping/matrixEQTL_results/'%&%p%&%'_best_cisQTL_sumstats_noB4.txt', sep=' ') %>% 
        mutate(id=p) %>% separate(id, c('condition', 'celltype', 'PCs'), '_') %>% 
    select(gene, snps, condition, celltype, beta, qvals) %>%
    rename(beta_noB4=beta, qvals_noB4=qvals)
  if (exists('noB4_bestMatrixeQTL')){
    noB4_bestMatrixeQTL <- rbind(noB4_bestMatrixeQTL, tmp)
    } else {noB4_bestMatrixeQTL <- tmp}
}

for (p in withB4_prefix){
  tmp <- fread('QTLmapping/matrixEQTL_results/'%&%p%&%'_best_cisQTL_sumstats.txt.gz', sep=' ') %>% 
    mutate(id=p) %>% separate(id, c('condition', 'celltype', 'PCs'), '_') %>% 
    select(gene, snps, condition, celltype, beta, qvals) %>%
    rename(beta_withB4=beta, qvals_withB4=qvals)
  if (exists('withB4_bestMatrixeQTL')){
    withB4_bestMatrixeQTL <- rbind(withB4_bestMatrixeQTL, tmp)
  } else {withB4_bestMatrixeQTL <- tmp}
}
rm(tmp)

# compare effect sizes between the same top gene=snp pair before MASH
beta_MatrixeQTLs_joint <- inner_join(noB4_bestMatrixeQTL, withB4_bestMatrixeQTL) %>% 
  mutate(sig = case_when(
    qvals_withB4 < 0.1 & qvals_noB4 < 0.1 ~ 'shared',
    qvals_withB4 < 0.1 & 
      (qvals_noB4 >= 0.1 | is.na(qvals_noB4)) ~ 'withB4_only',
    (qvals_withB4 >= 0.1 | is.na(qvals_withB4)) & 
      qvals_noB4 < 0.1 ~ 'noB4_only',
    TRUE ~ 'neither'
    ))

beta_MatrixeQTLs_joint$condition <- factor(beta_MatrixeQTLs_joint$condition, levels=c('NI','IVA','RV'))
beta_MatrixeQTLs_joint$sig <- factor(beta_MatrixeQTLs_joint$sig, levels=c('shared', 'withB4_only', 'noB4_only', 'neither'))
beta_MatrixeQTLs_joint %>% 
  ggplot(., aes(x=beta_withB4, y=beta_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(slope=1) 
ggsave('QTLmapping/withB4_noB4_eGenes_MatrixeQTL.betas.png', height=5, width=10)
beta_MatrixeQTLs_joint %>% filter(sig!='neither') %>%
  ggplot(., aes(x=beta_withB4, y=beta_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(slope=1) 
ggsave('QTLmapping/withB4_noB4_eGenes_sig.MatrixeQTL.betas.png', height=5, width=10)


#################
# Two Sample MR #
#################
TSMR_withB4 <- fread('QTLmapping/twosampleMR/compiled_mr_results.txt') %>%
  select(exposure, SNP, b, pval, gwas, condition, celltype, effect_allele.exposure) %>%
  rename(b_withB4=b, effect_allele.exposure_withB4=effect_allele.exposure, pval_withB4=pval)
TSMR_noB4 <- fread('QTLmapping/twosampleMR/compiled_mr_results_noB4.txt') %>%
  select(exposure, SNP, b, pval, gwas, condition, celltype, effect_allele.exposure) %>%
  rename(b_noB4=b, effect_allele.exposure_noB4=effect_allele.exposure, pval_noB4=pval)

harmonizedTSMR <- full_join(TSMR_withB4, TSMR_noB4) %>% 
  mutate(sig = case_when(
    pval_withB4 < 0.05 & pval_noB4 < 0.05 ~ 'shared',
    pval_withB4 < 0.05 & 
      (pval_noB4 >= 0.05 | is.na(pval_noB4)) ~ 'withB4_only',
    (pval_withB4 >= 0.05 | is.na(pval_withB4)) & 
      pval_noB4 < 0.05 ~ 'noB4_only',
    TRUE ~ 'neither'
  ))

summary_TSMR <- harmonizedTSMR %>%  
  group_by(condition, celltype, gwas, sig) %>% summarise(n_hits=n())
summary_TSMR$sig <- factor(summary_TSMR$sig, levels=c('shared', 'withB4_only', 'noB4_only', 'neither'))
summary_TSMR$condition <- factor(summary_TSMR$condition, levels=c('NI','IVA','RV'))
summary_TSMR %>% filter(sig!='neither') %>%
  ggplot(., aes(x=celltype, y=n_hits, fill=sig)) + geom_col(position='dodge') +
  theme_bw() + facet_wrap(~condition)
ggsave('QTLmapping/twosampleMR/withB4_noB4_n.sig.2SMR_split.png', height=3, width=9)

harmonizedTSMR <- inner_join(TSMR_withB4, TSMR_noB4) %>% 
  mutate(sig = case_when(
    pval_withB4 < 0.05 & pval_noB4 < 0.05 ~ 'shared',
    pval_withB4 < 0.05 & 
      (pval_noB4 >= 0.05 | is.na(pval_noB4)) ~ 'withB4_only',
    (pval_withB4 >= 0.05 | is.na(pval_withB4)) & 
      pval_noB4 < 0.05 ~ 'noB4_only',
    TRUE ~ 'neither'
  ))
harmonizedTSMR$condition <- factor(harmonizedTSMR$condition, levels=c('NI','IVA','RV'))
harmonizedTSMR$sig <- factor(harmonizedTSMR$sig, levels=c('shared', 'withB4_only', 'noB4_only', 'neither'))
harmonizedTSMR %>% 
  ggplot(., aes(x=b_withB4, y=b_noB4, color=sig)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(slope=1) 
ggsave('QTLmapping/twosampleMR/withB4_noB4_sig.2SMR.betas.png', height=5, width=10)

#################
# Colocalization #
#################
Coloc_withB4 <- fread('QTLmapping/colocalization/best_coloc_results.txt') %>%
  select(gene, snp, PP.H4.abf, condition, celltype, gwas) %>%
  rename(PP.H4.abf_withB4=PP.H4.abf)
Coloc_noB4 <- fread('QTLmapping/colocalization/best_coloc_results_noB4.txt') %>%
  select(gene, snp, PP.H4.abf, condition, celltype, gwas) %>%
  rename(PP.H4.abf_noB4=PP.H4.abf)

Coloc_joint <- full_join(Coloc_withB4, Coloc_noB4)
Coloc_joint$condition <- factor(Coloc_joint$condition, levels=c('NI','IVA','RV'))
Coloc_joint %>% 
  ggplot(., aes(x=PP.H4.abf_withB4, y=PP.H4.abf_noB4, color=gwas)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(slope=1) 
ggsave('QTLmapping/colocalization/withB4_noB4_PPH4s.png', height=5, width=10)
