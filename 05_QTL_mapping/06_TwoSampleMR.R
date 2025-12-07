library(tidyverse)
library(data.table)
library(TwoSampleMR)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/twosampleMR')

# files
input_prefix <- c('IVA_B_7', 'NI_B_7', 'RV_B_16', 'IVA_Mono_1', 'NI_Mono_2', 'RV_Mono_1', 
                  'IVA_NK_2', 'NI_NK_3', 'RV_NK_4', 'IVA_T-CD4_3', 'NI_T-CD4_2', 'RV_T-CD4_10',
                  'IVA_T-CD8_1', 'NI_T-CD8_6', 'RV_T-CD8_7')
gwas_files <- c('FerreiraMAR_COA.h.tsv.gz', 'SakaueS_COA.h.tsv.gz')

# load MAFs for eQTLs
mafs <- fread('../../genotypes/imputed_vcfs/plink.frq') %>% select(-NCHROBS) %>%
  rename(chr=CHR, snp_id=SNP, effect_allele=A1, other_allele=A2, maf=MAF)

# load gwas sumstats output [outcome]
for (g in gwas_files){
  if (g=='FerreiraMAR_COA.h.tsv.gz'){
    gwas_out <- fread('../colocalization/'%&%g) %>% mutate(snps=hm_chrom%&%':'%&%hm_pos,
                                                           hm_beta=log(hm_odds_ratio), ncase=13962, ncontrol=300671, phenotype='COA') %>% 
      group_by(snps) %>% slice_min(p_value, with_ties=FALSE) %>% ungroup() %>% as.data.frame() %>%
      format_data(type='outcome', snp_col='snps', chr_col='hm_chrom', beta_col='hm_beta', 
                  se_col='standard_error', eaf_col='hm_effect_allele_frequency', 
                  effect_allele_col='hm_effect_allele', other_allele_col='hm_other_allele', 
                  pval_col='p_value', samplesize_col='n', ncase_col='ncase', ncontrol_col='ncontrol', phenotype_col='phenotype')
  } else {
    gwas_out <- fread('../colocalization/'%&%g) %>% rename(hm_chrom=chromosome, hm_pos=base_pair_location, hm_other_allele=other_allele, 
                  hm_effect_allele=effect_allele, hm_effect_allele_frequency=effect_allele_frequency, hm_beta=beta) %>% 
      select(hm_chrom, hm_pos, hm_other_allele, hm_effect_allele, hm_effect_allele_frequency, hm_beta, standard_error, p_value) %>% 
      mutate(ncase=28259, ncontrol=572934, phenotype='COA', snps=hm_chrom%&%':'%&%hm_pos, n=601193) %>% group_by(snps) %>% 
      slice_min(p_value, with_ties=FALSE) %>% ungroup() %>% as.data.frame() %>% format_data(type='outcome', snp_col='snps', 
        chr_col='hm_chrom', beta_col='hm_beta', se_col='standard_error', eaf_col='hm_effect_allele_frequency', 
        effect_allele_col='hm_effect_allele', other_allele_col='hm_other_allele', pval_col='p_value', samplesize_col='n', 
        ncase_col='ncase', ncontrol_col='ncontrol', phenotype_col='phenotype')
  }
  
  for (e in input_prefix){
    # load best matrixeqtl output [exposure]
    matrix_out <- fread('../matrixEQTL_results/'%&%e%&%'_best_cisQTL_sumstats.txt') %>% 
      filter(!is.na(SE)) %>% separate(snps, into=c('snp_id', 'effect_allele'), sep='_') %>%
      inner_join(mafs, by=c('chr', 'snp_id', 'effect_allele')) %>% mutate(n=45) %>% 
      select(chr, pos, snp_id, effect_allele, gene, beta, SE, pvalue, other_allele, maf, n) %>% 
      rename(se=SE, SNP=snp_id, pval=pvalue, eaf=maf, samplesize=n) %>% format_data(type='exposure', phenotype_col='gene')
    
    # harmonise data
    harmonised_df <- harmonise_data(exposure_dat=matrix_out, outcome_dat=gwas_out, action=2)
    
    # perform MR
    mr_res <- mr(harmonised_df) 
    
    # append additional info
    eqtl_base <- sub('_[^_]+$', '', e) %>% str_split_1('_')
    gwas_base <- sub('_[^_]+$', '', g)
    mr_res <- mr_res %>% mutate(gwas=gwas_base, condition=eqtl_base[1], celltype=eqtl_base[2]) %>% 
      inner_join(harmonised_df, by=c('exposure')) %>% select(exposure, method, b, se, pval, gwas,
                                                             condition, celltype, SNP, effect_allele.exposure,
                                                             effect_allele.outcome, beta.exposure, beta.outcome, 
                                                             pval.exposure, pval.outcome)
    
    # merge results
    if (exists('mr.compiled')){
      mr.compiled <- rbind(mr.compiled, mr_res)
    } else {mr.compiled <- mr_res}
  }
}

# save results
mr.compiled$condition <- factor(mr.compiled$condition, levels=c('NI', 'IVA', 'RV'))
fwrite(mr.compiled, 'compiled_mr_results.txt', sep=' ')
mr.compiled$celltype <- gsub('T-CD4', 'CD4-T', mr.compiled$celltype)
mr.compiled$celltype <- gsub('T-CD8', 'CD8-T', mr.compiled$celltype)
## nts: the beta of the Wald ratio method is the ratio of the Boutcome/Bexposure

# find top SNP per facet (condition × celltype × gwas)
top_hits <- mr.compiled %>% group_by(condition, celltype, gwas) %>%
  slice_min(pval, n=1, with_ties=FALSE) %>% ungroup()

# volcano plots
ggplot(mr.compiled, aes(x=b, y=-log10(pval), color=celltype)) + geom_point(size=1, alpha = 0.8) +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', color='grey40') +
  labs(x='Boutcome/Bexposure', y='-log10(p-value)') + theme_bw() + 
  facet_grid(cols=vars(condition), rows=vars(gwas)) +
  geom_text_repel(data=top_hits, aes(label=paste0(exposure,' (',celltype,')')),
                  size=2.8, color='black', segment.color='gray60', box.padding=0.4, point.padding=0.3,
                  max.overlaps=Inf)
ggsave('MR_volcanoplot_allgenes.pdf', height=4, width=8)
ggsave('MR_volcanoplot_allgenes.png', height=4, width=8)

# qqplot
exp_p <- (1:length(mr.compiled$pval))/(length(mr.compiled$pval)+1)
chisq <- qchisq(1-mr.compiled$pval, df=1)
lambda <- median(chisq) / qchisq(0.5, df = 1)

png('MR_pvalues_qqplot_allgenes.png', width=1200, height=1200, res=300)
qqplot(x=-log10(exp_p), y=-log10(mr.compiled$pval), xlab='Expected -log10(p)', ylab='Observed -log10(p)')
abline(0, 1, col='red', lwd=2)
text(x=par()$usr[1] + 0.05 * diff(par()$usr[1:2]), y=par()$usr[4] - 0.05 * diff(par()$usr[3:4]),
  labels=paste0('lambda = ', round(lambda, 3)), adj=c(0, 1), cex=1.2)
dev.off()

# subset to mash-significant eGenes
mash_sig <- fread('../mashr/mashr_out_allstats_df.txt') %>% group_by(condition, celltype) %>%
  filter(lfsr<0.05)
mash_sig$celltype <- gsub('T-CD4', 'CD4-T', mash_sig$celltype)
mash_sig$celltype <- gsub('T-CD8', 'CD8-T', mash_sig$celltype)
sub_mr.compiled <- mr.compiled %>% inner_join(mash_sig, by=c('exposure'='gene', 'condition', 'celltype'))
sub_mr.compiled$condition <- factor(sub_mr.compiled$condition, levels=c('NI', 'IVA', 'RV'))

# find top SNP per facet (condition × celltype × gwas)
sub_top_hits <- sub_mr.compiled %>% group_by(condition, celltype, gwas) %>%
  slice_min(pval, n=1, with_ties=FALSE) %>% ungroup()

# volcano plots
ggplot(sub_mr.compiled, aes(x=b, y=-log10(pval), color=celltype)) + geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', color='grey40') +
  labs(x='Boutcome/Bexposure', y='-log10(p-value)') + theme_bw() + 
  facet_grid(cols=vars(condition), rows=vars(gwas)) +
  geom_text_repel(data=sub_top_hits, aes(label=paste0(exposure,' (',celltype,')')),
                  size=2.8, color='black', segment.color='gray60', box.padding=0.4, point.padding=0.3,
                  max.overlaps=Inf)
ggsave('MR_volcanoplot_mashgenes.pdf', height=4, width=8)
ggsave('MR_volcanoplot_mashgenes.png', height=4, width=8)

# qqplot
exp_p <- (1:length(sub_mr.compiled$pval))/(length(sub_mr.compiled$pval)+1)
chisq <- qchisq(1-sub_mr.compiled$pval, df=1)
lambda <- median(chisq) / qchisq(0.5, df = 1)

png('MR_pvalues_qqplot_mashgenes.png', width=1200, height=1200, res=300)
qqplot(x=-log10(exp_p), y=-log10(sub_mr.compiled$pval), xlab='Expected -log10(p)', ylab='Observed -log10(p)')
abline(0, 1, col='red', lwd=2)
text(x=par()$usr[1] + 0.05 * diff(par()$usr[1:2]), y=par()$usr[4] - 0.05 * diff(par()$usr[3:4]),
     labels=paste0('lambda = ', round(lambda, 3)), adj=c(0, 1), cex=1.2)
dev.off()
