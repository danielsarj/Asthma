library(tidyverse)
library(data.table)
library(coloc)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/colocalization')

# arguments
parser <- ArgumentParser()
parser$add_argument('--gwas')
parser$add_argument('--eqtl')
args <- parser$parse_args()

# define a helper function to get complement base
comp <- function(allele) {
  chartr('ATCG', 'TAGC', allele)
}

# load gene annotation from ensembl
annotations <- fread('../../DEanalysis/ensembl_genes.txt')

# read GWAS file
gwas_in <- fread(args$gwas) %>% select(hm_chrom, hm_pos, hm_other_allele, hm_effect_allele, hm_effect_allele_frequency, 
                            hm_odds_ratio, standard_error, p_value, n)
  
# add case/control fraction 
if (args$gwas=='FerreiraMAR_COA.h.tsv.gz'){
  gwas_in <- gwas_in %>% mutate(s=13962/300671)
} else {
  gwas_in <- gwas_in %>% mutate(s=28259/572934)
}
  
# load matrixeqtl output
matrix_out <- fread('../matrixEQTL_results/'%&%args$eqtl%&%'PCs_cisQTL_sumstats.txt') %>% 
  filter(!is.na(SE)) %>% separate(snps, into=c('snp_id', 'effectallele'), sep='_')
  
# find sdY per gene
exp_matrix <- fread('../'%&%args$eqtl%&%'PCs.txt') %>% pivot_longer(cols=c(-GENES)) %>% group_by(GENES) %>% 
  summarise(sdY=sd(value)) 
  
# merge info
matrix_out <- matrix_out %>% left_join(exp_matrix, by=c('gene'='GENES'))
rm(exp_matrix)
  
for (wk_gene in (unique(matrix_out$gene))){
  print(c(g, f, wk_gene))
    
  # get gene position 
  gene_pos <- annotations %>% filter(hgnc_symbol==wk_gene)
    
  # subset matrixeqtl output
  sub_matrix <- matrix_out %>% filter(gene==wk_gene)
    
  # subset gwas sumstats, and transform odds ratio into betas
  sub_gwas <- gwas_in %>% filter(hm_chrom==gene_pos$chromosome_name, hm_pos %in% intersect(hm_pos, sub_matrix$pos)) %>%
    mutate(snp_id=hm_chrom%&%':'%&%hm_pos, hm_beta=log(hm_odds_ratio)) %>% arrange(hm_pos)
    
  # re-filter sub_matrix (how to do this better?)
  sub_matrix <- sub_matrix %>% filter(pos %in% intersect(sub_gwas$hm_pos, pos))
      
  if (nrow(sub_gwas)==0 | nrow(sub_matrix)==0){
    next
  } else {
    # harmonize alleles
    sub_gwas <- sub_gwas %>% left_join(sub_matrix[, c('snp_id', 'effectallele')], by = 'snp_id') %>%
      mutate(
        # compute complements
        hm_effect_comp = comp(hm_effect_allele),
        hm_other_comp  = comp(hm_other_allele),
        
        # adjust hm_beta based on allele matching logic
        hm_beta = case_when(
          # perfect match with effect allele → do nothing
          effectallele == hm_effect_allele ~ hm_beta,
          # matches the *other* allele → flip beta
          effectallele == hm_other_allele ~ -hm_beta,
          # matches complement of effect allele (strand flip) → do nothing
          effectallele == hm_effect_comp ~ hm_beta,
          # matches complement of other allele → flip beta
          effectallele == hm_other_comp ~ -hm_beta,
          # otherwise → keep as is for now (will flag below)
          TRUE ~ hm_beta
          )
        )
    
      # identify SNPs with mismatched alleles
      mismatch_snps <- sub_gwas %>% filter(!(effectallele %in% c(hm_effect_allele, hm_other_allele,
        hm_effect_comp, hm_other_comp))) %>% pull(snp_id)
    
      # print any problematic SNPs to check later
      if (length(mismatch_snps) > 0) {
        cat('SNPs with allele mismatches:\n')
        print(mismatch_snps)
      } else {
        cat('All alleles are harmonized correctly.\n')
      }
    
      # create coloc input for eQTL data
      coloc_eqtl_obj <- list(beta=sub_matrix$beta, varbeta=sub_matrix$SE^2, snp=sub_matrix$snp_id,
                             position=sub_matrix$pos, pvalues=sub_matrix$pvalue, type='quant', sdY=sub_matrix$sdY[1]) 
      # create coloc input for GWAS data
      coloc_gwas_obj <- list(beta=sub_gwas$hm_beta, varbeta=sub_gwas$standard_error^2, snp=sub_gwas$snp_id,
                             position=sub_gwas$hm_pos, pvalues=sub_gwas$p_value, type='cc', s=sub_gwas$s[1]) 
      # run coloc
      my.res <- coloc.abf(dataset1=coloc_eqtl_obj, dataset2=coloc_gwas_obj)[['results']] %>% mutate(gene=wk_gene) %>% 
        select(snp, position, gene, SNP.PP.H4)

      # save results
      if (exists('coloc.results')){
        coloc.results <- rbind(coloc.results, my.res)
    } else {coloc.results <- my.res}
  }
}

# save results
eqtl_base <- sub('_[^_]+$', '', args$eqtl)
gwas_base <- sub('_[^_]+$', '', args$gwas)

fwrite(coloc.results, eqtl_base%&%'_'%&%gwas_base%&%'_coloc_results.txt', sep=' ')