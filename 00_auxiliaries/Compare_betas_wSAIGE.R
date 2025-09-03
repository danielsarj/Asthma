library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpointdensity)
library(viridis)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')
parser <- ArgumentParser()
parser$add_argument('-c', '--celltype', help='cell type to compile betas if mode is set to any "compile"')
parser$add_argument('-m', '--mode', help='"compile", "compile_best", "compile_all", "analyze", "analyze_best')
args <- parser$parse_args()
celltypes <- c('B_NI','CD4_T_NI','CD8_T_NI','monocytes_NI','NK_NI')

# COMPILE BETAS PER CELL TYPE 
if (args$mode == 'compile'){
  # read haleys betas
  haley_b <- fread('haley_betas.txt') %>% separate(gene_SNP, 
                                                   into=c('gene', 
                                                          'Chromosome', 
                                                          'Position', 'Ref', 'Alt'), 
                                                   sep='_') %>% 
    rename(CD4_T_NI=CD4T_NI, CD8_T_NI=CD8T_NI) %>% mutate(snps=Chromosome%&%':'%&%Position) %>% 
    drop_na() %>% select(gene, snps, Ref, Alt, all_of(args$celltype))
  
  # only keep SAIGE genes that are also in haleys results
  genes <- unique(haley_b$gene) %&% '.SAIGE.txt'
  saige_files <- list.files('Saige/step2/outputs/'%&%args$celltype, pattern='\\.txt$')
  saige_files <- saige_files[basename(saige_files) %in% genes]
  rm(genes)
  
  # get SAIGE beta for each gene-snp pair
  for (f in saige_files){
    # counter to keep track of progress
    (which(saige_files==f)/length(saige_files))*100
    
    # get gene name
    g_name <-  sub('\\.SAIGE\\.txt$', '', f)
    
    # get snp from haley
    snp_id <- haley_b %>% filter(gene==g_name) %>% pull(snps)
    
    # read file
    tmp_f <- fread('Saige/step2/outputs/'%&%args$celltype%&%'/'%&%f) %>% mutate(gene=g_name) %>% 
      filter(MarkerID==snp_id) %>% select(gene, MarkerID, Allele1, Allele2, BETA) %>% 
      rename(snps=MarkerID, Ref=Allele1, Alt=Allele2)
    
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
  
  # save file
  fwrite(saige_b, 'Saige/'%&%args$celltype%&%'_compiled_betas.txt', col.names=T)
  
  # COMPILE BEST SAIGE'S QTLs
} else if (args$mode == 'compile_best'){
  # read haleys betas
  haley_genes <- fread('haley_betas.txt') %>% separate(gene_SNP, 
                                                   into=c('gene', 
                                                          'Chromosome', 
                                                          'Position', 'Ref', 'Alt'), 
                                                   sep='_') %>% drop_na() %>% pull(gene)
  
  # only keep SAIGE genes that are also in haleys results
  genes <- unique(haley_genes) %&% '.SAIGE.txt'
  saige_files <- list.files('Saige/step2/outputs/'%&%args$celltype, pattern='\\.txt$')
  saige_files <- saige_files[basename(saige_files) %in% genes]
  rm(genes, haley_genes)
  
  # get SAIGE beta for each best gene-snp pair (lowest pval)
  for (f in saige_files){
    # counter to keep track of progress
    (which(saige_files==f)/length(saige_files))*100
    
    # get gene name
    g_name <-  sub('\\.SAIGE\\.txt$', '', f)
    
    # read file
    tmp_f <- fread('Saige/step2/outputs/'%&%args$celltype%&%'/'%&%f) %>% mutate(gene=g_name) %>% 
      slice_min(p.value) %>% select(gene, MarkerID, Allele1, Allele2, BETA) %>% 
      rename(snps=MarkerID, Ref=Allele1, Alt=Allele2)
    
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
  
  # save file
  fwrite(saige_b, 'Saige/'%&%args$celltype%&%'_compiled_best.betas.txt', col.names=T)
  
} else if (args$mode == 'compile_all'){
  # list all saige results
  saige_files <- list.files('Saige/step2/outputs/'%&%args$celltype, pattern='\\.txt$')

  # get SAIGE beta for each best gene-snp pair (lowest pval)
  for (f in saige_files){
    # counter to keep track of progress
    (which(saige_files==f)/length(saige_files))*100
    
    # get gene name
    g_name <-  sub('\\.SAIGE\\.txt$', '', f)
    
    # read file
    tmp_f <- fread('Saige/step2/outputs/'%&%args$celltype%&%'/'%&%f) %>% mutate(gene=g_name) %>% 
      slice_min(p.value, with_ties=F) %>% select(gene, MarkerID, Allele1, Allele2, BETA, p.value) %>% 
      rename(snps=MarkerID, Ref=Allele1, Alt=Allele2)
    
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
  
  # save file
  fwrite(saige_b, 'Saige/'%&%args$celltype%&%'_compiled_all.eQTLs.txt', col.names=T)
  
} else if (args$mode == 'analyze'){
  # read haleys betas
  haley_b <- readRDS('HALEYs/eQTL_QN_combined.rds')[[1]] %>% rownames_to_column()
  haley_cols <- colnames(haley_b)
  haley_cols <- gsub('CD8T', 'CD8_T', haley_cols)
  haley_cols <- gsub('CD4T', 'CD4_T', haley_cols)
  colnames(haley_b) <- haley_cols
  haley_b <- haley_b %>% separate(rowname, into=c('gene','chr','pos','Ref','Alt'), sep='_') %>%
    mutate(snps=chr%&%':'%&%pos) %>% select(gene, snps, Ref, Alt, contains('NI')) %>%
    pivot_longer(cols=starts_with('beta_'), names_to='celltype', values_to='H_beta',
                 names_transform=list(celltype = ~ sub('^beta_', '', .)))
  
  # read saige's betas
  for (ct in celltypes){
    tmp_f <- fread('Saige/'%&%ct%&%'_compiled_betas.txt') %>% rename(S_beta=BETA) %>% mutate(celltype=ct)
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
  rm(tmp_f)
  
  # join data frames
  joint_b <- inner_join(haley_b, saige_b, by=c('gene', 'snps', 'Ref', 'Alt', 'celltype'))
  
  # plot 
  ggplot(joint_b, aes(x=H_beta, y=S_beta)) + geom_pointdensity(show.legend=F) +
    stat_smooth(method='lm', geom='smooth', formula=y~x) +
    geom_abline(slope=1, color='red') + theme_bw() + 
    stat_cor(method='pearson', label.x=-1, label.y=1.5, size=3) +
    facet_wrap(~celltype) + scale_color_viridis()
  ggsave('Saige/Haleys_pseudobulk.vs.saige_betas_premash.pdf', height=5, width=8)
  
} else if (args$mode == 'analyze_best'){
  # read haleys betas
  haley_b <- readRDS('HALEYs/eQTL_QN_combined.rds')[[1]] %>% rownames_to_column()
  haley_cols <- colnames(haley_b)
  haley_cols <- gsub('CD8T', 'CD8_T', haley_cols)
  haley_cols <- gsub('CD4T', 'CD4_T', haley_cols)
  colnames(haley_b) <- haley_cols
  haley_b <- haley_b %>% separate(rowname, into=c('gene','chr','pos','Ref','Alt'), sep='_') %>%
    mutate(snps=chr%&%':'%&%pos) %>% select(gene, snps, contains('NI')) %>%
    pivot_longer(cols=starts_with('beta_'), names_to='celltype', values_to='H_beta',
                 names_transform=list(celltype = ~ sub('^beta_', '', .)))
  
  # read saige's best betas
  for (ct in celltypes){
    tmp_f <- fread('Saige/'%&%ct%&%'_compiled_best.betas.txt') %>% rename(S_beta=BETA) %>% mutate(celltype=ct) %>%
      select(-c(Ref, Alt))
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
  rm(tmp_f)
  
  # join data frames
  anti_b <- anti_join(saige_b, haley_b, by=c('gene', 'snps', 'celltype')) %>% group_by(gene, celltype) %>%
    slice(1) %>% ungroup()
  joint_b <- inner_join(saige_b, haley_b, by=c('gene', 'snps', 'celltype')) 
  
  # plot intersection of best betas and Haley's
  ggplot(joint_b, aes(x=H_beta, y=S_beta)) + geom_point() +
    stat_smooth(method='lm', geom='smooth', formula=y~x) +
    geom_abline(slope=1, color='red') + theme_bw() + 
    stat_cor(method='pearson', label.x=-1, label.y=1.5, size=3) +
    facet_wrap(~celltype) 
  ggsave('Saige/Haleys_pseudobulk.vs.saige_best.betas_premash_intersection.pdf', height=5, width=8)
  
} else if (args$mode == 'analyze_all'){

  # read saige's betas
  for (ct in celltypes){
    tmp_f <- fread('Saige/'%&%ct%&%'_compiled_all.eQTLs.txt') %>% rename(S_beta=BETA) %>% mutate(celltype=ct) %>%
      select(-c(Ref, Alt))
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
  rm(tmp_f)
  
  # read matrixeqtl's betas
  for (ct in c('B','Mono','NK','T-CD4','T-CD8')){
    tmp_f <- fread('HALEYs/matrixEQTL_results/NI_'%&%ct%&%'_cisQTL_sumstats.txt') %>% rename(M_beta=beta) %>% 
      group_by(gene, celltype) %>% slice_min(pvalue, with_ties=F) %>% select(gene, snps, M_beta, pvalue, celltype)
    if (exists('matrixeqtl_b')){
      matrixeqtl_b <- rbind(matrixeqtl_b, tmp_f)
    } else {matrixeqtl_b <- tmp_f}
  }
  rm(tmp_f)
  
  # make sure labels match
  saige_b$celltype <- gsub('monocytes_NI', 'Mono', saige_b$celltype)
  saige_b$celltype <- gsub('CD4_T_NI', 'T-CD4', saige_b$celltype)
  saige_b$celltype <- gsub('CD8_T_NI', 'T-CD8', saige_b$celltype)
  saige_b$celltype <- gsub('B_NI', 'B', saige_b$celltype)
  saige_b$celltype <- gsub('NK_NI', 'NK', saige_b$celltype)
  
  # combine dfs
  combined_dfs <- inner_join(saige_b, matrixeqtl_b, by=c('gene', 'celltype'))
  ggplot(combined_dfs, aes(x=-log10(pvalue), y=-log10(p.value))) + geom_point() + theme_bw() +
    facet_wrap(~celltype) + xlab('MatrixeQTL -log10(p)') + ylab('SAIGE-QTL -log10(p)') +
    geom_abline(slope=1, color='blue')
  ggsave('Saige/SAIGE_vs_MatrixeQTL_pvals_bestallQTLs_scatterplot.pdf', height=5, width=8)
  
  qqplot(-log10(saige_b$p.value), -log10(matrixeqtl_b$pvalue), xlab='SAIGE-QTL -log10(p)', ylab='MatrixeQTL -log10(p)')
  abline(0, 1, col='red')
    
  # histogram of pvalues
  tmp_m <- combined_dfs %>% select(pvalue, celltype) %>% mutate(method='matrixeqtl')
  tmp_s <- combined_dfs %>% select(p.value, celltype) %>% mutate(method='saige-qtl') %>%
    rename(pvalue=p.value)
  tmp_combined <- rbind(tmp_m, tmp_s)
  ggplot(tmp_combined, aes(x=pvalue)) + geom_histogram() + theme_bw() +
    facet_grid(cols=vars(celltype), rows=vars(method))
  ggsave('Saige/SAIGE_vs_MatrixeQTL_pvals_bestallQTLs_histogram.pdf', height=5, width=8)
  
  # summary table
  tmp_combined %>% group_by(celltype, method) %>% summarise(n_sigs=sum(pvalue<0.05), prop_sigs=sum(pvalue<0.05)/n())
  tmp_combined %>% group_by(celltype, method) %>% summarise(n_sigs=sum(pvalue<1e-5), prop_sigs=sum(pvalue<1e-5)/n())
  
}
