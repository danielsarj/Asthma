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
parser$add_argument('-m', '--mode', help='"compile", "compile_best", or analyze"')
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
}