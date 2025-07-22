library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpointdensity)
library(viridis)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')
parser <- ArgumentParser()
parser$add_argument('-c', '--celltype', help='cell type to compile betas if mode is set to "compile"')
parser$add_argument('-m', '--mode', help='"compile" to compile betas or "analyze" to assess correlation')
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
  mutate(snps=Chromosome%&%':'%&%Position) %>% drop_na() %>% select(gene, snps, Ref, Alt, args$celltype)

  # only keep SAIGE genes that are also in haleys results
  genes <- unique(haley_b$gene) %&% '.SAIGE.txt'
  sage_files <- list.files('Saige/step2/outputs/'%&%args$celltype, pattern='\\.txt$')
  sage_files <- sage_files[basename(sage_files) %in% genes]
  rm(genes)

  # get SAIGE beta for each gene-snp pair
  for (f in sage_files){
    # get gene name
    g_name <-  sub('\\.SAIGE\\.txt$', '', f)
  
    # get snp from haley
    snp_id <- haley_b %>% filter(gene==g_name) %>% pull(snps)
  
    # read file
    tmp_f <- fread('Saige/step2/outputs/'%&%args$celltype%&%'/'%&%f) %>% mutate(gene=g_name) %>% 
      filter(MarkerID==snp_id) %>% select(gene, MarkerID, Allele1, Allele2, BETA) %>% 
      rename(snps=MarkerID, Ref=Allele1, Alt=Allele2, args$celltype=BETA)
  
    if (exists('saige_b')){
      saige_b <- rbind(saige_b, tmp_f)
    } else {saige_b <- tmp_f}
  }
} else if (args$mode == 'analyze'){
  
  # # join dfs
  # joint_betas <- inner_join(haley_b, saige_b, by=c('gene', 'snps')) 
  # joint_betas <- joint_betas %>% rename(H_ref=Ref.x, H_alt=Alt.x, H_beta=NK_NI.x, S_ref=Ref.y, S_alt=Alt.y, S_beta=NK_NI.y) %>%
  #   mutate(celltype='NK')
  # 
  # # check if refs and alts match
  # which(joint_betas$H_ref!=joint_betas$S_ref)
  # which(joint_betas$H_alt!=joint_betas$S_alt)
  # 
  # # scatter plot with point density 
  # ggplot(joint_betas, aes(x=H_beta, y=S_beta)) + geom_pointdensity(show.legend=F) +
  #   stat_smooth(method='lm', geom='smooth', formula=y~x) +
  #   geom_abline(slope=1, color='red') + theme_bw() + 
  #   stat_cor(method='pearson', label.x=-1, label.y=1.5, size=3) +
  #   facet_wrap(~celltype) + scale_color_viridis()
  # 
  # ggsave('mashr/eQTLmapping_DanielvsHaley_betas_absbetas_postmash.pdf', height=5, width=10)
  
}