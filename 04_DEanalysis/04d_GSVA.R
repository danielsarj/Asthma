library(GSVA)
library(Seurat)
library(msigdbr)
library(janitor)
library(ggpubr)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')

# get human hallmark gene sets
hallmark_genes <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# load seurat object
objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')

# extract count data
expr_matrix <- GetAssayData(objs, assay='RNA', slot='data') 

# run GSVA
params_gsva <- gsvaParam(expr_matrix, hallmark_genes)
gsva_results <- gsva(params_gsva, verbose=TRUE) %>% as.data.frame() %>% rownames_to_column('Pathways')

# analyze GSVA results
subset_gsva <- gsva_results %>% filter(grepl("INTERFERON", Pathways)) %>% t() %>% 
  as.data.frame() %>% rownames_to_column() %>% row_to_names(1) %>% 
  separate(Pathways, into=c('ID', 'condition', 'celltype'), sep='_')
subset_gsva[,4] <- as.numeric(subset_gsva[,4])
subset_gsva[,5] <- as.numeric(subset_gsva[,5])

for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)

    # subset object
    sub_gsva <- subset_gsva %>% filter(celltype==ctype, (condition=='NI' | condition==cond))
    
    # remove unique indvs
    sub_gsva <- sub_gsva %>% group_by(ID) %>% mutate(is_unique = n() == 1)
    sub_gsva <- sub_gsva %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(sub_gsva$ID)){
      # select indv
      subset_score <- sub_gsva %>% filter(ID==id) 
      
      # confirm position of NI vs infection
      NI_row <- which(subset_score$condition=='NI')
      Inf_row <- which(subset_score$condition==cond)
      
      # compute inf score
      subset_score$IFNa_score <- subset_score$HALLMARK_INTERFERON_ALPHA_RESPONSE[Inf_row] - subset_score$HALLMARK_INTERFERON_ALPHA_RESPONSE[NI_row]
      subset_score$IFNy_score <- subset_score$HALLMARK_INTERFERON_GAMMA_RESPONSE[Inf_row] - subset_score$HALLMARK_INTERFERON_GAMMA_RESPONSE[NI_row]
      
      # edit df
      subset_score <- subset_score %>% filter(condition==cond) %>% 
        select(ID, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.gsva.ifn.scores')){
        paired.gsva.ifn.scores <- rbind(paired.gsva.ifn.scores, subset_score)
      } else {paired.gsva.ifn.scores <- subset_score}
    }
  }
}
# merge w metadata
paired.gsva.ifn.scores <- paired.gsva.ifn.scores %>% full_join(sample_m, by=c('ID'='ID'))

# plot
paired.gsva.ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=IFNa_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('GSVA_IFNascore_pairedID_asthma_boxplots.pdf', height=4, width=7)

paired.gsva.ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=IFNy_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('GSVA_IFNyscore_pairedID_asthma_boxplots.pdf', height=4, width=7)

# run ssGSEA
params_ssgsea <- ssgseaParam(expr_matrix, hallmark_genes)
ssgsea_results <- gsva(params_ssgsea, verbose=TRUE) %>% as.data.frame() %>% rownames_to_column('Pathways')

# analyze GSVA results
subset_ssgsea <- ssgsea_results %>% filter(grepl("INTERFERON", Pathways)) %>% t() %>% 
  as.data.frame() %>% rownames_to_column() %>% row_to_names(1) %>% 
  separate(Pathways, into=c('ID', 'condition', 'celltype'), sep='_')
subset_ssgsea[,4] <- as.numeric(subset_ssgsea[,4])
subset_ssgsea[,5] <- as.numeric(subset_ssgsea[,5])

for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)
    
    # subset object
    sub_ssgsea <- subset_ssgsea %>% filter(celltype==ctype, (condition=='NI' | condition==cond))
    
    # remove unique indvs
    sub_ssgsea <- sub_ssgsea %>% group_by(ID) %>% mutate(is_unique = n() == 1)
    sub_ssgsea <- sub_ssgsea %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(sub_ssgsea$ID)){
      # select indv
      subset_score <- sub_ssgsea %>% filter(ID==id) 
      
      # confirm position of NI vs infection
      NI_row <- which(subset_score$condition=='NI')
      Inf_row <- which(subset_score$condition==cond)
      
      # compute inf score
      subset_score$IFNa_score <- subset_score$HALLMARK_INTERFERON_ALPHA_RESPONSE[Inf_row] - subset_score$HALLMARK_INTERFERON_ALPHA_RESPONSE[NI_row]
      subset_score$IFNy_score <- subset_score$HALLMARK_INTERFERON_GAMMA_RESPONSE[Inf_row] - subset_score$HALLMARK_INTERFERON_GAMMA_RESPONSE[NI_row]
      
      # edit df
      subset_score <- subset_score %>% filter(condition==cond) %>% 
        select(ID, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.ssgsea.ifn.scores')){
        paired.ssgsea.ifn.scores <- rbind(paired.ssgsea.ifn.scores, subset_score)
      } else {paired.ssgsea.ifn.scores <- subset_score}
    }
  }
}
# merge w metadata
paired.ssgsea.ifn.scores <- paired.ssgsea.ifn.scores %>% full_join(sample_m, by=c('ID'='ID'))

# plot
paired.ssgsea.ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=IFNa_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('ssGSEA_IFNascore_pairedID_asthma_boxplots.pdf', height=4, width=7)

paired.ssgsea.ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=IFNy_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('ssGSEA_IFNyscore_pairedID_asthma_boxplots.pdf', height=4, width=7)
