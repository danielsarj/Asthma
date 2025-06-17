library(Seurat)
library(msigdbr)
library(data.table)
library(tidyverse)
library(ggpubr)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('RV', 'IVA', 'NI')
celltypes <- c('B','T-CD4','T-CD8','Mono','NK')

# load sc seurat object
sc_obj <- readRDS('NI_IVA_RV.integrated.w_celltype.rds') %>% 
  SetIdent(value='celltype')
sc_obj@meta.data$condition <- factor(sc_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# get human hallmark gene sets
ifn_genes <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)

# filter for IFN modules 
ifn_genes <- ifn_genes[grep('INTERFERON', names(ifn_genes), value=T)]

# add module score
sc_obj <- AddModuleScore(sc_obj, features=ifn_genes[1], name='IFNa_response')
sc_obj <- AddModuleScore(sc_obj, features=ifn_genes[2], name='IFNy_response')

# viz
#DimPlot(sc_obj, reduction='rna.umap', group.by='celltype', split.by='condition', label=TRUE, label.size=5, repel=TRUE)
FeaturePlot(sc_obj, features='IFNa_response1', label=TRUE, repel=TRUE, split.by='condition') /
FeaturePlot(sc_obj, features='IFNy_response1', label=TRUE, repel=TRUE, split.by='condition')
ggsave('UMAP_IFNresponse.pdf', height=6, width=8)

# retrieve scores per cell type 
for (cond in conditions){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    
    # subset seurat object
    subset_obj <- subset(sc_obj, celltype==ctype & condition==cond)
    
    # retrieve metadata and merge with sample-level metadata
    sub_mdata <- subset_obj@meta.data %>% inner_join(sample_m, by=c('IDs'='ID'))
    
    # subset mdata 
    sub_mdata <- sub_mdata %>% select(IDs, condition, celltype, IFNa_response1, IFNy_response1, asthma, income)
    
    if (exists('ifn.scores')){
      ifn.scores <- rbind(ifn.scores, sub_mdata)
    } else {ifn.scores <- sub_mdata}
  }
}
  
ifn.scores %>% ggplot(., aes(x=condition, y=IFNa_response1, fill=asthma)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNascore_asthma_boxplots.pdf', height=4, width=6)
ifn.scores %>% ggplot(., aes(x=condition, y=IFNy_response1, fill=asthma)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw()   
ggsave('scGSEA_IFNyscore_asthma_boxplots.pdf', height=4, width=6)
ifn.scores %>% mutate(across(where(is.character), ~na_if(.x, ''))) %>% drop_na() %>% 
  ggplot(., aes(x=condition, y=IFNa_response1, fill=income)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNascore_income_boxplots.pdf', height=6, width=10)
ifn.scores %>% mutate(across(where(is.character), ~na_if(.x, ''))) %>% drop_na() %>% 
  ggplot(., aes(x=condition, y=IFNy_response1, fill=income)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNyscore_income_boxplots.pdf', height=6, width=10)


# now indv level 
summary_scores_asthma <- ifn.scores %>% group_by(IDs, condition, celltype, asthma) %>% 
  summarise(IFNa_response=mean(IFNa_response1), IFNy_response=mean(IFNy_response1)) %>% ungroup()
summary_scores_asthma %>% ggplot(., aes(x=condition, y=IFNa_response, fill=asthma)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNascore_byID_asthma_boxplots.pdf', height=4, width=6)
summary_scores_asthma %>% ggplot(., aes(x=condition, y=IFNy_response, fill=asthma)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw()   
ggsave('scGSEA_IFNyscore_byID_asthma_boxplots.pdf', height=4, width=6)

summary_scores_income <- ifn.scores %>% group_by(IDs, condition, celltype, income) %>% 
  summarise(IFNa_response=mean(IFNa_response1), IFNy_response=mean(IFNy_response1)) %>% ungroup()
summary_scores_income %>% mutate(across(where(is.character), ~na_if(.x, ''))) %>% drop_na() %>% 
  ggplot(., aes(x=condition, y=IFNa_response, fill=income)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNascore_byID_income_boxplots.pdf', height=6, width=10)
summary_scores_income %>% mutate(across(where(is.character), ~na_if(.x, ''))) %>% drop_na() %>% 
  ggplot(., aes(x=condition, y=IFNy_response, fill=income)) + geom_boxplot() + 
  facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNyscore_byID_income_boxplots.pdf', height=6, width=10)

### PAIRED IND LEVEL
# retrieve scores per cell type 
for (cond in c('RV','IVA')){
  print(cond)
  for (ctype in celltypes){
    print(ctype)
    
    # subset seurat object
    subset_obj <- subset(sc_obj, celltype==ctype & (condition=='NI' | condition==cond))
    
    # retrieve metadata and merge with sample-level metadata
    sub_mdata <- subset_obj@meta.data %>% inner_join(sample_m, by=c('IDs'='ID'))
    
    # get avg score per individual
    sub_indv <- sub_mdata %>% group_by(IDs, condition, celltype) %>% 
      summarise(IFNa_response=mean(IFNa_response1), IFNy_response=mean(IFNy_response1)) %>% ungroup()
    
    # remove unique indvs
    sub_indv <- sub_indv %>% group_by(IDs) %>% mutate(is_unique = n() == 1)
    sub_indv <- sub_indv %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(sub_indv$IDs)){
      # select indv
      subset_score <- sub_indv %>% filter(IDs==id) 
      
      # confirm position of NI vs infection
      NI_row <- which(subset_score$condition=='NI')
      Inf_row <- which(subset_score$condition==cond)
      
      # compute inf score
      subset_score$IFNa_score <- subset_score$IFNa_response[Inf_row] - subset_score$IFNa_response[NI_row]
      subset_score$IFNy_score <- subset_score$IFNy_response[Inf_row] - subset_score$IFNy_response[NI_row]
      
      # edit df
      subset_score <- subset_score %>% filter(condition==cond) %>% 
        select(IDs, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.ifn.scores')){
        paired.ifn.scores <- rbind(paired.ifn.scores, subset_score)
      } else {paired.ifn.scores <- subset_score}
    }
  }
}
# merge w metadata
paired.ifn.scores <- paired.ifn.scores %>% full_join(sample_m, by=c('IDs'='ID'))

# plot
paired.ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=IFNa_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNascore_pairedID_asthma_boxplots.pdf', height=4, width=6)
paired.ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=IFNy_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() 
ggsave('scGSEA_IFNyscore_pairedID_asthma_boxplots.pdf', height=4, width=6)


paired.ifn.scores %>% filter(condition=='RV') %>% 
  drop_na() %>% ggplot(., aes(x=condition, y=IFNy_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() 



### SAME THING AGAIN BUT SCORING ON ADJUSTED RESIDUALS
# load pseudobulk seurat object
bulk_obj <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds') 
bulk_obj@meta.data$condition <- factor(bulk_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# merge metadata
mdata <- bulk_obj@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
bulk_obj@meta.data <- mdata

for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)
    # subset
    subset_asthma_obj <- subset(bulk_obj, celltype==ctype & 
                                  (condition=='NI' | condition==cond))
  
    # adjust for confounding things for asthma interaction
    subset_asthma_obj <- ScaleData(subset_asthma_obj, 
                        vars.to.regress=c('batch','age','gender','n','albuterol'), 
                        features=rownames(subset_asthma_obj), verbose=FALSE)

    # create an assay using residuals
    residual_assay <- CreateAssayObject(data=GetAssayData(subset_asthma_obj, slot='scale.data'))
    subset_asthma_obj[['residuals']] <- residual_assay

    # run AddModuleScore on the new assay
    subset_asthma_obj <- AddModuleScore(subset_asthma_obj, features=ifn_genes[1], 
                                        name='IFNa_response', assay='residuals')
    subset_asthma_obj <- AddModuleScore(subset_asthma_obj, features=ifn_genes[2], 
                                        name='IFNy_response', assay='residuals')
    
    # retrieve results
    subset_asthma_results <- subset_asthma_obj@meta.data %>% 
      select(IDs, condition, celltype, IFNa_response1, IFNy_response1)
    
    # remove unique indvs
    subset_asthma_results <- subset_asthma_results %>% group_by(IDs) %>% mutate(is_unique = n() == 1)
    subset_asthma_results <- subset_asthma_results %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(subset_asthma_results$IDs)){
      # select indv
      subset_score <- subset_asthma_results %>% filter(IDs==id) 
      
      # confirm position of NI vs infection
      NI_row <- which(subset_score$condition=='NI')
      Inf_row <- which(subset_score$condition==cond)
      
      # compute inf score
      subset_score$IFNa_score <- subset_score$IFNa_response1[Inf_row] - subset_score$IFNa_response1[NI_row]
      subset_score$IFNy_score <- subset_score$IFNy_response1[Inf_row] - subset_score$IFNy_response1[NI_row]
      
      # edit df
      subset_score <- subset_score %>% filter(condition==cond) %>% 
        select(IDs, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, subset_score)
      } else {paired.bulk.ifn.scores <- subset_score}
    }
  }
}
# merge w metadata
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% full_join(sample_m, by=c('IDs'='ID'))

asthma_ifn_scores <- paired.bulk.ifn.scores %>% select(celltype, condition, IFNa_score, IFNy_score, asthma, albuterol) %>%
  drop_na()

# plot
asthma_ifn_scores %>% ggplot(., aes(x=condition, y=IFNa_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('scGSEA_IFNascore_pairedID_adjusted.exp_asthma_boxplots.pdf', height=4, width=6)
asthma_ifn_scores %>% ggplot(., aes(x=condition, y=IFNy_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='wilcox.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('scGSEA_IFNyscore_pairedID_adjusted.exp_asthma_boxplots.pdf', height=4, width=6)
