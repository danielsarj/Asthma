library(Seurat)
library(msigdbr)
library(data.table)
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
