library(tidyverse)
library(Seurat)
library(data.table)
library(scProportionTest)
library(speckle)
library(miloR)
"%&%" <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('NI', 'IVA', 'RV')
celltypes <- c('B','CD4-T','CD8-T','Mono','NK')

# load seurat object
obj <- readRDS('NI_IVA_RV.integrated.w_celltype.rds')

# adjust metadata
sample_m <- fread('../sample_metadata.txt')
sample_m$income <- ifelse(sample_m$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Low', 'High')
mdata <- left_join(obj@meta.data, sample_m, by=c('IDs'='ID'))
rownames(mdata) <- rownames(obj@meta.data)
obj@meta.data <- mdata
rm(mdata, sample_m)

# create object without batch 4
obj_nob4 <- subset(obj, subset=batch!='B4')  

########################
### scProportionTest ###
########################
for (cond in conditions){
  # create scProp objects
  scprop_obj <- subset(obj, subset=condition==cond) %>% sc_utils()
  scprop_obj_nob4 <- subset(obj_nob4, subset=condition==cond) %>% sc_utils()
  
  # perform test for income
  prop_test <- permutation_test(
    scprop_obj, cluster_identity='celltype',
    sample_1='Low', sample_2='High',
    sample_identity='income')
  permutation_plot(prop_test, order_clusters=F)
  ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_income.png', height=5, width=7)
  prop_test <- prop_test@results %>% bind_rows() %>% 
    rename(reference=Low, alternative=High) %>% 
    mutate(identity='income', condition=cond, batch4='yes')
  if (exists('compiled_scProp')){
    compiled_scProp <- rbind(compiled_scProp, prop_test)
  } else {compiled_scProp <- prop_test}
  
  prop_test <- permutation_test(
    scprop_obj_nob4, cluster_identity='celltype',
    sample_1='Low', sample_2='High',
    sample_identity='income')
  permutation_plot(prop_test, order_clusters=F)
  ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_income_nob4.png', height=5, width=7)
  prop_test <- prop_test@results %>% bind_rows() %>% 
    rename(reference=Low, alternative=High) %>% 
    mutate(identity='income', condition=cond, batch4='no')
  compiled_scProp <- rbind(compiled_scProp, prop_test)
  
  # perform test for asthma
  prop_test <- permutation_test(
    scprop_obj, cluster_identity='celltype',
    sample_1='No', sample_2='Yes',
    sample_identity='asthma')
  permutation_plot(prop_test, order_clusters=F)
  ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_asthma.png', height=5, width=7)
  prop_test <- prop_test@results %>% bind_rows() %>% 
    rename(reference=No, alternative=Yes) %>% 
    mutate(identity='asthma', condition=cond, batch4='yes')
  compiled_scProp <- rbind(compiled_scProp, prop_test)
  
  prop_test <- permutation_test(
    scprop_obj_nob4, cluster_identity='celltype',
    sample_1='No', sample_2='Yes',
    sample_identity='asthma')
  permutation_plot(prop_test, order_clusters=F)
  ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_asthma_nob4.png', height=5, width=7)
  prop_test <- prop_test@results %>% bind_rows() %>% 
    rename(reference=No, alternative=Yes) %>% 
    mutate(identity='asthma', condition=cond, batch4='no')
  compiled_scProp <- rbind(compiled_scProp, prop_test)
  
  # perform test for infection
  if (cond!='NI'){
    # recreate scProp objects
    scprop_obj <- subset(obj, subset= (condition==cond | condition=='NI')) %>% sc_utils()
    scprop_obj_nob4 <- subset(obj_nob4, subset= (condition==cond | condition=='NI')) %>% sc_utils()
    
    prop_test <- permutation_test(
      scprop_obj, cluster_identity='celltype',
      sample_1='NI', sample_2=cond,
      sample_identity='condition')
    permutation_plot(prop_test, order_clusters=F)
    ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_infection.png', height=5, width=7)
    prop_test <- prop_test@results %>% bind_rows() %>% 
      rename(reference=NI, alternative=cond) %>% 
      mutate(identity='infection', condition=cond, batch4='yes')
    compiled_scProp <- rbind(compiled_scProp, prop_test)
    
    prop_test <- permutation_test(
      scprop_obj_nob4, cluster_identity='celltype',
      sample_1='NI', sample_2=cond,
      sample_identity='condition')
    permutation_plot(prop_test, order_clusters=F)
    ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_infection_nob4.png', height=5, width=7)
    prop_test <- prop_test@results %>% bind_rows() %>% 
      rename(reference=NI, alternative=cond) %>% 
      mutate(identity='infection', condition=cond, batch4='no')
    compiled_scProp <- rbind(compiled_scProp, prop_test)
  }
}
# look at results
compiled_scProp$condition <- factor(compiled_scProp$condition, levels=conditions)
sig_compiled_scProp <- compiled_scProp %>% filter(FDR<0.05) %>% 
  mutate(direction=ifelse(obs_log2FD>0, 'Up', 'Down'))

sig_compiled_scProp %>% filter(batch4=='yes') %>% 
  ggplot(., aes(x=clusters, y=identity, size=-log10(FDR), color=abs(obs_log2FD), shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='-log10(FDR)', color='log2FD', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/scProportionTest_sigresults.png', height=5, width=9)

sig_compiled_scProp %>% filter(batch4=='no') %>% 
  ggplot(., aes(x=clusters, y=identity, size=-log10(FDR), color=abs(obs_log2FD), shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='-log10(FDR)', color='log2FD', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/scProportionTest_sigresults_nob4.png', height=5, width=9)


