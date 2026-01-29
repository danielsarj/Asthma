library(tidyverse)
library(Seurat)
library(data.table)
library(scProportionTest)
library(speckle)
library(miloR)
library(SingleCellExperiment)
library(scater)
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
rm(sample_m)

########################
### scProportionTest ###
########################
for (cond in conditions){
  for (b in c('yes', 'no')){
    
    # create scProp objects
    if (b=='yes'){
      scprop_obj <- subset(obj, subset= condition==cond) %>% sc_utils()
    } else {
      scprop_obj <- subset(obj, subset= condition==cond & batch!='B4') %>% sc_utils()
    }

    # perform test for income
    prop_test <- permutation_test(
      scprop_obj, cluster_identity='celltype',
      sample_1='Low', sample_2='High',
      sample_identity='income')
    permutation_plot(prop_test, order_clusters=F)
    ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_income_'%&%b%&%'b4.png', height=5, width=7)
    prop_test <- prop_test@results %>% bind_rows() %>% 
      rename(reference=Low, alternative=High) %>% 
      mutate(identity='income', condition=cond, batch4=b)
    if (exists('compiled_scProp')){
      compiled_scProp <- rbind(compiled_scProp, prop_test)
    } else {compiled_scProp <- prop_test}
  
    # perform test for asthma
    prop_test <- permutation_test(
      scprop_obj, cluster_identity='celltype',
      sample_1='No', sample_2='Yes',
      sample_identity='asthma')
    permutation_plot(prop_test, order_clusters=F)
    ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_asthma_'%&%b%&%'b4.png', height=5, width=7)
    prop_test <- prop_test@results %>% bind_rows() %>% 
      rename(reference=No, alternative=Yes) %>% 
      mutate(identity='asthma', condition=cond, batch4=b)
    compiled_scProp <- rbind(compiled_scProp, prop_test)
  
    # perform test for infection
    if (cond!='NI'){
      # recreate scProp objects
      if (b=='yes'){
        scprop_obj <- subset(obj, subset= (condition==cond | condition=='NI')) %>% sc_utils()
      } else {
        scprop_obj <- subset(obj, subset= (condition==cond | condition=='NI') & batch!='B4') %>% sc_utils()
      }

      prop_test <- permutation_test(
        scprop_obj, cluster_identity='celltype',
        sample_1='NI', sample_2=cond,
        sample_identity='condition')
      permutation_plot(prop_test, order_clusters=F)
      ggsave('proportion_plots/scProportionTest_'%&%cond%&%'_infection_'%&%b%&%'b4.png', height=5, width=7)
      prop_test <- prop_test@results %>% bind_rows() %>% 
        rename(reference=NI, alternative=cond) %>% 
        mutate(identity='infection', condition=cond, batch4=b)
      compiled_scProp <- rbind(compiled_scProp, prop_test)
    }
  }
}
rm(scprop_obj, prop_test)
# look at results
compiled_scProp$condition <- factor(compiled_scProp$condition, levels=conditions)
sig_compiled_scProp <- compiled_scProp %>% filter(FDR<0.05) %>% 
  mutate(direction=ifelse(obs_log2FD>0, 'Up', 'Down'))

sig_compiled_scProp %>% filter(batch4=='yes') %>% 
  ggplot(., aes(x=clusters, y=identity, size=-log10(FDR), color=abs(obs_log2FD), shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='-log10(FDR)', color='log2FD', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/scProportionTest_sigresults_yesb4.png', height=5, width=9)

sig_compiled_scProp %>% filter(batch4=='no') %>% 
  ggplot(., aes(x=clusters, y=identity, size=-log10(FDR), color=abs(obs_log2FD), shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='-log10(FDR)', color='log2FD', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/scProportionTest_sigresults_nob4.png', height=5, width=9)

###############
### speckle ###
###############
for (cond in conditions){
  for (b in c('yes', 'no')){
    
    # create objects
    if (b=='yes'){
      subset_obj <- subset(obj, subset=condition==cond) 
    } else {
      subset_obj <- subset(obj, subset= condition==cond & batch!='B4') 
    }
    
    # transform proportions 
    props <- getTransformedProps(subset_obj@meta.data$celltype, subset_obj@meta.data$IDs, transform='logit')
    # subset metadata
    subset_mdata <- mdata %>% filter(condition==cond, IDs %in% colnames(props[[1]])) %>%
      select(IDs, condition, batch, age, gender, asthma, income) %>% unique() %>%
      arrange(IDs)
  
    # perform test for income
    ## create design matrix
    design_matrix <- model.matrix(~0+income+batch+age+gender, data=subset_mdata)
    mycontr <- makeContrasts(incomeHigh-incomeLow, levels=design_matrix)
    ## run propeller
    results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% rename(alternative=PropMean.incomeHigh, reference=PropMean.incomeLow) %>%
      mutate(identity='income', condition=cond, batch4=b)
    ## save results
    if (exists('compiled_speckle')){
      compiled_speckle <- rbind(compiled_speckle, results)
    } else {compiled_speckle <- results}
  
    # perform test for asthma
    ## create design matrix
    design_matrix <- model.matrix(~0+asthma+batch+age+gender, data=subset_mdata)
    mycontr <- makeContrasts(asthmaYes-asthmaNo, levels=design_matrix)
    ## run propeller
    results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% rename(alternative=PropMean.asthmaYes, reference=PropMean.asthmaNo) %>%
      mutate(identity='asthma', condition=cond, batch4=b)
    ## save results
    compiled_speckle <- rbind(compiled_speckle, results)

    # perform test for infection
    if (cond!='NI'){
      # recreate objects
      if (b=='yes'){
        subset_obj <- subset(obj, subset= (condition==cond | condition=='NI')) 
      } else {
        subset_obj <- subset(obj, subset= (condition==cond | condition=='NI') & batch!='B4') 
      }
      
      # transform proportions 
      subset_obj_mdata <- subset_obj@meta.data %>% mutate(identifier=IDs%&%'_'%&%condition)
      props <- getTransformedProps(subset_obj_mdata$celltype, subset_obj_mdata$identifier, transform='logit')
      # subset metadata
      subset_mdata <- subset_obj_mdata %>% select(identifier, condition, batch, age, gender) %>% unique() %>%
        arrange(identifier)
      subset_mdata$condition <- factor(subset_mdata$condition, levels=c('NI', cond))
      # create design matrix
      design_matrix <- model.matrix(~0+condition+batch+age+gender, data=subset_mdata)
      
      if (cond=='IVA'){
        mycontr <- makeContrasts(conditionIVA-conditionNI, levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% rename(alternative=PropMean.conditionIVA, reference=PropMean.conditionNI) %>%
          mutate(identity='infection', condition=cond, batch4=b)
      } else {
        mycontr <- makeContrasts(conditionRV-conditionNI, levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% rename(alternative=PropMean.conditionRV, reference=PropMean.conditionNI) %>%
          mutate(identity='infection', condition=cond, batch4=b)
      }
      
      # save results
      compiled_speckle <- rbind(compiled_speckle, results)
    }
  }
}
rm(subset_obj, results, props, subset_mdata, subset_obj_mdata, mycontr, design_matrix)
# look at results
compiled_speckle$condition <- factor(compiled_speckle$condition, levels=conditions)
sig_compiled_speckle <- compiled_speckle %>% filter(FDR<0.05) %>% 
  mutate(direction=ifelse(PropRatio>1, 'Up', 'Down'))

sig_compiled_speckle %>% filter(batch4=='yes') %>% 
  ggplot(., aes(x=clusters, y=identity, size=-log10(FDR), color=PropRatio, shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='-log10(FDR)', color='PropRatio', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/speckle_sigresults_yesb4.png', height=5, width=9)

sig_compiled_speckle %>% filter(batch4=='no') %>% 
  ggplot(., aes(x=clusters, y=identity, size=-log10(FDR), color=PropRatio, shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='-log10(FDR)', color='PropRatio', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/speckle_sigresults_nob4.png', height=5, width=9)

#############
### miloR ###
#############
for (cond in conditions){
  for (b in c('yes', 'no')){
    
    # create objects
    if (b=='yes'){
      subset_obj <- subset(obj, subset=condition==cond) %>% as.SingleCellExperiment() %>% Milo()
    } else {
      subset_obj <- subset(obj, subset= condition==cond & batch!='B4') %>% as.SingleCellExperiment() %>% Milo()
    }
    
    # construct KNN graph and nhoods
    subset_obj <- buildGraph(subset_obj, k=25, d=30) %>%
      makeNhoods(k=25, d=30, refined=TRUE, prop=0.2)
    ## todo: pick cond-b specific k
    
    # count cells in nhoods
    subset_obj <- countCells(subset_obj, meta.data=data.frame(colData(subset_obj)), sample='IDs')

    # create design matrix
    milo_design <- data.frame(colData(subset_obj))[,c('IDs','asthma')] %>% distinct()
    rownames(milo_design) <- milo_design$IDs
  
    # compute nhood distance
    subset_obj <- calcNhoodDistance(subset_obj, d=30)    
    
    # differential abundance in nhoods
    da_results <- testNhoods(subset_obj, design=~asthma, design.df=milo_design)
    da_results <- annotateNhoods(subset_obj, da_results, 'celltype')
    
    # viz
    subset_obj <- buildNhoodGraph(subset_obj)
    reducedDim(subset_obj, 'UMAP') <- reducedDim(subset_obj, 'RNA.UMAP')
    plotUMAP(subset_obj, colour_by='celltype') + plotNhoodGraphDA(subset_obj, da_results, alpha=0.05) 
    plotDAbeeswarm(da_results, group.by='celltype')
  }
}
