library(tidyverse)
library(Seurat)
library(data.table)
library(scProportionTest)
library(speckle)
library(miloR)
library(SingleCellExperiment)
library(ggbeeswarm)
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
      dplyr::rename(reference=Low, alternative=High) %>% 
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
      dplyr::rename(reference=No, alternative=Yes) %>% 
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
        dplyr::rename(reference=NI, alternative=cond) %>% 
        mutate(identity='infection', condition=cond, batch4=b)
      compiled_scProp <- rbind(compiled_scProp, prop_test)
    }
  }
}
rm(scprop_obj, prop_test)
# look at results
compiled_scProp$condition <- factor(compiled_scProp$condition, levels=conditions)
compiled_scProp$identity <- factor(compiled_scProp$identity, levels=c('infection', 'asthma', 'income'))
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
      rownames_to_column(var='clusters') %>% dplyr::rename(alternative=PropMean.incomeHigh, reference=PropMean.incomeLow) %>%
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
      rownames_to_column(var='clusters') %>% dplyr::rename(alternative=PropMean.asthmaYes, reference=PropMean.asthmaNo) %>%
      mutate(identity='asthma', condition=cond, batch4=b)
    ## save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # perform test for infection and interactions
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
      subset_mdata <- subset_obj_mdata %>% select(identifier, condition, batch, age, gender, asthma, income) %>% unique() %>%
        arrange(identifier)
      subset_mdata$condition <- factor(subset_mdata$condition, levels=c('NI', cond))
      subset_mdata$asthma <- factor(subset_mdata$asthma, levels=c('No', 'Yes'))
      subset_mdata$income <- factor(subset_mdata$income, levels=c('Low', 'High'))
      
      # create design matrix for infection
      design_matrix <- model.matrix(~0+condition+batch+age+gender, data=subset_mdata)
      if (cond=='IVA'){
        mycontr <- makeContrasts(conditionIVA-conditionNI, levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% dplyr::rename(alternative=PropMean.conditionIVA, reference=PropMean.conditionNI) %>%
          mutate(identity='infection', condition=cond, batch4=b)
      } else {
        mycontr <- makeContrasts(conditionRV-conditionNI, levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% dplyr::rename(alternative=PropMean.conditionRV, reference=PropMean.conditionNI) %>%
          mutate(identity='infection', condition=cond, batch4=b)
      }
      # save results
      compiled_speckle <- rbind(compiled_speckle, results)
      
      # create design matrix for asthma interaction
      subset_mdata$group <- with(subset_mdata, paste(condition, asthma, sep='_'))
      subset_mdata$group <- factor(subset_mdata$group, levels = c('NI_No', 'NI_Yes', cond%&%'_No', cond%&%'_Yes'))
      design_matrix <- model.matrix(~0+group+batch+age+gender, data=subset_mdata)
      if (cond=='IVA'){
        mycontr <- makeContrasts(interaction=(groupIVA_Yes-groupIVA_No)-(groupNI_Yes-groupNI_No), levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% 
          mutate(alternative=PropMean.groupIVA_Yes-PropMean.groupIVA_No, reference=PropMean.groupNI_Yes-PropMean.groupNI_No,
                 identity='infection:asthma', condition=cond, batch4=b) %>% 
          select(clusters, alternative, reference, PropRatio, Tstatistic, P.Value, FDR, identity, condition, batch4)
      } else {
        mycontr <- makeContrasts(interaction=(groupRV_Yes-groupRV_No)-(groupNI_Yes-groupNI_No), levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% 
          mutate(alternative=PropMean.groupRV_Yes-PropMean.groupRV_No, reference=PropMean.groupNI_Yes-PropMean.groupNI_No,
                 identity='infection:asthma', condition=cond, batch4=b) %>% 
          select(clusters, alternative, reference, PropRatio, Tstatistic, P.Value, FDR, identity, condition, batch4)
      }
      # save results
      compiled_speckle <- rbind(compiled_speckle, results)
      
      # create design matrix for income interaction
      subset_mdata$group <- with(subset_mdata, paste(condition, income, sep='_'))
      subset_mdata$group <- factor(subset_mdata$group, levels = c('NI_Low', 'NI_High', cond%&%'_Low', cond%&%'_High'))
      design_matrix <- model.matrix(~0+group+batch+age+gender, data=subset_mdata)
      if (cond=='IVA'){
        mycontr <- makeContrasts(interaction=(groupIVA_High-groupIVA_Low)-(groupNI_High-groupNI_Low), levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% 
          mutate(alternative=PropMean.groupIVA_High-PropMean.groupIVA_Low, reference=PropMean.groupNI_High-PropMean.groupNI_Low,
                 identity='infection:income', condition=cond, batch4=b) %>% 
          select(clusters, alternative, reference, PropRatio, Tstatistic, P.Value, FDR, identity, condition, batch4)
      } else {
        mycontr <- makeContrasts(interaction=(groupRV_High-groupRV_Low)-(groupNI_High-groupNI_Low), levels=design_matrix)
        # run propeller
        results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
          rownames_to_column(var='clusters') %>% 
          mutate(alternative=PropMean.groupRV_High-PropMean.groupRV_Low, reference=PropMean.groupNI_High-PropMean.groupNI_Low,
                 identity='infection:income', condition=cond, batch4=b) %>% 
          select(clusters, alternative, reference, PropRatio, Tstatistic, P.Value, FDR, identity, condition, batch4)
      }
      # save results
      compiled_speckle <- rbind(compiled_speckle, results)
    }
  }
}
rm(subset_obj, results, props, subset_mdata, subset_obj_mdata, mycontr, design_matrix)
# look at results
compiled_speckle$condition <- factor(compiled_speckle$condition, levels=conditions)
compiled_speckle$identity <- factor(compiled_speckle$identity, levels=c('infection', 'asthma', 'income', 'infection:asthma', 'infection:income'))
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
    subset_obj <- buildGraph(subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA')
    
    # count cells in nhoods
    subset_obj <- countCells(subset_obj, meta.data=data.frame(colData(subset_obj)), sample='IDs')
    
    # create design matrix
    milo_design <- data.frame(colData(subset_obj))[,c('IDs','batch','age','gender','asthma','income')] %>% distinct()
    rownames(milo_design) <- milo_design$IDs
    milo_design$asthma <- factor(milo_design$asthma, levels=c('No', 'Yes'))
    milo_design$income <- factor(milo_design$income, levels=c('Low', 'High'))
    
    # compute nhood distance
    subset_obj <- calcNhoodDistance(subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    
    # differential abundance in nhoods for income
    da_results <- testNhoods(subset_obj, design=~batch+age+gender+income, design.df=milo_design)
    da_results <- annotateNhoods(subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='income', condition=cond, batch4=b,
                                        is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + geom_quasirandom(data=filter(da_results, !is.sig),
           aes(x=logFC), color='grey', alpha=0.5, orientation='y') + geom_quasirandom(data=filter(da_results, is.sig),
           aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_income_beeswarm_'%&%b%&%'b4.png', height=5, width=7)
    ## save results
    if (exists('compiled_milo')){
      compiled_milo <- rbind(compiled_milo, da_results)
    } else {compiled_milo <- da_results}
    
    # differential abundance in nhoods for asthma
    da_results <- testNhoods(subset_obj, design=~batch+age+gender+asthma, design.df=milo_design)
    da_results <- annotateNhoods(subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='asthma', condition=cond, batch4=b,
                                        is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + geom_quasirandom(data=filter(da_results, !is.sig),
           aes(x=logFC), color='grey', alpha=0.5, orientation='y') + geom_quasirandom(data=filter(da_results, is.sig),
           aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_asthma_beeswarm_'%&%b%&%'b4.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # perform test for infection and interactions
    if (cond!='NI'){
      # recreate objects
      if (b=='yes'){
        subset_obj <- subset(obj, subset= (condition==cond | condition=='NI')) %>% as.SingleCellExperiment() %>% Milo()
      } else {
        subset_obj <- subset(obj, subset= (condition==cond | condition=='NI') & batch!='B4') %>% as.SingleCellExperiment() %>% Milo()
      }
      
      # construct KNN graph and nhoods
      subset_obj <- buildGraph(subset_obj, k=50, d=30) %>%
        makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA')
      
      # count cells in nhoods per experimental condition (aka infection status)
      colData(subset_obj)$group <- colData(subset_obj)$IDs%&%'_'%&%colData(subset_obj)$condition
      subset_obj <- countCells(subset_obj, meta.data=data.frame(colData(subset_obj)), sample='group')
      
      # create design matrix
      milo_design <- data.frame(colData(subset_obj))[,c('IDs','group','batch','age','gender','condition','asthma','income')] %>% distinct()
      rownames(milo_design) <- milo_design$group
      milo_design$condition <- factor(milo_design$condition, levels=c('NI', cond))
      milo_design$asthma <- factor(milo_design$asthma, levels=c('No', 'Yes'))
      milo_design$income <- factor(milo_design$income, levels=c('Low', 'High'))
      
      # compute nhood distance
      subset_obj <- calcNhoodDistance(subset_obj, d=30, reduced.dim='INTEGRATED.RPCA') 
      
      # differential abundance in nhoods for infection
      da_results <- testNhoods(subset_obj, design=~batch+age+gender+condition, design.df=milo_design)
      da_results <- annotateNhoods(subset_obj, da_results, 'celltype')
      da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
      da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
      da_results <- da_results %>% mutate(identity='infection', condition=cond, batch4=b,
                                          is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
      ggplot(da_results, aes(y=celltype)) + geom_quasirandom(data=filter(da_results, !is.sig),
             aes(x=logFC), color='grey', alpha=0.5, orientation='y') + geom_quasirandom(data=filter(da_results, is.sig),
             aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
      ggsave('proportion_plots/miloR_'%&%cond%&%'_infection_beeswarm_'%&%b%&%'b4.png', height=5, width=7)
      ## save results
      compiled_milo <- rbind(compiled_milo, da_results)
      
      # differential abundance in nhoods for infection and income
      da_results <- testNhoods(subset_obj, design=~batch+age+gender+condition*income, design.df=milo_design)
      da_results <- annotateNhoods(subset_obj, da_results, 'celltype')
      da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
      da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
      da_results <- da_results %>% mutate(identity='infection:income', condition=cond, batch4=b,
                                          is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
      ggplot(da_results, aes(y=celltype)) + geom_quasirandom(data=filter(da_results, !is.sig),
             aes(x=logFC), color='grey', alpha=0.5, orientation='y') + geom_quasirandom(data=filter(da_results, is.sig),
             aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
      ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_income_beeswarm_'%&%b%&%'b4.png', height=5, width=7)
      ## save results
      compiled_milo <- rbind(compiled_milo, da_results)
      
      # differential abundance in nhoods for infection and asthma
      da_results <- testNhoods(subset_obj, design=~batch+age+gender+condition*asthma, design.df=milo_design)
      da_results <- annotateNhoods(subset_obj, da_results, 'celltype')
      da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
      da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
      da_results <- da_results %>% mutate(identity='infection:asthma', condition=cond, batch4=b,
                                          is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
      ggplot(da_results, aes(y=celltype)) + geom_quasirandom(data=filter(da_results, !is.sig),
             aes(x=logFC), color='grey', alpha=0.5, orientation='y') + geom_quasirandom(data=filter(da_results, is.sig),
             aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
      ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_asthma_beeswarm_'%&%b%&%'b4.png', height=5, width=7)
      ## save results
      compiled_milo <- rbind(compiled_milo, da_results)
    }
  }
}
rm(subset_obj, milo_design, da_results)
# look at results
compiled_milo$condition <- factor(compiled_milo$condition, levels=conditions)
compiled_milo$identity <- factor(compiled_milo$identity, levels=c('infection', 'asthma', 'income', 'infection:asthma', 'infection:income'))

# summarize over condition, celltype, identity, and batch4
summary_sig_compiled_milo <- compiled_milo %>%
  group_by(condition, celltype, identity, batch4) %>%
  summarise(
    n_nhoods = n(),
    n_sig = sum(is.sig),
    frac_sig = n_sig / n_nhoods,
    median_logFC_sig = median(logFC[is.sig], na.rm=TRUE),
    direction = sign(sum(logFC[is.sig], na.rm=TRUE)),
    .groups = 'drop') %>% drop_na() %>%
  mutate(direction=ifelse(direction==1, 'Up', 'Down'))

summary_sig_compiled_milo %>% filter(batch4=='yes') %>% 
  ggplot(., aes(x=celltype, y=identity, size=frac_sig, color=median_logFC_sig, shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='Fraction of sig. nhoods', color='Median LogFC', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/miloR_sigresults_yesb4.png', height=5, width=10)

summary_sig_compiled_milo %>% filter(batch4=='no') %>% 
  ggplot(., aes(x=celltype, y=identity, size=frac_sig, color=median_logFC_sig, shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='Fraction of sig. nhoods', color='Median LogFC', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/miloR_sigresults_nob4.png', height=5, width=10)
