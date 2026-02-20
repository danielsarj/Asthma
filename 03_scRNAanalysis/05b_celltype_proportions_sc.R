library(tidyverse)
library(Seurat)
library(data.table)
library(speckle)
library(miloR)
library(SingleCellExperiment)
library(ggbeeswarm)
"%&%" <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('NI', 'IVA', 'RV')
celltypes <- c('B','CD4-T','CD8-T','Mono','NK')
interactions <- c('none','asthma', 'income', 'ACT', 'ACE', 'resilience', 's_support', 'total_r',
          'year_r', 'life_r', 'stress_r', 'kid24h_r', 'kid_d', 'infection_c')

# load seurat object
obj <- readRDS('NI_IVA_RV.integrated.w_celltype_new.rds')

# load sample metadata (age, sex, asthma, income, albuterol) and add to seurat's metadata
sample_m <- fread('../sample_metadata.txt')
sample_m$income <- na_if(sample_m$income, '')
sample_m$income <- ifelse(is.na(sample_m$income), NA, 
              ifelse(sample_m$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'), 'Low', 'High'))
sample_m$albuterol <- na_if(sample_m$albuterol, '')
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='ID'))
rownames(mdata) <- rownames(obj@meta.data)
obj@meta.data <- mdata
rm(sample_m, mdata)

# load sample metadata (act, ace, resilience, social support, racism, albuterol) and add to seurat's metadata
sample_m <- fread('../Sample_test_scores_Araujo.tsv') %>% 
  mutate(ACE_result = coalesce(ACE_Percent_Self, ACE_Percent_Parent)) %>%
  select(Study_ID, Recorded_Diagnosis, ACT_score, ACE_result, Parent_Resilience_Score, Parents_Score_Avg, 
     Total_Racist_Events, Year_Racist_events, Life_Racist_events, Racist_stress, Racism_child_24hr, Experience_Discrimination_child)
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='Study_ID'))
rownames(mdata) <- rownames(obj@meta.data)
obj@meta.data <- mdata
rm(sample_m, mdata)

# load sample metadata (infection at time of collectiong) and add to seurat's metadata
sample_m <- fread('../SEA_Metadata_Pathogen_Araujo.tsv') %>% 
  dplyr::rename(infection_agent = Results, infection_status = Comment)
sample_m$infection_status <- gsub('infection', 'Positive', sample_m$infection_status)
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='Study.ID.'))
rownames(mdata) <- rownames(obj@meta.data)
obj@meta.data <- mdata
rm(sample_m, mdata)

# remove batch 4
obj <- subset(obj, subset= batch!='B4')

###############
### speckle ###
###############
for (cond in conditions){
  # create object
  subset_obj <- subset(obj, subset=condition==cond) 
  # transform proportions 
  props <- getTransformedProps(subset_obj@meta.data$celltype, subset_obj@meta.data$IDs, transform='logit')
  # subset metadata
  subset_mdata <- subset_obj@meta.data %>% filter(condition==cond, IDs %in% colnames(props[[1]])) %>%
  select(IDs, condition, batch, age, gender, Recorded_Diagnosis, income, ACT_score,
       ACE_result, Parent_Resilience_Score, Parents_Score_Avg, Total_Racist_Events,
       Year_Racist_events, Life_Racist_events, Racist_stress, Racism_child_24hr, 
       Experience_Discrimination_child, infection_status) %>% unique() %>% arrange(IDs)
  
  subset_mdata$gender <- factor(subset_mdata$gender, levels=c('Male','Female'))
  subset_mdata$batch <- factor(subset_mdata$batch, levels=c('B1','B2','B3'))
  subset_mdata$income <- factor(subset_mdata$income, levels=c('Low','High'))
  subset_mdata$Recorded_Diagnosis <- factor(subset_mdata$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
  subset_mdata$infection_status <- factor(subset_mdata$infection_status, levels=c('Negative', 'Positive'))

  # perform test for asthma
  ## create design matrix
  design_matrix <- model.matrix(~0+Recorded_Diagnosis+batch+age+gender, data=subset_mdata)
  mycontr <- makeContrasts(Recorded_DiagnosisRecorded_Asthma_Diagnosis-Recorded_DiagnosisNo_Diagnosis, levels=design_matrix)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Recorded_Diagnosis)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
    rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
    mutate(identity='asthma', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
  ## save results
  if (exists('compiled_speckle')){
    compiled_speckle <- rbind(compiled_speckle, results)
  } else {compiled_speckle <- results}
  
  # perform test for income
  ## create design matrix
  design_matrix <- model.matrix(~0+income+batch+age+gender, data=subset_mdata)
  mycontr <- makeContrasts(incomeHigh-incomeLow, levels=design_matrix)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(income)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
    rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
    mutate(identity='income', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for ACT
  ## create design matrix
  design_matrix <- model.matrix(~0+ACT_score+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(ACT_score)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='ACT', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for ACE
  ## create design matrix
  design_matrix <- model.matrix(~0+ACE_result+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(ACE_result)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='ACE', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for resilience
  ## create design matrix
  design_matrix <- model.matrix(~0+Parent_Resilience_Score+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Parent_Resilience_Score)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='resilience', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for social support
  ## create design matrix
  design_matrix <- model.matrix(~0+Parents_Score_Avg+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Parents_Score_Avg)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='s_support', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for total racism
  ## create design matrix
  design_matrix <- model.matrix(~0+Total_Racist_Events+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Total_Racist_Events)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='total_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for year racism
  ## create design matrix
  design_matrix <- model.matrix(~0+Year_Racist_events+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Year_Racist_events)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='year_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for life racism
  ## create design matrix
  design_matrix <- model.matrix(~0+Life_Racist_events+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Life_Racist_events)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='life_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for stress racism
  ## create design matrix
  design_matrix <- model.matrix(~0+Racist_stress+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Racist_stress)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='stress_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for kid24h racism
  ## create design matrix
  design_matrix <- model.matrix(~0+Racism_child_24hr+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Racism_child_24hr)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='kid24h_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for kid discrimination
  ## create design matrix
  design_matrix <- model.matrix(~0+Experience_Discrimination_child+batch+age+gender, data=subset_mdata)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(Experience_Discrimination_child)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
    topTable(coef=1, number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
    dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='kid_d', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for infection during collection
  ## create design matrix
  design_matrix <- model.matrix(~0+infection_status+batch+age+gender, data=subset_mdata)
  mycontr <- makeContrasts(infection_statusPositive-infection_statusNegative, levels=design_matrix)
  # subset props object
  idx <- subset_mdata %>% filter(!is.na(infection_status)) %>% select(IDs) %>% pull()
  sub_props <- props
  sub_props$Counts <- props$Counts[,idx, drop = FALSE]
  sub_props$TransformedProps <- props$TransformedProps[,idx, drop = FALSE]
  sub_props$Proportions <- props$Proportions[,idx, drop = FALSE]
  ## run propeller
  results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
    rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
    mutate(identity='infection_c', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
  ## save results
  compiled_speckle <- rbind(compiled_speckle, results)
  
  # perform test for infection and interactions
  if (cond!='NI'){
    # recreate objects
    subset_obj <- subset(obj, subset= (condition==cond | condition=='NI')) 
    
    # transform proportions 
    subset_obj_mdata <- subset_obj@meta.data %>% mutate(identifier=IDs%&%'_'%&%condition)
    props <- getTransformedProps(subset_obj_mdata$celltype, subset_obj_mdata$identifier, transform='logit')
    # subset metadata
    subset_mdata <- subset_obj_mdata %>% select(identifier, condition, batch, age, gender, Recorded_Diagnosis, income, ACT_score,
             ACE_result, Parent_Resilience_Score, Parents_Score_Avg, Total_Racist_Events,
             Year_Racist_events, Life_Racist_events, Racist_stress, Racism_child_24hr, 
             Experience_Discrimination_child, infection_status) %>% unique() %>% arrange(identifier)
    subset_mdata$condition <- factor(subset_mdata$condition, levels=c('NI', cond))
    subset_mdata$gender <- factor(subset_mdata$gender, levels=c('Male','Female'))
    subset_mdata$batch <- factor(subset_mdata$batch, levels=c('B1','B2','B3'))
    subset_mdata$income <- factor(subset_mdata$income, levels=c('Low','High'))
    subset_mdata$Recorded_Diagnosis <- factor(subset_mdata$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
    subset_mdata$infection_status <- factor(subset_mdata$infection_status, levels=c('Negative', 'Positive'))
    
    # create design matrix for infection
    design_matrix <- model.matrix(~0+condition+batch+age+gender, data=subset_mdata)
    if (cond=='IVA'){
    mycontr <- makeContrasts(conditionIVA-conditionNI, levels=design_matrix)
    # run propeller
    results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
      mutate(identity='infection', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    } else {
    mycontr <- makeContrasts(conditionRV-conditionNI, levels=design_matrix)
    # run propeller
    results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
      mutate(identity='infection', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    }
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for asthma interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Recorded_Diagnosis))
    sub_subset_mdata$group <- with(sub_subset_mdata, paste(condition, Recorded_Diagnosis, sep='_'))
    sub_subset_mdata$group <- factor(sub_subset_mdata$group, 
                                 levels = c('NI_No_Diagnosis', 'NI_Recorded_Asthma_Diagnosis', cond%&%'_No_Diagnosis', cond%&%'_Recorded_Asthma_Diagnosis'))
    design_matrix <- model.matrix(~0+group+batch+age+gender, data=sub_subset_mdata)
    if (cond=='IVA'){
      mycontr <- makeContrasts(interaction=(groupIVA_Recorded_Asthma_Diagnosis-groupIVA_No_Diagnosis)-(groupNI_Recorded_Asthma_Diagnosis-groupNI_No_Diagnosis), levels=design_matrix)
    # run propeller
    results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
      mutate(identity='infection:asthma', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    } else {
    mycontr <- makeContrasts(interaction=(groupRV_Recorded_Asthma_Diagnosis-groupRV_No_Diagnosis)-(groupNI_Recorded_Asthma_Diagnosis-groupNI_No_Diagnosis), levels=design_matrix)
    # run propeller
    results <- propeller.ttest(props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
      mutate(identity='infection:asthma', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    }
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for income interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(income))
    sub_subset_mdata$group <- with(sub_subset_mdata, paste(condition, income, sep='_'))
    sub_subset_mdata$group <- factor(sub_subset_mdata$group, levels = c('NI_Low', 'NI_High', cond%&%'_Low', cond%&%'_High'))
    design_matrix <- model.matrix(~0+group+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    if (cond=='IVA'){
      mycontr <- makeContrasts(interaction=(groupIVA_High-groupIVA_Low)-(groupNI_High-groupNI_Low), levels=design_matrix)
      # run propeller
      results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
        rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
        mutate(identity='infection:income', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    } else {
    mycontr <- makeContrasts(interaction=(groupRV_High-groupRV_Low)-(groupNI_High-groupNI_Low), levels=design_matrix)
    # run propeller
    results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
      rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
      mutate(identity='infection:income', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    }
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for ACT interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(ACT_score))
    design_matrix <- model.matrix(~0+ACT_score*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
     results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
       topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:ACT', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for ACE interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(ACE_result))
    design_matrix <- model.matrix(~0+ACE_result*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:ACE', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for resilience interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Parent_Resilience_Score))
    design_matrix <- model.matrix(~0+Parent_Resilience_Score*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:resilience', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for social support interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Parents_Score_Avg))
    design_matrix <- model.matrix(~0+Parents_Score_Avg*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:s_support', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for total racism interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Total_Racist_Events))
    design_matrix <- model.matrix(~0+Total_Racist_Events*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:total_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for year racism interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Year_Racist_events))
    design_matrix <- model.matrix(~0+Year_Racist_events*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:year_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for life racism interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Life_Racist_events))
    design_matrix <- model.matrix(~0+Life_Racist_events*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:life_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for stress racism interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Racist_stress))
    design_matrix <- model.matrix(~0+Racist_stress*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:stress_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for kid24h racism interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Racism_child_24hr))
    design_matrix <- model.matrix(~0+Racism_child_24hr*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:kid24h_r', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for kid discrimination interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(Experience_Discrimination_child))
    design_matrix <- model.matrix(~0+Experience_Discrimination_child*condition+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    # run propeller
    results <- lmFit(sub_props$Proportions, design_matrix) %>% eBayes(robust=TRUE) %>%
      topTable(coef=ncol(.), number=Inf) %>% rownames_to_column(var='clusters') %>% select(clusters, logFC, adj.P.Val) %>%
      dplyr::rename(FDR=adj.P.Val, change=logFC) %>% mutate(identity='infection:kid_d', condition=cond, direction=ifelse(change>0, 'Up', 'Down'))
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
    
    # create design matrix for infection at collection interaction
    sub_subset_mdata <- subset_mdata %>% filter(!is.na(infection_status))
    sub_subset_mdata$group <- with(sub_subset_mdata, paste(condition, infection_status, sep='_'))
    sub_subset_mdata$group <- factor(sub_subset_mdata$group, levels = c('NI_Negative', 'NI_Positive', cond%&%'_Negative', cond%&%'_Positive'))
    design_matrix <- model.matrix(~0+group+batch+age+gender, data=sub_subset_mdata)
    # subset props object
    sub_props <- props
    sub_props$Counts <- props$Counts[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$TransformedProps <- props$TransformedProps[,sub_subset_mdata$identifier, drop = FALSE]
    sub_props$Proportions <- props$Proportions[,sub_subset_mdata$identifier, drop = FALSE]
    if (cond=='IVA'){
      mycontr <- makeContrasts(interaction=(groupIVA_Positive-groupIVA_Negative)-(groupNI_Positive-groupNI_Negative), levels=design_matrix)
      # run propeller
      results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
        rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
        mutate(identity='infection:infection_c', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    } else {
      mycontr <- makeContrasts(interaction=(groupRV_Positive-groupRV_Negative)-(groupNI_Positive-groupNI_Negative), levels=design_matrix)
      # run propeller
      results <- propeller.ttest(sub_props, design_matrix, mycontr, robust=TRUE, trend=FALSE, sort=FALSE) %>%
        rownames_to_column(var='clusters') %>% select(clusters, PropRatio, FDR) %>% dplyr::rename(change=PropRatio) %>% 
        mutate(identity='infection:infection_c', condition=cond, direction=ifelse(change>1, 'Up', 'Down'))
    }
    # save results
    compiled_speckle <- rbind(compiled_speckle, results)
  }
}
rm(subset_obj, results, props, sub_props, subset_mdata, sub_subset_mdata, subset_obj_mdata, mycontr, design_matrix)
fwrite(compiled_speckle, 'compiled_speckle_results.txt', sep=' ')

# look at results
compiled_speckle$condition <- factor(compiled_speckle$condition, levels=conditions)
mod_interactions <- gsub('none', 'infection', interactions)
prefixed_interactions <- paste0("infection:", interactions[-1])
compiled_speckle$identity <- factor(compiled_speckle$identity, levels=c(mod_interactions, prefixed_interactions))

sig_compiled_speckle <- compiled_speckle %>% filter(FDR<0.05) 

ggplot(sig_compiled_speckle, aes(x=clusters, y=identity, size=-log10(FDR), shape=direction)) +
geom_point() + scale_color_gradient(low='blue', high='red') +  
labs(x='Cell type', y=NULL, size='-log10(FDR)', shape='Direction') +
theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/speckle_sigresults_new.png', height=4, width=6)

#############
### miloR ###
#############
for (cond in conditions){
  subset_obj <- subset(obj, subset=condition==cond) %>% as.SingleCellExperiment() %>% Milo()
  # create design matrix
  milo_design <- data.frame(colData(subset_obj))[,c('IDs','batch','age','gender','income','Recorded_Diagnosis',
                                                    'ACT_score','ACE_result','Parent_Resilience_Score','Parents_Score_Avg', 
                                                    'Total_Racist_Events','Year_Racist_events','Life_Racist_events', 
                                                    'Racist_stress','Racism_child_24hr','Experience_Discrimination_child', 
                                                    'infection_status')] %>% distinct()
  rownames(milo_design) <- milo_design$IDs
  milo_design$gender <- factor(milo_design$gender, levels=c('Male','Female'))
  milo_design$batch <- factor(milo_design$batch, levels=c('B1','B2','B3'))
  milo_design$income <- factor(milo_design$income, levels=c('Low','High'))
  milo_design$Recorded_Diagnosis <- factor(milo_design$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
  milo_design$infection_status <- factor(milo_design$infection_status, levels=c('Negative', 'Positive'))

  # income
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(income))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+income, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='income', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_income_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_income_nhoodgraph_new.png', height=5, width=7)
  ## save results
  if (exists('compiled_milo')){
    compiled_milo <- rbind(compiled_milo, da_results)
  } else {compiled_milo <- da_results}
  
  # asthma
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Recorded_Diagnosis))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Recorded_Diagnosis, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='asthma', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_asthma_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_asthma_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # ACT
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(ACT_score))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+ACT_score, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='ACT', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_ACT_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_ACT_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # ACE
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(ACE_result))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+ACE_result, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='ACE', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_ACE_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_ACE_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # resilience
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Parent_Resilience_Score))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Parent_Resilience_Score, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='resilience', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_resilience_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_resilience_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # social support
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Parents_Score_Avg))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Parents_Score_Avg, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='s_support', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_social_support_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_social_support_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # total racism
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Total_Racist_Events))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Total_Racist_Events, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='total_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_total_racism_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_total_racism_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # year racism
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Year_Racist_events))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Year_Racist_events, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='year_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_year_racism_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_year_racism_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # life racism
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Life_Racist_events))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Life_Racist_events, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='life_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_life_racism_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_life_racism_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # stress racism
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Racist_stress))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Racist_stress, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='stress_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_stress_racism_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_stress_racism_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # kid24h racism
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Racism_child_24hr))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Racism_child_24hr, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='kid24h_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_kid24h_racism_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_kid24h_racism_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # kid discrimination
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(Experience_Discrimination_child))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+Experience_Discrimination_child, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='kid_d', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_kid_discrimination_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_kid_discrimination_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # infection at collection
  # create design matrix 
  sub_milo_design <- milo_design %>% filter(!is.na(infection_status))
  no_NA_idx <- as.vector(sub_milo_design$IDs)
  # subset milo object
  sub_subset_obj <- subset_obj[,colData(subset_obj)$IDs %in% no_NA_idx]
  # construct KNN graph and nhoods
  sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
    makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
    buildNhoodGraph()
  reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
  # count cells in nhoods
  sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='IDs')
  # compute nhood distance
  sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
  # differential abundance in nhoods 
  da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+infection_status, design.df=sub_milo_design)
  da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
  da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
  da_results <- da_results %>% mutate(identity='infection_c', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
  da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
  ggplot(da_results, aes(y=celltype)) + 
    geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
    geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
  ggsave('proportion_plots/miloR_'%&%cond%&%'_infection_collection_beeswarm_new.png', height=5, width=7)
  plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
  ggsave('proportion_plots/miloR_'%&%cond%&%'_infection_collection_nhoodgraph_new.png', height=5, width=7)
  ## save results
  compiled_milo <- rbind(compiled_milo, da_results)
  
  # perform test for infection and interactions
  if (cond!='NI'){
    # recreate objects
    subset_obj <- subset(obj, subset= (condition==cond | condition=='NI')) %>% as.SingleCellExperiment() %>% Milo()
    colData(subset_obj)$group <- colData(subset_obj)$IDs%&%'_'%&%colData(subset_obj)$condition
    
    # create design matrix
    milo_design <- data.frame(colData(subset_obj))[,c('group','IDs','batch','age','gender','condition','income','Recorded_Diagnosis',
                                                      'ACT_score','ACE_result','Parent_Resilience_Score','Parents_Score_Avg', 
                                                      'Total_Racist_Events','Year_Racist_events','Life_Racist_events', 
                                                      'Racist_stress','Racism_child_24hr','Experience_Discrimination_child', 
                                                      'infection_status')] %>% distinct()
    rownames(milo_design) <- milo_design$group
    milo_design$condition <- factor(milo_design$condition, levels=c('NI', cond))
    milo_design$gender <- factor(milo_design$gender, levels=c('Male','Female'))
    milo_design$batch <- factor(milo_design$batch, levels=c('B1','B2','B3'))
    milo_design$income <- factor(milo_design$income, levels=c('Low','High'))
    milo_design$Recorded_Diagnosis <- factor(milo_design$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
    milo_design$infection_status <- factor(milo_design$infection_status, levels=c('Negative', 'Positive'))
    
    # infection
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(condition))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_infection_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_infection_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:income
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(income))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*income, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:income', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_income_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_income_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:asthma
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Recorded_Diagnosis))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Recorded_Diagnosis, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:asthma', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_asthma_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_asthma_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:ACT
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(ACT_score))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*ACT_score, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:ACT', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_ACT_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_ACT_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:ACE
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(ACE_result))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*ACE_result, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:ACE', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_ACE_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_ACE_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:resilience
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Parent_Resilience_Score))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Parent_Resilience_Score, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:resilience', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_resilience_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_resilience_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:social_support
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Parents_Score_Avg))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Parents_Score_Avg, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:s_support', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_social_support_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_social_support_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:total_racism
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Total_Racist_Events))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Total_Racist_Events, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:total_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_total_racism_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_total_racism_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:year_racism
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Year_Racist_events))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Year_Racist_events, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:year_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_year_racism_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_year_racism_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:life_racism
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Life_Racist_events))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Life_Racist_events, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:life_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_life_racism_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_life_racism_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:stress_racism
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Racist_stress))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Racist_stress, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:stress_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_stress_racism_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_stress_racism_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:kid24h_racism
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Racism_child_24hr))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Racism_child_24hr, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:kid24h_r', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_kid24h_racism_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_kid24h_racism_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:kid_discrimination
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(Experience_Discrimination_child))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*Experience_Discrimination_child, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:kid_d', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_kid_discrimination_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_kid_discrimination_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
    
    # infection:infection_collection
    # create design matrix 
    sub_milo_design <- milo_design %>% filter(!is.na(infection_status))
    no_NA_idx <- as.vector(sub_milo_design$group)
    # subset milo object
    sub_subset_obj <- subset_obj[,colData(subset_obj)$group %in% no_NA_idx]
    # construct KNN graph and nhoods
    sub_subset_obj <- buildGraph(sub_subset_obj, k=50, d=30) %>%
      makeNhoods(k=50, d=30, refined=TRUE, prop=0.1, reduced_dims='INTEGRATED.RPCA') %>%
      buildNhoodGraph()
    reducedDim(sub_subset_obj, 'UMAP') <- reducedDim(sub_subset_obj, 'RNA.UMAP')
    # count cells in nhoods
    sub_subset_obj <- countCells(sub_subset_obj, meta.data=data.frame(colData(sub_subset_obj)), sample='group')
    # compute nhood distance
    sub_subset_obj <- calcNhoodDistance(sub_subset_obj, d=30, reduced.dim='INTEGRATED.RPCA')   
    # differential abundance in nhoods 
    da_results <- testNhoods(sub_subset_obj, design=~batch+age+gender+condition*infection_status, design.df=sub_milo_design)
    da_results <- annotateNhoods(sub_subset_obj, da_results, 'celltype')
    da_results$celltype <- ifelse(da_results$celltype_fraction<0.7, 'Mixed', da_results$celltype)
    da_results <- da_results %>% mutate(identity='infection:infection_c', condition=cond, is.sig=ifelse(SpatialFDR<0.05, TRUE, FALSE))
    da_results$celltype <- factor(da_results$celltype, levels=c(celltypes, 'Mixed'))
    ggplot(da_results, aes(y=celltype)) + 
      geom_quasirandom(data=filter(da_results, !is.sig), aes(x=logFC), color='grey', alpha=0.5, orientation='y') + 
      geom_quasirandom(data=filter(da_results, is.sig), aes(x=logFC), color='red', alpha=1, orientation='y') + theme_bw()
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_infection_collection_beeswarm_new.png', height=5, width=7)
    plotNhoodGraphDA(sub_subset_obj, da_results, alpha=0.05)
    ggsave('proportion_plots/miloR_'%&%cond%&%'_interaction_infection_collection_nhoodgraph_new.png', height=5, width=7)
    ## save results
    compiled_milo <- rbind(compiled_milo, da_results)
  }
}
rm(subset_obj, milo_design, da_results)
fwrite(compiled_milo, 'compiled_miloR_results.txt', sep=' ')

# look at results
compiled_milo$condition <- factor(compiled_milo$condition, levels=conditions)
compiled_milo$identity <- factor(compiled_milo$identity, levels=c(mod_interactions, prefixed_interactions))

# summarize over condition, celltype, identity, and batch4
summary_sig_compiled_milo <- compiled_milo %>%
  group_by(condition, celltype, identity) %>%
  summarise(
    n_nhoods = n(),
    n_sig = sum(is.sig),
    frac_sig = n_sig / n_nhoods,
    median_logFC_sig = median(logFC[is.sig], na.rm=TRUE),
    direction = sign(sum(logFC[is.sig], na.rm=TRUE)),
    .groups = 'drop') %>% drop_na() %>%
  mutate(direction=ifelse(direction==1, 'Up', 'Down')) %>%
  filter(celltype!='Mixed')

ggplot(summary_sig_compiled_milo, aes(x=celltype, y=identity, size=frac_sig, color=median_logFC_sig, shape=direction)) +
  geom_point() + scale_color_gradient(low='blue', high='red') +  
  labs(x='Cell type', y=NULL, size='Fraction of sig. nbhds', color='Median LogFC', shape='Direction') +
  theme_bw() + facet_wrap(~condition)
ggsave('proportion_plots/miloR_sigresults_new.png', height=5, width=12)
