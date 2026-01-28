library(Seurat)
library(tidyverse)
library(edgeR)
library(limma)
library(broom)
library(data.table)
library(betareg)
library(DirichletReg)
library(ggpubr)
"%&%" <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B','CD4-T','CD8-T','Mono','NK')

# function to compute RMSE
rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

# load seurat object
obj <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')

# low level comparisons
mdata <- obj@meta.data %>% select(IDs, batch, condition, celltype, prop, asthma, income) 
sub_iva_mdata <- mdata %>% filter(condition!='RV') %>% mutate(infection='IVA') 
sub_iva_mdata$condition <- gsub('IVA', 'Inf', sub_iva_mdata$condition)
sub_rv_mdata <- mdata %>% filter(condition!='IVA') %>% mutate(infection='RV') 
sub_rv_mdata$condition <- gsub('RV', 'Inf', sub_rv_mdata$condition)
joint_df <- rbind(sub_iva_mdata, sub_rv_mdata)
rm(sub_iva_mdata, sub_rv_mdata)
joint_df$condition <- factor(joint_df$condition, levels=c('NI', 'Inf'))

# proportion of each cell type per NI or Inf group
avg_df <- joint_df %>% group_by(condition, celltype, infection) %>%
  summarise(prop = mean(prop), .groups = 'drop')
ggplot(joint_df, aes(x=condition, y=prop, group=IDs)) + geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=1, color='red'), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop), color='red', size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() + theme(legend.position='none')
ggsave('proportion_plots/celltype_prop_simple.png', height=5, width=5)
#ggsave('proportion_plots/celltype_prop_simple.pdf', height=5, width=5)

## same thing, but without batch 4
avg_df <- joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection) %>%
  summarise(prop = mean(prop), .groups = 'drop')
joint_df %>% filter(batch!='B4') %>% ggplot(., aes(x=condition, y=prop, group=IDs)) + 
  geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=1, color='red'), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop), color='red', size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() + theme(legend.position='none')
ggsave('proportion_plots/celltype_prop_simple_noB4.png', height=5, width=5)
#ggsave('proportion_plots/celltype_prop_simple_noB4.pdf', height=5, width=5)

# proportion of each cell type per NI or Inf group per asthma status
avg_df <- joint_df %>% group_by(condition, celltype, infection, asthma) %>%
  summarise(prop = mean(prop), .groups = 'drop')
ggplot(joint_df, aes(x=condition, y=prop, group=IDs, color=asthma)) + geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=asthma, color=asthma), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop, color=asthma), size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() 
ggsave('proportion_plots/celltype_prop_simple_asthma.png', height=5, width=5)
#ggsave('proportion_plots/celltype_prop_simple_asthma.pdf', height=5, width=5)

# test if there are differences in proportion
joint_df %>% group_by(condition, celltype, infection) %>% filter(n_distinct(asthma) == 2) %>% 
  do(tidy(t.test(prop ~ asthma, data = .))) %>% ungroup() %>% mutate(p_adj=p.adjust(p.value, method='fdr')) %>%
  filter(p.value<0.05)

## same thing, but without batch 4
# proportion of each cell type per NI or Inf group per asthma status
avg_df <- joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection, asthma) %>%
  summarise(prop = mean(prop), .groups = 'drop')
joint_df %>% filter(batch!='B4') %>% ggplot(., aes(x=condition, y=prop, group=IDs, color=asthma)) + 
  geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=asthma, color=asthma), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop, color=asthma), size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() 
ggsave('proportion_plots/celltype_prop_simple_asthma_noB4.png', height=5, width=5)
#ggsave('proportion_plots/celltype_prop_simple_asthma_noB4.pdf', height=5, width=5)

# test if there are differences in proportion
joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection) %>% filter(n_distinct(asthma)==2) %>% 
  do(tidy(t.test(prop ~ asthma, data = .))) %>% ungroup() %>% mutate(p_adj=p.adjust(p.value, method='fdr')) %>%
  filter(p.value<0.05)

# proportion of each cell type per NI or Inf group per income status
joint_df$income <- ifelse(joint_df$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'), 'Low', 'High')
joint_df$income <- factor(joint_df$income, levels=c('Low','High'))

avg_df <- joint_df %>% group_by(condition, celltype, infection, income) %>%
  summarise(prop = mean(prop), .groups = 'drop')
ggplot(joint_df, aes(x=condition, y=prop, group=IDs, color=income)) + geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=income, color=income), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop, color=income), size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() 
ggsave('proportion_plots/celltype_prop_simple_income.png', height=5, width=5)
#ggsave('proportion_plots/celltype_prop_simple_income.pdf', height=5, width=5)

# test if there are differences in proportion
joint_df %>% group_by(condition, celltype, infection) %>% filter(n_distinct(income) == 2) %>% 
  do(tidy(t.test(prop ~ income, data = .))) %>% ungroup() %>% mutate(p_adj=p.adjust(p.value, method='fdr')) %>%
  filter(p.value<0.05)

## same thing, but without batch 4
# proportion of each cell type per NI or Inf group per asthma status
avg_df <- joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection, income) %>%
  summarise(prop = mean(prop), .groups = 'drop')
joint_df %>% filter(batch!='B4') %>% ggplot(., aes(x=condition, y=prop, group=IDs, color=income)) + 
  geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=income, color=income), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop, color=income), size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() 
ggsave('proportion_plots/celltype_prop_simple_income_noB4.png', height=5, width=5)
#ggsave('proportion_plots/celltype_prop_simple_income_noB4.pdf', height=5, width=5)

# test if there are differences in proportion
joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection) %>% filter(n_distinct(income)==2) %>% 
  do(tidy(t.test(prop ~ income, data = .))) %>% ungroup() %>% mutate(p_adj=p.adjust(p.value, method='fdr')) %>%
  filter(p.value<0.05)

########################################
### MULTIPLE LINEAR REGRESSION MODEL ###
########################################
# first, test asthma and income per infection
for (i in 1:length(conditions)){
  for (ctype in celltypes){

    # extract and format metadata
    mdata <- obj@meta.data %>% filter(celltype==ctype, (condition==conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                       'Low', 'High')
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    mdata_nob4 <- mdata %>% filter(batch!='B4')
    
    # fit linear model in income, with batch 4 
    linear_m <- lm(prop~batch+age+gender+income, data=mdata) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='incomeHigh')
    if (exists('compiled_linear_income')){
      compiled_linear_income <- rbind(compiled_linear_income, results_linear_m)
    } else {compiled_linear_income <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+income, data=mdata)) 
    rmse_lin <- rmse(mdata$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_income')){
      rmse_compiled_linear_income <- rbind(rmse_compiled_linear_income, rmse_lin)
    } else {rmse_compiled_linear_income <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(celltype, income, residuals)
    ggplot(linear_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_income.png', height=3, width=4)
    
    # fit linear model in income, without batch 4 
    linear_m <- lm(prop~batch+age+gender+income, data=mdata_nob4) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='incomeHigh')
    if (exists('compiled_linear_income_nob4')){
      compiled_linear_income_nob4 <- rbind(compiled_linear_income_nob4, results_linear_m)
    } else {compiled_linear_income_nob4 <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+income, data=mdata_nob4)) 
    rmse_lin <- rmse(mdata_nob4$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_income_nob4')){
      rmse_compiled_linear_income_nob4 <- rbind(rmse_compiled_linear_income_nob4, rmse_lin)
    } else {rmse_compiled_linear_income_nob4 <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata_nob4) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(celltype, income, residuals)
    ggplot(linear_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_income_nob4.png', height=3, width=4)
    
    # fit linear model in asthma, with batch 4 
    linear_m <- lm(prop~batch+age+gender+asthma, data=mdata) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='asthmaYes')
    if (exists('compiled_linear_asthma')){
      compiled_linear_asthma <- rbind(compiled_linear_asthma, results_linear_m)
    } else {compiled_linear_asthma <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+asthma, data=mdata)) 
    rmse_lin <- rmse(mdata$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_asthma')){
      rmse_compiled_linear_asthma <- rbind(rmse_compiled_linear_asthma, rmse_lin)
    } else {rmse_compiled_linear_asthma <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(celltype, asthma, residuals)
    ggplot(linear_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_asthma.png', height=3, width=4)
    
    # fit linear model in asthma, without batch 4 
    linear_m <- lm(prop~batch+age+gender+asthma, data=mdata_nob4) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>%
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='asthmaYes')
    if (exists('compiled_linear_asthma_nob4')){
      compiled_linear_asthma_nob4 <- rbind(compiled_linear_asthma_nob4, results_linear_m)
    } else {compiled_linear_asthma_nob4 <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+asthma, data=mdata_nob4)) 
    rmse_lin <- rmse(mdata_nob4$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_asthma_nob4')){
      rmse_compiled_linear_asthma_nob4 <- rbind(rmse_compiled_linear_asthma_nob4, rmse_lin)
    } else {rmse_compiled_linear_asthma_nob4 <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata_nob4) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(celltype, asthma, residuals)
    ggplot(linear_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_asthma_nob4.png', height=3, width=4)
  }
}

########################################
### MULTIPLE LINEAR REGRESSION MODEL ###
########################################
# now, test asthma and income interacting with infection
conditions <- c('RV', 'IVA')
for (i in 1:length(conditions)){
  for (ctype in celltypes){
    
    # extract and format metadata
    mdata <- obj@meta.data %>% filter(celltype==ctype, (condition==conditions[i]|condition=='NI'))
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Low', 'High')
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    mdata <- mdata %>% filter(IDs %in% IDs[duplicated(IDs)])    # remove IDs that are not paired
    mdata_nob4 <- mdata %>% filter(batch!='B4')
    
    # fit linear model in infection, with batch 4 
    linear_m <- lm(prop~batch+age+gender+condition, data=mdata) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='condition'%&%conditions[i]) %>% mutate(terms='infection')
    if (exists('compiled_linear_infection')){
      compiled_linear_infection <- rbind(compiled_linear_infection, results_linear_m)
    } else {compiled_linear_infection <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+condition, data=mdata)) 
    rmse_lin <- rmse(mdata$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='infection') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_infection')){
      rmse_compiled_linear_infection <- rbind(rmse_compiled_linear_infection, rmse_lin)
    } else {rmse_compiled_linear_infection <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(celltype, condition, residuals)
    ggplot(linear_m, aes(x=condition, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('NI', conditions[i]))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_infection.png', height=3, width=4)
    
    # fit linear model in infection, without batch 4 
    linear_m <- lm(prop~batch+age+gender+condition, data=mdata_nob4) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='condition'%&%conditions[i]) %>% mutate(terms='infection')
    if (exists('compiled_linear_infection_nob4')){
      compiled_linear_infection_nob4 <- rbind(compiled_linear_infection_nob4, results_linear_m)
    } else {compiled_linear_infection_nob4 <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+condition, data=mdata_nob4)) 
    rmse_lin <- rmse(mdata_nob4$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='infection') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_infection_nob4')){
      rmse_compiled_linear_infection_nob4 <- rbind(rmse_compiled_linear_infection_nob4, rmse_lin)
    } else {rmse_compiled_linear_infection_nob4 <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata_nob4) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(celltype, condition, residuals)
    ggplot(linear_m, aes(x=condition, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('NI', conditions[i]))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_infection_nob4.png', height=3, width=4)
    
    # fit linear model in income, with batch 4 
    linear_m <- lm(prop~batch+age+gender+condition*income, data=mdata) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='condition'%&%conditions[i]%&%':incomeHigh')
    if (exists('compiled_linear_interaction_income')){
      compiled_linear_interaction_income <- rbind(compiled_linear_interaction_income, results_linear_m)
    } else {compiled_linear_interaction_income <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+condition*income, data=mdata)) 
    rmse_lin <- rmse(mdata$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='condition'%&%conditions[i]%&%':incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_interaction_income')){
      rmse_compiled_linear_interaction_income <- rbind(rmse_compiled_linear_interaction_income, rmse_lin)
    } else {rmse_compiled_linear_interaction_income <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(IDs, celltype, condition, income, residuals) %>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, income, residuals)
    ggplot(linear_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_infection.income.png', height=3, width=4)
    
    # fit linear model in income, without batch 4 
    linear_m <- lm(prop~batch+age+gender+condition*income, data=mdata_nob4) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='condition'%&%conditions[i]%&%':incomeHigh')
    if (exists('compiled_linear_interaction_income_nob4')){
      compiled_linear_interaction_income_nob4 <- rbind(compiled_linear_interaction_income_nob4, results_linear_m)
    } else {compiled_linear_interaction_income_nob4 <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+condition*income, data=mdata_nob4)) 
    rmse_lin <- rmse(mdata_nob4$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='condition'%&%conditions[i]%&%':incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_interaction_income_nob4')){
      rmse_compiled_linear_interaction_income_nob4 <- rbind(rmse_compiled_linear_interaction_income_nob4, rmse_lin)
    } else {rmse_compiled_linear_interaction_income_nob4 <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata_nob4) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(IDs, celltype, condition, income, residuals) %>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, income, residuals)
    ggplot(linear_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_infection.income_nob4.png', height=3, width=4)
    
    # fit linear model in asthma, with batch 4 
    linear_m <- lm(prop~batch+age+gender+condition*asthma, data=mdata)%>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='condition'%&%conditions[i]%&%':asthmaYes')
    if (exists('compiled_linear_interaction_asthma')){
      compiled_linear_interaction_asthma <- rbind(compiled_linear_interaction_asthma, results_linear_m)
    } else {compiled_linear_interaction_asthma <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+condition*asthma, data=mdata)) 
    rmse_lin <- rmse(mdata$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='condition'%&%conditions[i]%&%':asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_interaction_asthma')){
      rmse_compiled_linear_interaction_asthma <- rbind(rmse_compiled_linear_interaction_asthma, rmse_lin)
    } else {rmse_compiled_linear_interaction_asthma <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(IDs, celltype, condition, asthma, residuals) %>%
    pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, asthma, residuals)
    ggplot(linear_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_infection.asthma.png', height=3, width=4)
    
    # fit linear model in asthma, without batch 4 
    linear_m <- lm(prop~batch+age+gender+condition*asthma, data=mdata_nob4) %>% summary()
    betas <- linear_m$coefficients[,1]
    pvals <- linear_m$coefficients[,4]
    terms <- rownames(linear_m$coefficients)
    results_linear_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='condition'%&%conditions[i]%&%':asthmaYes')
    if (exists('compiled_linear_interaction_asthma_nob4')){
      compiled_linear_interaction_asthma_nob4 <- rbind(compiled_linear_interaction_asthma_nob4, results_linear_m)
    } else {compiled_linear_interaction_asthma_nob4 <- results_linear_m}
    # compute RMSE
    linear_m <- predict(lm(prop~batch+age+gender+condition*asthma, data=mdata_nob4)) 
    rmse_lin <- rmse(mdata_nob4$prop, linear_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='condition'%&%conditions[i]%&%':asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_linear_interaction_asthma_nob4')){
      rmse_compiled_linear_interaction_asthma_nob4 <- rbind(rmse_compiled_linear_interaction_asthma_nob4, rmse_lin)
    } else {rmse_compiled_linear_interaction_asthma_nob4 <- rmse_lin}
    ## plot residuals
    linear_m <- lm(prop~batch+age+gender, data=mdata_nob4) %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(IDs, celltype, condition, asthma, residuals) %>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, asthma, residuals)
    ggplot(linear_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_linear_infection.asthma_nob4.png', height=3, width=4)
  }
}

######################################
### MULTIPLE BETA REGRESSION MODEL ###
######################################
# first, test asthma and income per infection
conditions <- c('NI', 'RV', 'IVA')
for (i in 1:length(conditions)){
  for (ctype in celltypes){
    
    # extract and format metadata
    mdata <- obj@meta.data %>% filter(celltype==ctype, (condition==conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Low', 'High')
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    mdata_nob4 <- mdata %>% filter(batch!='B4')
    
    # fit beta model in income, with batch 4 
    beta_m <- betareg(prop~batch+age+gender+income, data=mdata, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='incomeHigh')
    if (exists('compiled_beta_income')){
      compiled_beta_income <- rbind(compiled_beta_income, results_beta_m)
    } else {compiled_beta_income <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+income, data=mdata, link='logit'), type='response')
    rmse_beta <- rmse(mdata$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_income')){
      rmse_compiled_beta_income <- rbind(rmse_compiled_beta_income, rmse_beta)
    } else {rmse_compiled_beta_income <- rmse_beta}
    ## plot residuals
    beta_m <- betareg(prop~batch+age+gender, data=mdata, link='logit')%>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(celltype, income, residuals)
    ggplot(beta_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_income.png', height=3, width=4)
    
    # fit beta model in income, without batch 4 
    beta_m <- betareg(prop~batch+age+gender+income, data=mdata_nob4, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='incomeHigh')
    if (exists('compiled_beta_income_nob4')){
      compiled_beta_income_nob4 <- rbind(compiled_beta_income_nob4, results_beta_m)
    } else {compiled_beta_income_nob4 <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+income, data=mdata_nob4, link='logit'), type='response')
    rmse_beta <- rmse(mdata_nob4$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_income_nob4')){
      rmse_compiled_beta_income_nob4 <- rbind(rmse_compiled_beta_income_nob4, rmse_beta)
    } else {rmse_compiled_beta_income_nob4 <- rmse_beta}
    ## plot residuals    
    beta_m <- betareg(prop~batch+age+gender, data=mdata_nob4, link='logit')%>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(celltype, income, residuals)
    ggplot(beta_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_income_nob4.png', height=3, width=4)
    
    # fit beta model in asthma, with batch 4 
    beta_m <- betareg(prop~batch+age+gender+asthma, data=mdata, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='asthmaYes')
    if (exists('compiled_beta_asthma')){
      compiled_beta_asthma <- rbind(compiled_beta_asthma, results_beta_m)
    } else {compiled_beta_asthma <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+asthma, data=mdata, link='logit'), type='response')
    rmse_beta <- rmse(mdata$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_asthma')){
      rmse_compiled_beta_asthma <- rbind(rmse_compiled_beta_asthma, rmse_beta)
    } else {rmse_compiled_beta_asthma <- rmse_beta}
    ## plot residuals    
    beta_m <- betareg(prop~batch+age+gender, data=mdata, link='logit')%>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(celltype, asthma, residuals)
    ggplot(beta_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_asthma.png', height=3, width=4)
    
    # fit beta model in asthma, without batch 4 
    beta_m <- betareg(prop~batch+age+gender+asthma, data=mdata_nob4, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='asthmaYes')
    if (exists('compiled_beta_asthma_nob4')){
      compiled_beta_asthma_nob4 <- rbind(compiled_beta_asthma_nob4, results_beta_m)
    } else {compiled_beta_asthma_nob4 <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+asthma, data=mdata_nob4, link='logit'), type='response')
    rmse_beta <- rmse(mdata_nob4$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_asthma_nob4')){
      rmse_compiled_beta_asthma_nob4 <- rbind(rmse_compiled_beta_asthma_nob4, rmse_beta)
    } else {rmse_compiled_beta_asthma_nob4 <- rmse_beta}
    ## plot residuals    
    beta_m <- betareg(prop~batch+age+gender, data=mdata_nob4, link='logit')%>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(celltype, asthma, residuals)
    ggplot(beta_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_asthma_nob4.png', height=3, width=4)
  }
}

######################################
### MULTIPLE BETA REGRESSION MODEL ###
######################################
# now, test asthma and income interacting with infection
conditions <- c('RV', 'IVA')
for (i in 1:length(conditions)){
  for (ctype in celltypes){
    
    # extract and format metadata
    mdata <- obj@meta.data %>% filter(celltype==ctype, (condition==conditions[i]|condition=='NI'))
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Low', 'High')
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    mdata <- mdata %>% filter(IDs %in% IDs[duplicated(IDs)])    # remove IDs that are not paired
    mdata_nob4 <- mdata %>% filter(batch!='B4')
    
    # fit beta model in infection, with batch 4 
    beta_m <- betareg(prop~batch+age+gender+condition, data=mdata, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='condition'%&%conditions[i]) %>% mutate(terms='infection')
    if (exists('compiled_beta_infection')){
      compiled_beta_infection <- rbind(compiled_beta_infection, results_beta_m)
    } else {compiled_beta_infection <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+condition, data=mdata, link='logit'), type='response')
    rmse_beta <- rmse(mdata$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='infection') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_infection')){
      rmse_compiled_beta_infection <- rbind(rmse_compiled_beta_infection, rmse_beta)
    } else {rmse_compiled_beta_infection <- rmse_beta}
    ## plot residuals    
    beta_m <- betareg(prop~batch+age+gender, data=mdata, link='logit') %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(celltype, condition, residuals)
    ggplot(beta_m, aes(x=condition, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('NI', conditions[i]))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_infection.png', height=3, width=4)
    
    # fit beta model in infection, without batch 4 
    beta_m <- betareg(prop~batch+age+gender+condition, data=mdata_nob4, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='condition'%&%conditions[i]) %>% mutate(terms='infection')
    if (exists('compiled_beta_infection_nob4')){
      compiled_beta_infection_nob4 <- rbind(compiled_beta_infection_nob4, results_beta_m)
    } else {compiled_beta_infection_nob4 <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+condition, data=mdata_nob4, link='logit'), type='response')
    rmse_beta <- rmse(mdata_nob4$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='infection') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_infection_nob4')){
      rmse_compiled_beta_infection_nob4 <- rbind(rmse_compiled_beta_infection_nob4, rmse_beta)
    } else {rmse_compiled_beta_infection_nob4 <- rmse_beta}
    ## plot residuals    
    beta_m <- betareg(prop~batch+age+gender, data=mdata_nob4, link='logit') %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(celltype, condition, residuals)
    ggplot(beta_m, aes(x=condition, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('NI', conditions[i]))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_infection_nob4.png', height=3, width=4)
    
    # fit beta model in income, with batch 4 
    beta_m <- betareg(prop~batch+age+gender+condition*income, data=mdata, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='condition'%&%conditions[i]%&%':incomeHigh')
    if (exists('compiled_beta_interaction_income')){
      compiled_beta_interaction_income <- rbind(compiled_beta_interaction_income, results_beta_m)
    } else {compiled_beta_interaction_income <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+condition*income, data=mdata, link='logit'), type='response')
    rmse_beta <- rmse(mdata$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='condition'%&%conditions[i]%&%':incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_interaction_income')){
      rmse_compiled_beta_interaction_income <- rbind(rmse_compiled_beta_interaction_income, rmse_beta)
    } else {rmse_compiled_beta_interaction_income <- rmse_beta}
    ## plot residuals   
    beta_m <- betareg(prop~batch+age+gender, data=mdata, link='logit') %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(IDs, condition, celltype, income, residuals)%>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, income, residuals)
    ggplot(beta_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_infection.income.png', height=3, width=4)
    
    # fit beta model in income, without batch 4 
    beta_m <- betareg(prop~batch+age+gender+condition*income, data=mdata_nob4, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='condition'%&%conditions[i]%&%':incomeHigh')
    if (exists('compiled_beta_interaction_income_nob4')){
      compiled_beta_interaction_income_nob4 <- rbind(compiled_beta_interaction_income_nob4, results_beta_m)
    } else {compiled_beta_interaction_income_nob4 <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+condition*income, data=mdata_nob4, link='logit'), type='response')
    rmse_beta <- rmse(mdata_nob4$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='condition'%&%conditions[i]%&%':incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_interaction_income_nob4')){
      rmse_compiled_beta_interaction_income_nob4 <- rbind(rmse_compiled_beta_interaction_income_nob4, rmse_beta)
    } else {rmse_compiled_beta_interaction_income_nob4 <- rmse_beta}
    ## plot residuals   
    beta_m <- betareg(prop~batch+age+gender, data=mdata_nob4, link='logit') %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(IDs, condition, celltype, income, residuals)%>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, income, residuals)
    ggplot(beta_m, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_infection.income_nob4.png', height=3, width=4)
    
    # fit beta model in asthma, with batch 4 
    beta_m <- betareg(prop~batch+age+gender+condition*asthma, data=mdata, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes') %>%
      filter(terms=='condition'%&%conditions[i]%&%':asthmaYes')
    if (exists('compiled_beta_interaction_asthma')){
      compiled_beta_interaction_asthma <- rbind(compiled_beta_interaction_asthma, results_beta_m)
    } else {compiled_beta_interaction_asthma <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+condition*asthma, data=mdata, link='logit'), type='response')
    rmse_beta <- rmse(mdata$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='yes', terms='condition'%&%conditions[i]%&%':asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_interaction_asthma')){
      rmse_compiled_beta_interaction_asthma <- rbind(rmse_compiled_beta_interaction_asthma, rmse_beta)
    } else {rmse_compiled_beta_interaction_asthma <- rmse_beta}
    ## plot residuals   
    beta_m <- betareg(prop~batch+age+gender, data=mdata, link='logit') %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata, by=c('orig.ident')) %>% select(IDs, condition, celltype, asthma, residuals)%>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, asthma, residuals)
    ggplot(beta_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_infection.asthma.png', height=3, width=4)
    
    # fit beta model in asthma, without batch 4 
    beta_m <- betareg(prop~batch+age+gender+condition*asthma, data=mdata_nob4, link='logit') %>% summary()
    betas <- beta_m$coefficients$mean[,1]
    pvals <- beta_m$coefficients$mean[,4]
    terms <- rownames(beta_m$coefficients$mean)
    results_beta_m <- data.frame(terms, betas, pvals) %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no') %>%
      filter(terms=='condition'%&%conditions[i]%&%':asthmaYes')
    if (exists('compiled_beta_interaction_asthma_nob4')){
      compiled_beta_interaction_asthma_nob4 <- rbind(compiled_beta_interaction_asthma_nob4, results_beta_m)
    } else {compiled_beta_interaction_asthma_nob4 <- results_beta_m}
    # compute RMSE
    beta_m <- predict(betareg(prop~batch+age+gender+condition*asthma, data=mdata_nob4, link='logit'), type='response')
    rmse_beta <- rmse(mdata_nob4$prop, beta_m) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=ctype, batch4='no', terms='condition'%&%conditions[i]%&%':asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_beta_interaction_asthma_nob4')){
      rmse_compiled_beta_interaction_asthma_nob4 <- rbind(rmse_compiled_beta_interaction_asthma_nob4, rmse_beta)
    } else {rmse_compiled_beta_interaction_asthma_nob4 <- rmse_beta}
    ## plot residuals   
    beta_m <- betareg(prop~batch+age+gender, data=mdata_nob4, link='logit') %>% residuals(type='pearson') %>% 
      as.data.frame() %>% rename(residuals='.') %>% rownames_to_column('orig.ident') %>% 
      left_join(mdata_nob4, by=c('orig.ident')) %>% select(IDs, condition, celltype, asthma, residuals) %>%
      pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
      select(celltype, asthma, residuals)
    ggplot(beta_m, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y='Delta '%&%ctype%&%' proportion (residualsINF - residualsNI)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_beta_infection.asthma_nob4.png', height=3, width=4)
  }
}

###########################################
### MULTIPLE DIRICHLET REGRESSION MODEL ###
###########################################
# first, test asthma and income per infection
conditions <- c('NI', 'RV', 'IVA')
eps <- 1e-6
for (i in 1:length(conditions)){

    # extract and format metadata
    mdata <- obj@meta.data %>% filter((condition==conditions[i])) %>% select(-c(orig.ident, n, avg_mt))
    mdata <- mdata %>% pivot_wider(names_from='celltype', values_from='prop', values_fill=0)
    # add small pseudocount to 0s and renormalize
    mdata[celltypes] <- mdata[celltypes] + eps
    mdata[celltypes] <- mdata[celltypes] / rowSums(mdata[celltypes])
    mdata$Y <- DR_data(mdata[,10:ncol(mdata)])
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$albuterol <- na_if(mdata$albuterol, '')
    mdata$albuterol <- factor(mdata$albuterol, levels=c('No', 'Yes'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Low', 'High')
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    mdata_nob4 <- mdata %>% filter(batch!='B4')

    # fit Dirichlet model in income, with batch 4 
    dr_m <- DirichReg(Y~batch+age+gender+income, data=mdata) %>% summary()
    dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
      filter(str_detect(rowname, 'incomeHigh')) %>% 
      mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='yes',
             terms='incomeHigh') %>% select(terms, Estimate, `Pr(>|z|)`,
                                            condition, celltype, batch4) %>%
      rename(betas=Estimate, pvals=`Pr(>|z|)`)
    if (exists('compiled_dr_income')){
      compiled_dr_income <- rbind(compiled_dr_income, dr_m_coeffs)
    } else {compiled_dr_income <- dr_m_coeffs}
    # compute RMSE
    dr_m <- predict(DirichReg(Y~batch+age+gender+income, data=mdata)) %>% as.data.frame()
    for (j in 1:length(celltypes)){
      rmse_dr <- rmse(mdata[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
        mutate(condition=conditions[i], celltype=celltypes[j], batch4='yes', terms='incomeHigh') %>%
        rename(rmse='.')
      if (exists('rmse_compiled_dr_income')){
        rmse_compiled_dr_income <- rbind(rmse_compiled_dr_income, rmse_dr)
      } else {rmse_compiled_dr_income <- rmse_dr}
    }
    ## plot residuals   
    dr_res <- DirichReg(Y~batch+age+gender, data=mdata) %>% residuals(type='standardized') %>% 
      unclass() %>% as.data.frame() %>% mutate(IDs=mdata$IDs) %>%
      pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata, by=c('IDs')) %>% 
      select(condition, celltype, income, residuals) 
    for (ctype in celltypes){
      dr_tmp <- dr_res %>% filter(celltype==ctype)
      ggplot(dr_tmp, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
        stat_compare_means(comparisons=list(c('Low', 'High'))) + 
        scale_y_continuous(expand=expansion(mult=.3)) + 
        labs(y=ctype%&%' proportion (residuals)')
      ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_income.png', height=3, width=4)
    }

    # fit Dirichlet model in income, without batch 4 
    dr_m <- DirichReg(Y~batch+age+gender+income, data=mdata_nob4) %>% summary()
    dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
      filter(str_detect(rowname, 'incomeHigh')) %>% 
      mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='no',
             terms='incomeHigh') %>% select(terms, Estimate, `Pr(>|z|)`,
                                            condition, celltype, batch4) %>%
      rename(betas=Estimate, pvals=`Pr(>|z|)`)
    if (exists('compiled_dr_income_nob4')){
      compiled_dr_income_nob4 <- rbind(compiled_dr_income_nob4, dr_m_coeffs)
    } else {compiled_dr_income_nob4 <- dr_m_coeffs}
    # compute RMSE
    dr_m <- predict(DirichReg(Y~batch+age+gender+income, data=mdata_nob4)) %>% as.data.frame()
    for (j in 1:length(celltypes)){
      rmse_dr <- rmse(mdata_nob4[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
        mutate(condition=conditions[i], celltype=celltypes[j], batch4='no', terms='incomeHigh') %>%
        rename(rmse='.')
      if (exists('rmse_compiled_dr_income_nob4')){
        rmse_compiled_dr_income_nob4 <- rbind(rmse_compiled_dr_income_nob4, rmse_dr)
      } else {rmse_compiled_dr_income_nob4 <- rmse_dr}
    }
    ## plot residuals 
    dr_res <- DirichReg(Y~batch+age+gender, data=mdata_nob4) %>% residuals(type='standardized') %>% 
      unclass() %>% as.data.frame() %>% mutate(IDs=mdata_nob4$IDs) %>%
      pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata_nob4, by=c('IDs')) %>% 
      select(condition, celltype, income, residuals) 
    for (ctype in celltypes){
      dr_tmp <- dr_res %>% filter(celltype==ctype)
      ggplot(dr_tmp, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
        stat_compare_means(comparisons=list(c('Low', 'High'))) + 
        scale_y_continuous(expand=expansion(mult=.3)) + 
        labs(y=ctype%&%' proportion (residuals)')
      ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_income_nob4.png', height=3, width=4)
    }

    # fit Dirichlet model in asthma, with batch 4 
    dr_m <- DirichReg(Y~batch+age+gender+asthma, data=mdata) %>% summary()
    dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
      filter(str_detect(rowname, 'asthmaYes')) %>% 
      mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='yes',
             terms='asthmaYes') %>% select(terms, Estimate, `Pr(>|z|)`,
                                            condition, celltype, batch4) %>%
      rename(betas=Estimate, pvals=`Pr(>|z|)`)
    if (exists('compiled_dr_asthma')){
      compiled_dr_asthma <- rbind(compiled_dr_asthma, dr_m_coeffs)
    } else {compiled_dr_asthma <- dr_m_coeffs}
    # compute RMSE
    dr_m <- predict(DirichReg(Y~batch+age+gender+asthma, data=mdata)) %>% as.data.frame()
    for (j in 1:length(celltypes)){
      rmse_dr <- rmse(mdata[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
        mutate(condition=conditions[i], celltype=celltypes[j], batch4='yes', terms='asthmaYes') %>%
        rename(rmse='.')
      if (exists('rmse_compiled_dr_asthma')){
        rmse_compiled_dr_asthma <- rbind(rmse_compiled_dr_asthma, rmse_dr)
      } else {rmse_compiled_dr_asthma <- rmse_dr}
    }
    ## plot residuals     
    dr_res <- DirichReg(Y~batch+age+gender, data=mdata) %>% residuals(type='standardized') %>% 
      unclass() %>% as.data.frame() %>% mutate(IDs=mdata$IDs) %>%
      pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata, by=c('IDs')) %>% 
      select(condition, celltype, asthma, residuals) 
    for (ctype in celltypes){
      dr_tmp <- dr_res %>% filter(celltype==ctype)
      ggplot(dr_tmp, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
        stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
        scale_y_continuous(expand=expansion(mult=.3)) + 
        labs(y=ctype%&%' proportion (residuals)')
      ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_asthma.png', height=3, width=4)
    }
    
    # fit Dirichlet model in asthma, without batch 4 
    dr_m <- DirichReg(Y~batch+age+gender+asthma, data=mdata_nob4) %>% summary()
    dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
      filter(str_detect(rowname, 'asthmaYes')) %>% 
      mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='no',
             terms='asthmaYes') %>% select(terms, Estimate, `Pr(>|z|)`,
                                            condition, celltype, batch4) %>%
      rename(betas=Estimate, pvals=`Pr(>|z|)`)
    if (exists('compiled_dr_asthma_nob4')){
      compiled_dr_asthma_nob4 <- rbind(compiled_dr_asthma_nob4, dr_m_coeffs)
    } else {compiled_dr_asthma_nob4 <- dr_m_coeffs}
    # compute RMSE
    dr_m <- predict(DirichReg(Y~batch+age+gender+asthma, data=mdata_nob4)) %>% as.data.frame()
    for (j in 1:length(celltypes)){
      rmse_dr <- rmse(mdata_nob4[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
        mutate(condition=conditions[i], celltype=celltypes[j], batch4='no', terms='asthmaYes') %>%
        rename(rmse='.')
      if (exists('rmse_compiled_dr_asthma_nob4')){
        rmse_compiled_dr_asthma_nob4 <- rbind(rmse_compiled_dr_asthma_nob4, rmse_dr)
      } else {rmse_compiled_dr_asthma_nob4 <- rmse_dr}
    }
    ## plot residuals    
    dr_res <- DirichReg(Y~batch+age+gender, data=mdata_nob4) %>% residuals(type='standardized') %>% 
      unclass() %>% as.data.frame() %>% mutate(IDs=mdata_nob4$IDs) %>%
      pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata_nob4, by=c('IDs')) %>% 
      select(condition, celltype, asthma, residuals) 
    for (ctype in celltypes){
      dr_tmp <- dr_res %>% filter(celltype==ctype)
      ggplot(dr_tmp, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
        stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
        scale_y_continuous(expand=expansion(mult=.3)) + 
        labs(y=ctype%&%' proportion (residuals)')
      ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_asthma_nob4.png', height=3, width=4)
    }
}

###########################################
### MULTIPLE DIRICHLET REGRESSION MODEL ###
###########################################
# now, test asthma and income interacting with infection
conditions <- c('RV', 'IVA')
eps <- 1e-6
for (i in 1:length(conditions)){
  
  # extract and format metadata
  mdata <- obj@meta.data %>% filter((condition==conditions[i]|condition=='NI')) %>% select(-c(orig.ident, n, avg_mt))
  mdata <- mdata %>% pivot_wider(names_from='celltype', values_from='prop', values_fill=0)
  # add small pseudocount to 0s and renormalize
  mdata[celltypes] <- mdata[celltypes] + eps
  mdata[celltypes] <- mdata[celltypes] / rowSums(mdata[celltypes])
  mdata$Y <- DR_data(mdata[,10:ncol(mdata)])
  mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
  mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
  mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
  mdata$income <- na_if(mdata$income, '')
  mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                         'Low', 'High')
  mdata$income <- factor(mdata$income, levels=c('Low','High'))
  mdata <- mdata %>% filter(IDs %in% IDs[duplicated(IDs)])    # remove IDs that are not paired
  mdata_nob4 <- mdata %>% filter(batch!='B4')
  
  # fit Dirichlet model in infection, with batch 4 
  dr_m <- DirichReg(Y~batch+age+gender+condition, data=mdata) %>% summary()
  dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
    filter(str_detect(rowname, 'condition'%&%conditions[i])) %>% 
    mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='yes',
           terms='infection') %>% select(terms, Estimate, `Pr(>|z|)`,
                                          condition, celltype, batch4) %>%
    rename(betas=Estimate, pvals=`Pr(>|z|)`)
  if (exists('compiled_dr_infection')){
    compiled_dr_infection <- rbind(compiled_dr_infection, dr_m_coeffs)
  } else {compiled_dr_infection <- dr_m_coeffs}
  # compute RMSE
  dr_m <- predict(DirichReg(Y~batch+age+gender+condition, data=mdata)) %>% as.data.frame()
  for (j in 1:length(celltypes)){
    rmse_dr <- rmse(mdata[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=celltypes[j], batch4='yes', terms='infection') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_dr_infection')){
      rmse_compiled_dr_infection <- rbind(rmse_compiled_dr_infection, rmse_dr)
    } else {rmse_compiled_dr_infection <- rmse_dr}
  }
  ## plot residuals    
  dr_res <- DirichReg(Y~batch+age+gender, data=mdata) %>% residuals(type='standardized') %>% 
    unclass() %>% as.data.frame() %>% mutate(IDs=mdata$IDs, condition=mdata$condition) %>%
    pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata, by=c('IDs', 'condition')) %>% 
    select(condition, celltype, residuals) 
  for (ctype in celltypes){
    dr_tmp <- dr_res %>% filter(celltype==ctype)
    ggplot(dr_tmp, aes(x=condition, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('NI', conditions[i]))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_infection.png', height=3, width=4)
  }
  
  # fit Dirichlet model in infection, without batch 4 
  dr_m <- DirichReg(Y~batch+age+gender+condition, data=mdata_nob4) %>% summary()
  dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
    filter(str_detect(rowname, 'condition'%&%conditions[i])) %>% 
    mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='no',
           terms='infection') %>% select(terms, Estimate, `Pr(>|z|)`,
                                         condition, celltype, batch4) %>%
    rename(betas=Estimate, pvals=`Pr(>|z|)`)
  if (exists('compiled_dr_infection_nob4')){
    compiled_dr_infection_nob4 <- rbind(compiled_dr_infection_nob4, dr_m_coeffs)
  } else {compiled_dr_infection_nob4 <- dr_m_coeffs}
  # compute RMSE
  dr_m <- predict(DirichReg(Y~batch+age+gender+condition, data=mdata_nob4)) %>% as.data.frame()
  for (j in 1:length(celltypes)){
    rmse_dr <- rmse(mdata_nob4[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=celltypes[j], batch4='no', terms='infection') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_dr_infection_nob4')){
      rmse_compiled_dr_infection_nob4 <- rbind(rmse_compiled_dr_infection_nob4, rmse_dr)
    } else {rmse_compiled_dr_infection_nob4 <- rmse_dr}
  }
  ## plot residuals   
  dr_res <- DirichReg(Y~batch+age+gender, data=mdata_nob4) %>% residuals(type='standardized') %>% 
    unclass() %>% as.data.frame() %>% mutate(IDs=mdata_nob4$IDs, condition=mdata_nob4$condition) %>%
    pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata_nob4, by=c('IDs', 'condition')) %>% 
    select(condition, celltype, residuals) 
  for (ctype in celltypes){
    dr_tmp <- dr_res %>% filter(celltype==ctype)
    ggplot(dr_tmp, aes(x=condition, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('NI', conditions[i]))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_infection_nob4.png', height=3, width=4)
  }

  # fit Dirichlet model in income, with batch 4 
  dr_m <- DirichReg(Y~batch+age+gender+condition*income, data=mdata) %>% summary()
  dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
    filter(str_detect(rowname, 'condition'%&%conditions[i]%&%'.incomeHigh')) %>% 
    mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='yes',
           terms='condition'%&%conditions[i]%&%':incomeHigh') %>% select(terms, Estimate, `Pr(>|z|)`,
                                          condition, celltype, batch4) %>%
    rename(betas=Estimate, pvals=`Pr(>|z|)`)
  if (exists('compiled_dr_interaction_income')){
    compiled_dr_interaction_income <- rbind(compiled_dr_interaction_income, dr_m_coeffs)
  } else {compiled_dr_interaction_income <- dr_m_coeffs}
  # compute RMSE
  dr_m <- predict(DirichReg(Y~batch+age+gender+condition*income, data=mdata)) %>% as.data.frame()
  for (j in 1:length(celltypes)){
    rmse_dr <- rmse(mdata[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=celltypes[j], batch4='yes', terms='condition'%&%conditions[i]%&%':incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_dr_interaction_income')){
      rmse_compiled_dr_interaction_income <- rbind(rmse_compiled_dr_interaction_income, rmse_dr)
    } else {rmse_compiled_dr_interaction_income <- rmse_dr}
  }
  ## plot residuals  
  dr_res <- DirichReg(Y~batch+age+gender, data=mdata) %>% residuals(type='standardized') %>% 
    unclass() %>% as.data.frame() %>% mutate(IDs=mdata$IDs, condition=mdata$condition) %>%
    pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata, by=c('IDs', 'condition')) %>% 
    select(IDs, condition, celltype, income, residuals)%>%
    pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
    select(celltype, income, residuals)
  for (ctype in celltypes){
    dr_tmp <- dr_res %>% filter(celltype==ctype)
    ggplot(dr_tmp, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_infection.income.png', height=3, width=4)
  }
  
  # fit Dirichlet model in income, without batch 4 
  dr_m <- DirichReg(Y~batch+age+gender+condition*income, data=mdata_nob4) %>% summary()
  dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
    filter(str_detect(rowname, 'condition'%&%conditions[i]%&%'.incomeHigh')) %>% 
    mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='no',
           terms='condition'%&%conditions[i]%&%':incomeHigh') %>% select(terms, Estimate, `Pr(>|z|)`,
                                          condition, celltype, batch4) %>%
    rename(betas=Estimate, pvals=`Pr(>|z|)`)
  if (exists('compiled_dr_interaction_income_nob4')){
    compiled_dr_interaction_income_nob4 <- rbind(compiled_dr_interaction_income_nob4, dr_m_coeffs)
  } else {compiled_dr_interaction_income_nob4 <- dr_m_coeffs}
  # compute RMSE
  dr_m <- predict(DirichReg(Y~batch+age+gender+condition*income, data=mdata_nob4)) %>% as.data.frame()
  for (j in 1:length(celltypes)){
    rmse_dr <- rmse(mdata_nob4[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=celltypes[j], batch4='no', terms='condition'%&%conditions[i]%&%':incomeHigh') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_dr_interaction_income_nob4')){
      rmse_compiled_dr_interaction_income_nob4 <- rbind(rmse_compiled_dr_interaction_income_nob4, rmse_dr)
    } else {rmse_compiled_dr_interaction_income_nob4 <- rmse_dr}
  }
  ## plot residuals  
  dr_res <- DirichReg(Y~batch+age+gender, data=mdata_nob4) %>% residuals(type='standardized') %>% 
    unclass() %>% as.data.frame() %>% mutate(IDs=mdata_nob4$IDs, condition=mdata_nob4$condition) %>%
    pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata_nob4, by=c('IDs', 'condition')) %>% 
    select(IDs, condition, celltype, income, residuals)%>%
    pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
    select(celltype, income, residuals)
  for (ctype in celltypes){
    dr_tmp <- dr_res %>% filter(celltype==ctype)
    ggplot(dr_tmp, aes(x=income, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('Low', 'High'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_infection.income_nob4.png', height=3, width=4)
  }
  
  # fit Dirichlet model in asthma, with batch 4 
  dr_m <- DirichReg(Y~batch+age+gender+condition*asthma, data=mdata) %>% summary()
  dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
    filter(str_detect(rowname, 'condition'%&%conditions[i]%&%'.asthmaYes')) %>% 
    mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='yes',
           terms='condition'%&%conditions[i]%&%':asthmaYes') %>% select(terms, Estimate, `Pr(>|z|)`,
                                                                         condition, celltype, batch4) %>%
    rename(betas=Estimate, pvals=`Pr(>|z|)`)
  if (exists('compiled_dr_interaction_asthma')){
    compiled_dr_interaction_asthma <- rbind(compiled_dr_interaction_asthma, dr_m_coeffs)
  } else {compiled_dr_interaction_asthma <- dr_m_coeffs}
  # compute RMSE
  dr_m <- predict(DirichReg(Y~batch+age+gender+condition*asthma, data=mdata)) %>% as.data.frame()
  for (j in 1:length(celltypes)){
    rmse_dr <- rmse(mdata[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=celltypes[j], batch4='yes', terms='condition'%&%conditions[i]%&%':asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_dr_interaction_asthma')){
      rmse_compiled_dr_interaction_asthma <- rbind(rmse_compiled_dr_interaction_asthma, rmse_dr)
    } else {rmse_compiled_dr_interaction_asthma <- rmse_dr}
  }
  ## plot residuals  
  dr_res <- DirichReg(Y~batch+age+gender, data=mdata) %>% residuals(type='standardized') %>% 
    unclass() %>% as.data.frame() %>% mutate(IDs=mdata$IDs, condition=mdata$condition) %>%
    pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata, by=c('IDs', 'condition')) %>% 
    select(IDs, condition, celltype, asthma, residuals)%>%
    pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
    select(celltype, asthma, residuals)
  for (ctype in celltypes){
    dr_tmp <- dr_res %>% filter(celltype==ctype)
    ggplot(dr_tmp, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_infection.asthma.png', height=3, width=4)
  }
  
  # fit Dirichlet model in asthma, without batch 4 
  dr_m <- DirichReg(Y~batch+age+gender+condition*asthma, data=mdata_nob4) %>% summary()
  dr_m_coeffs <- dr_m$coef.mat %>% as.data.frame() %>% rownames_to_column() %>%
    filter(str_detect(rowname, 'condition'%&%conditions[i]%&%'.asthmaYes')) %>% 
    mutate(condition=conditions[i], celltype=dr_m$varnames, batch4='no',
           terms='condition'%&%conditions[i]%&%':asthmaYes') %>% select(terms, Estimate, `Pr(>|z|)`,
                                                                         condition, celltype, batch4) %>%
    rename(betas=Estimate, pvals=`Pr(>|z|)`)
  if (exists('compiled_dr_interaction_asthma_nob4')){
    compiled_dr_interaction_asthma_nob4 <- rbind(compiled_dr_interaction_asthma_nob4, dr_m_coeffs)
  } else {compiled_dr_interaction_asthma_nob4 <- dr_m_coeffs}
  # compute RMSE
  dr_m <- predict(DirichReg(Y~batch+age+gender+condition*asthma, data=mdata_nob4)) %>% as.data.frame()
  for (j in 1:length(celltypes)){
    rmse_dr <- rmse(mdata_nob4[[15]][,celltypes[j]], dr_m[,j]) %>% as.data.frame() %>% 
      mutate(condition=conditions[i], celltype=celltypes[j], batch4='no', terms='condition'%&%conditions[i]%&%':asthmaYes') %>%
      rename(rmse='.')
    if (exists('rmse_compiled_dr_interaction_asthma_nob4')){
      rmse_compiled_dr_interaction_asthma_nob4 <- rbind(rmse_compiled_dr_interaction_asthma_nob4, rmse_dr)
    } else {rmse_compiled_dr_interaction_asthma_nob4 <- rmse_dr}
  }
  ## plot residuals  
  dr_res <- DirichReg(Y~batch+age+gender, data=mdata_nob4) %>% residuals(type='standardized') %>% 
    unclass() %>% as.data.frame() %>% mutate(IDs=mdata_nob4$IDs, condition=mdata_nob4$condition) %>%
    pivot_longer(cols=celltypes, names_to='celltype', values_to='residuals') %>% left_join(mdata_nob4, by=c('IDs', 'condition')) %>% 
    select(IDs, condition, celltype, asthma, residuals)%>%
    pivot_wider(names_from=condition, values_from=residuals) %>% mutate(residuals=.data[[conditions[i]]]-NI) %>%
    select(celltype, asthma, residuals)
  for (ctype in celltypes){
    dr_tmp <- dr_res %>% filter(celltype==ctype)
    ggplot(dr_tmp, aes(x=asthma, y=residuals)) + geom_boxplot() + theme_bw() +
      stat_compare_means(comparisons=list(c('No', 'Yes'))) + 
      scale_y_continuous(expand=expansion(mult=.3)) + 
      labs(y=ctype%&%' proportion (residuals)')
    ggsave('proportion_plots/'%&%conditions[i]%&%'_'%&%ctype%&%'_dr_infection.asthma_nob4.png', height=3, width=4)
  }
}

# perform multiple testing correction on each model separately 
list_of_results <- mget(ls(pattern='^compiled'))
for (j in seq(1:length(list_of_results))){
  list_of_results[[j]] <- list_of_results[[j]] %>% mutate(p_adj=p.adjust(pvals, method='fdr'))
}
single_results <- do.call(rbind, list_of_results) %>% rownames_to_column('regression') %>% 
  mutate(regression = case_when(
    grepl('_beta_', regression)   ~ 'beta',
    grepl('_dr_', regression)     ~ 'DR',
    grepl('_linear_', regression) ~ 'linear',
    TRUE                          ~ NA_character_
  )
) 
single_results <- single_results %>% mutate(ultra_p_adj=p.adjust(pvals, method='fdr')) %>%
  arrange(terms, condition, celltype, batch4) 
sig_single_resutls <- single_results %>% filter(p_adj<0.05)
ultra_sig_single_resutls <- single_results %>% filter(ultra_p_adj<0.05)


# sig. in DR but not in beta or linear:
# asthma or income changes the relative composition of cell types (compositional dependencies captured by DR), 
# without producing a strong marginal change in any single celltype’s absolute proportion
# sig. in beta or linear but not in DR:
# ? 


sig_single_resutls %>% filter(batch4=='yes') %>% ggplot(., aes(x=celltype, y=betas, color=regression)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(terms), rows=vars(condition)) + ggtitle('batch4')
sig_single_resutls %>% filter(batch4=='no') %>% ggplot(., aes(x=celltype, y=betas, color=regression)) + geom_point() +
  theme_bw() + facet_grid(cols=vars(terms), rows=vars(condition)) + ggtitle('no_batch4')

# compare RMSE between models
list_of_rmse <- mget(ls(pattern='^rmse_compiled'))
all_RMSEs <- do.call(rbind, list_of_rmse) %>% rownames_to_column('regression') %>% 
  mutate(regression = case_when(
    grepl('_beta_', regression)   ~ 'beta',
    grepl('_dr_', regression)     ~ 'DR',
    grepl('_linear_', regression) ~ 'linear',
    TRUE                          ~ NA_character_
  ), terms = case_when(
    grepl(':asthmaYes', terms)   ~ 'infection:asthmaYes',
    grepl(':incomeHigh', terms)   ~ 'infection:incomeHigh',
    grepl('asthmaYes', terms) ~ 'asthmaYes',
    grepl('incomeHigh', terms) ~ 'incomeHigh',
    grepl('infection', terms) ~ 'infection',
    TRUE                          ~ NA_character_
  )
)
all_RMSEs$regression <- factor(all_RMSEs$regression, levels=c('linear', 'beta', 'DR'))
all_RMSEs$condition <- factor(all_RMSEs$condition, levels=c('NI', 'IVA', 'RV'))

all_RMSEs %>% filter(batch4=='yes') %>% ggplot(., aes(x=celltype, y=rmse, fill=regression)) + geom_col(position='dodge') + theme_bw() +
  facet_grid(cols=vars(terms), rows=vars(condition))
ggsave('proportion_plots/RMSE_yesb4.png', height=5, width=12)
all_RMSEs %>% filter(batch4=='no') %>% ggplot(., aes(x=celltype, y=rmse, fill=regression)) + geom_col(position='dodge') + theme_bw() +
  facet_grid(cols=vars(terms), rows=vars(condition))
ggsave('proportion_plots/RMSE_nob4.png', height=5, width=12)

all_RMSEs %>% filter(batch4=='yes', regression!='DR') %>% 
  ggplot(., aes(x=celltype, y=rmse, fill=regression)) + geom_col(position='dodge') + theme_bw() +
  facet_grid(cols=vars(terms), rows=vars(condition))
ggsave('proportion_plots/RMSE_noDR_yesb4.png', height=5, width=12)
all_RMSEs %>% filter(batch4=='no', regression!='DR') %>% 
  ggplot(., aes(x=celltype, y=rmse, fill=regression)) + geom_col(position='dodge') + theme_bw() +
  facet_grid(cols=vars(terms), rows=vars(condition))
ggsave('proportion_plots/RMSE_noDR_nob4.png', height=5, width=12)
