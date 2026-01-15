library(Seurat)
library(tidyverse)
library(edgeR)
library(limma)
library(broom)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('RV', 'IVA')
celltypes <- c('B','CD4-T','CD8-T','Mono','NK')

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
ggsave('celltype_prop_simple.png', height=5, width=5)
ggsave('celltype_prop_simple.pdf', height=5, width=5)

## same thing, but without batch 4
avg_df <- joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection) %>%
  summarise(prop = mean(prop), .groups = 'drop')
joint_df %>% filter(batch!='B4') %>% ggplot(., aes(x=condition, y=prop, group=IDs)) + 
  geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=1, color='red'), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop), color='red', size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() + theme(legend.position='none')
ggsave('celltype_prop_simple_noB4.png', height=5, width=5)
ggsave('celltype_prop_simple_noB4.pdf', height=5, width=5)

# proportion of each cell type per NI or Inf group per asthma status
avg_df <- joint_df %>% group_by(condition, celltype, infection, asthma) %>%
  summarise(prop = mean(prop), .groups = 'drop')
ggplot(joint_df, aes(x=condition, y=prop, group=IDs, color=asthma)) + geom_point(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=avg_df, aes(group=asthma, color=asthma), linewidth=2) +
  geom_point(data=avg_df, aes(x=condition, y=prop, color=asthma), size=3, inherit.aes=FALSE) +
  facet_grid(rows=vars(celltype), cols=vars(infection), scales='free') + theme_bw() 
ggsave('celltype_prop_simple_asthma.png', height=5, width=5)
ggsave('celltype_prop_simple_asthma.pdf', height=5, width=5)

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
ggsave('celltype_prop_simple_asthma_noB4.png', height=5, width=5)
ggsave('celltype_prop_simple_asthma_noB4.pdf', height=5, width=5)

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
ggsave('celltype_prop_simple_income.png', height=5, width=5)
ggsave('celltype_prop_simple_income.pdf', height=5, width=5)

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
ggsave('celltype_prop_simple_income_noB4.png', height=5, width=5)
ggsave('celltype_prop_simple_income_noB4.pdf', height=5, width=5)

# test if there are differences in proportion
joint_df %>% filter(batch!='B4') %>% group_by(condition, celltype, infection) %>% filter(n_distinct(income)==2) %>% 
  do(tidy(t.test(prop ~ income, data = .))) %>% ungroup() %>% mutate(p_adj=p.adjust(p.value, method='fdr')) %>%
  filter(p.value<0.05)

for (i in 1:length(conditions)){
  for (ctype in celltypes){

    # extract and format metadata
    mdata <- obj@meta.data %>% filter(celltype==ctype, (condition==conditions[i]|condition=='NI'))
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$albuterol <- na_if(mdata$albuterol, '')
    mdata$albuterol <- factor(mdata$albuterol, levels=c('No', 'Yes'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- ifelse(mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                       'Low', 'High')
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    
    # remove IDs that are not paired
    mdata <- mdata %>% filter(IDs %in% IDs[duplicated(IDs)])
    
    # define design matrix for infection
    design <- model.matrix(~batch+age+gender+avg_mt+condition, data=mdata)
    
    # fit linear model
    a <- lm(mdata$prop ~ design) %>% summary()
  }
}












