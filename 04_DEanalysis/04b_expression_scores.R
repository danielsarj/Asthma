library(tidyverse)
library(data.table)
library(Seurat)
library(limma)
library(msigdbr)
'%&%' <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('IVA', 'RV')
celltypes <- c('B','CD4-T','CD8-T','Mono','NK')
interactions <- c('none', 'asthma_alb', 'income', 'ACT', 'ACE', 'resilience', 'social_support', 'total_racism',
                  'year_racism', 'life_racism', 'stress_racism', 'kid_24h_racism', 'kid_discrimination', 'collection_infection')

# get human hallmark gene sets
genes_hallmark <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)
names(genes_hallmark) <- gsub('HALLMARK_', '', names(genes_hallmark))

# load pseudobulk seurat object
bulk_obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks_new.rds')

# remove batch 4
bulk_obj <- subset(bulk_obj, subset= batch!='B4')
bulk_obj@meta.data$condition <- factor(bulk_obj@meta.data$condition, levels=c('NI','IVA','RV'))
bulk_obj@meta.data$gender <- factor(bulk_obj@meta.data$gender, levels=c('Male','Female'))
bulk_obj@meta.data$batch <- factor(bulk_obj@meta.data$batch, levels=c('B1','B2','B3'))
bulk_obj@meta.data$income <- factor(bulk_obj@meta.data$income, levels=c('Low','High'))
bulk_obj@meta.data$albuterol <- factor(bulk_obj@meta.data$albuterol, levels=c('No','Yes'))
bulk_obj@meta.data$Recorded_Diagnosis <- factor(bulk_obj@meta.data$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
bulk_obj@meta.data$infection_status <- factor(bulk_obj@meta.data$infection_status, levels=c('Negative', 'Positive'))


compute_pathway_expression <- function(seurat_obj, cond, ctype, int) {
  # subset seurat object 
  tmp <- subset(seurat_obj, celltype==ctype & (condition=='NI' | condition==cond))
  
  # get samples/genes from voom expression dataframe (easier to filter stuff)
  if (int=='none'){
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_voom_expression_new.txt.gz') %>% column_to_rownames('Gene')
  } else {
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_'%&%int%&%'_voom_expression_new.txt.gz') %>% column_to_rownames('Gene')
  }
  v_samples <- colnames(v)
  v_genes <- rownames(v)
  rm(v)
  
  # subset metadata
  sub_mdata <- tmp@meta.data %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
  
  # extract and subset count matrix
  count <- tmp@assays$RNA$counts
  count <- count[rownames(count) %in% v_genes,]
  count <- count[,colnames(count) %in% v_samples]
  count <- DGEList(counts=count) %>% calcNormFactors()
  
  # create design matrix
  if (int=='asthma_alb'){
    design <- model.matrix(~batch+age+gender+n+avg_mt+prop+albuterol, data=sub_mdata)
  } else {
    design <- model.matrix(~batch+age+gender+n+avg_mt+prop, data=sub_mdata)
  }
  
  # fit linear model
  voom_obj <- voom(count, design, plot=F)
  l_model <- lmFit(voom_obj, design=design) %>% eBayes()
  
  # extract residuals and add the intercept
  residuals <- residuals.MArrayLM(l_model, voom_obj)
  intercept <- l_model$coefficients[,'(Intercept)']
  corrected_expression <- residuals + intercept
  
  # average expression per interaction 
  if (int=='none'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      select(Gene, IDs, condition, celltype, value) %>% group_by(Gene, condition, celltype) %>% 
      summarise(mean_value=mean(value)) %>% ungroup() 
    exp_long$condition <- gsub('NI', 'NI_'%&%tolower(cond), exp_long$condition)
  } else if (int=='asthma_alb'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Recorded_Diagnosis) %>% group_by(Gene, condition, celltype, Recorded_Diagnosis) %>% 
      summarise(mean_value=mean(value)) %>%
      ungroup() %>% group_by(Gene, celltype, Recorded_Diagnosis) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='income'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, income) %>% group_by(Gene, condition, celltype, income) %>% 
      summarise(mean_value=mean(value)) %>%
      ungroup() %>% group_by(Gene, celltype, income) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='ACT'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, ACT_score) %>% mutate(quartile=factor(ntile(ACT_score, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='ACE'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, ACE_result) %>% mutate(quartile=factor(ntile(ACE_result, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='resilience'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Parent_Resilience_Score) %>% mutate(quartile=factor(ntile(Parent_Resilience_Score, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='social_support'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Parents_Score_Avg) %>% mutate(quartile=factor(ntile(Parents_Score_Avg, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='total_racism'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Total_Racist_Events) %>% mutate(quartile=factor(ntile(Total_Racist_Events, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='year_racism'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Year_Racist_events) %>% mutate(quartile=factor(ntile(Year_Racist_events, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  }  else if (int=='life_racism'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Life_Racist_events) %>% mutate(quartile=factor(ntile(Life_Racist_events, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='stress_racism'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Racist_stress) %>% mutate(quartile=factor(ntile(Racist_stress, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='kid_24h_racism'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Racism_child_24hr) %>% mutate(quartile=factor(ntile(Racism_child_24hr, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='kid_discrimination'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, Experience_Discrimination_child) %>% mutate(quartile=factor(ntile(Experience_Discrimination_child, 4), levels=1:4)) %>%
      group_by(Gene, condition, celltype, quartile) %>% summarise(mean_value=mean(value)) %>% ungroup() %>% 
      group_by(Gene, celltype, quartile) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } else if (int=='collection_infection'){
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(bulk_obj@meta.data, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, infection_status) %>% group_by(Gene, condition, celltype, infection_status) %>% 
      summarise(mean_value=mean(value)) %>%
      ungroup() %>% group_by(Gene, celltype, infection_status) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>% ungroup() %>% mutate(condition=cond)
  } 
  
  return(exp_long)
}

# compute expression scores for all combinations of conditions, celltypes, and interactions
param_grid <- expand_grid(
  condition = conditions,
  celltype = celltypes,
  interaction = interactions)

results_pathway_scores <- pmap(
  param_grid,
  \(condition, celltype, interaction){
    compute_pathway_expression(
      bulk_obj,
      condition,
      celltype,
      interaction) %>% mutate(interaction=interaction)}) 

## VISUALIZATION of IFN PATHWAY SCORES
# no interaction
(results_pathway_scores[which(param_grid$interaction=='none')] %>% bind_rows() %>% mutate(condition=factor(condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))) %>%
  filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% ggplot(., aes(x=condition, y=mean_value, fill=condition)) + 
  geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('mean exp. IFNa genes') + xlab(NULL) + 
  stat_compare_means(aes(group=Gene), method='t.test', label='p.format', paired=T, comparisons=list(c('NI_iva','IVA'), c('NI_rv','RV'))) +
  facet_wrap(~celltype, nrow=1) + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
(results_pathway_scores[which(param_grid$interaction=='none')] %>% bind_rows() %>% mutate(condition=factor(condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))) %>%
  filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% ggplot(., aes(x=condition, y=mean_value, fill=condition)) + 
  geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('mean exp. IFNy genes') + xlab(NULL) + 
  stat_compare_means(aes(group=Gene), method='t.test', label='p.format', paired=T, comparisons=list(c('NI_iva','IVA'), c('NI_rv','RV'))) +
  facet_wrap(~celltype, nrow=1) + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))))
ggsave('IFN_genes_expression_infection_violinplots_new.pdf', height=5, width=10)

# asthma
(results_pathway_scores[which(param_grid$interaction=='asthma_alb')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    mutate(Recorded_Diagnosis=recode(Recorded_Diagnosis, 'No_Diagnosis'='non_asthmatic', 'Recorded_Asthma_Diagnosis'='asthmatic')) %>%
    ggplot(., aes(x=Recorded_Diagnosis, y=delta_value, fill=Recorded_Diagnosis)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('non_asthmatic','asthmatic'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
(results_pathway_scores[which(param_grid$interaction=='asthma_alb')] %>% bind_rows() %>% 
   filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
   mutate(Recorded_Diagnosis=recode(Recorded_Diagnosis, 'No_Diagnosis'='non_asthmatic', 'Recorded_Asthma_Diagnosis'='asthmatic')) %>%
   ggplot(., aes(x=Recorded_Diagnosis, y=delta_value, fill=Recorded_Diagnosis)) + 
   geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
   stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('non_asthmatic','asthmatic'))) +
   facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_asthma_violinplots_new.pdf', height=5, width=10)

# income
(results_pathway_scores[which(param_grid$interaction=='income')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=income, y=delta_value, fill=income)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('Low','High'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='income')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=income, y=delta_value, fill=income)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('Low','High'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_income_violinplots_new.pdf', height=5, width=10)

# ACT
(results_pathway_scores[which(param_grid$interaction=='ACT')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='ACT')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_ACT_violinplots_new.pdf', height=5, width=10)

# ACE
(results_pathway_scores[which(param_grid$interaction=='ACE')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='ACE')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_ACE_violinplots_new.pdf', height=5, width=10)

# resilience
(results_pathway_scores[which(param_grid$interaction=='resilience')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='resilience')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_resilience_violinplots_new.pdf', height=5, width=10)

# social support
(results_pathway_scores[which(param_grid$interaction=='social_support')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='social_support')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_social_support_violinplots_new.pdf', height=5, width=10)

# total racism
(results_pathway_scores[which(param_grid$interaction=='total_racism')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='total_racism')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_total_racism_violinplots_new.pdf', height=5, width=10)

# year racism
(results_pathway_scores[which(param_grid$interaction=='year_racism')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='year_racism')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_year_racism_violinplots_new.pdf', height=5, width=10)

# life racism
(results_pathway_scores[which(param_grid$interaction=='life_racism')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='life_racism')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_life_racism_violinplots_new.pdf', height=5, width=10)

# stress racism
(results_pathway_scores[which(param_grid$interaction=='stress_racism')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='stress_racism')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_stress_racism_violinplots_new.pdf', height=5, width=10)

# kid 24h racism
(results_pathway_scores[which(param_grid$interaction=='kid_24h_racism')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='kid_24h_racism')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_kid_24h_racism_violinplots_new.pdf', height=5, width=10)

# kid discrimination
(results_pathway_scores[which(param_grid$interaction=='kid_discrimination')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='kid_discrimination')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=quartile, y=delta_value, fill=quartile)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('1','4'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_kid_discrimination_violinplots_new.pdf', height=5, width=10)

# infection at collection
(results_pathway_scores[which(param_grid$interaction=='collection_infection')] %>% bind_rows() %>% 
    filter(Gene %in% genes_hallmark$INTERFERON_ALPHA_RESPONSE) %>% 
    ggplot(., aes(x=infection_status, y=delta_value, fill=infection_status)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNa genes') + xlab(NULL) + 
    stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('Negative','Positive'))) +
    facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) /
  (results_pathway_scores[which(param_grid$interaction=='collection_infection')] %>% bind_rows() %>% 
     filter(Gene %in% genes_hallmark$INTERFERON_GAMMA_RESPONSE) %>% 
     ggplot(., aes(x=infection_status, y=delta_value, fill=infection_status)) + 
     geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5, width=0.3) + theme_bw() + ylab('delta mean exp. IFNy genes') + xlab(NULL) + 
     stat_compare_means(aes(group=Gene, condition), method='t.test', label='p.format', comparisons=list(c('Negative','Positive'))) +
     facet_grid(cols=vars(celltype), rows=vars(condition), scales='free') + guides(fill='none') + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))) 
ggsave('IFN_genes_expression_infection_collection_violinplots_new.pdf', height=5, width=10)

