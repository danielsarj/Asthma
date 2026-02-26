library(Seurat)
library(msigdbr)
library(data.table)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(edgeR)
library(limma)
library(purrr)
library(ggpmisc)
"%&%" <- function(a,b) paste(a,b, sep = "")
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


compute_pathway_scores <- function(seurat_obj, hallmark, cond, ctype, int) {
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
  
  # compute scores per hallmark pathway
  for (p in seq(length(hallmark))){
    # subset genes
    subset_expression <- corrected_expression %>% as.data.frame() %>% filter(rownames(.) %in% hallmark[[p]]) %>% as.matrix()
    
    # scale and compute average pathway expression per individual
    subset_expression <- subset_expression %>% t() %>% scale() %>% t() %>% colMeans(na.rm=TRUE) %>% as.data.frame() %>% rownames_to_column('ID')
    colnames(subset_expression)[2] <- names(hallmark)[p]
    
    # group scores
    if (exists('tmp_scores')){
      tmp_scores <- inner_join(tmp_scores, subset_expression, by=c('ID')) 
    } else {tmp_scores <-  subset_expression}
  }
  
  # reformat output
  path_scores <- tmp_scores %>% separate(col=ID, into=c('ID', 'condition', 'celltype'), sep='_') %>%
    group_by(ID) %>% mutate(is_unique = n() == 1, interaction=int) %>% ungroup() %>% filter(is_unique==FALSE) %>% select(-is_unique) %>%
    relocate(ID, condition, celltype, interaction)
  path_scores$condition <- gsub('NI', 'NI_'%&%tolower(cond), path_scores$condition)
  
  return(path_scores)
}

# compute pathway scores for all combinations of conditions, celltypes, and interactions
param_grid <- expand_grid(
  condition = conditions,
  celltype = celltypes,
  interaction = interactions)
results_pathway_scores <- pmap(
  param_grid,
  \(condition, celltype, interaction){
    compute_pathway_scores(
      bulk_obj,
      genes_hallmark,
      condition,
      celltype,
      interaction)}) %>% bind_rows()


## VISUALIZATION of IFN PATHWAY SCORES
# no interaction
infection_scores <- results_pathway_scores %>% filter(interaction=='none') %>% 
  select(condition, celltype, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score')
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))

ggplot(infection_scores, aes(x=condition, y=score)) + geom_boxplot()+ ylab('pathway score') +
  stat_compare_means(aes(group=ID), method='t.test', label='p.signif', comparisons=list(c('NI_iva', 'IVA'), c('NI_rv', 'RV')), paired=T) + 
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_infection_boxplots_new.png', height=5, width=10)

# interactions
IVA_scores <- results_pathway_scores %>% 
  filter(!interaction %in% c('none') & condition %in% c('NI_iva', 'IVA'))
IVA_scores$condition <- gsub('NI_iva', 'NI', IVA_scores$condition)
IVA_scores <- IVA_scores %>% group_by(condition, celltype, interaction) %>% 
  pivot_wider(names_from=condition, values_from=-c(ID, celltype, condition, interaction)) %>%
  ungroup() %>% mutate(across(ends_with('IVA'), ~ . - get(sub(paste0('_', 'IVA', '$'), '_NI', cur_column())),
                              .names = 'delta_{.col}')) %>% rename_with(~ sub(paste0('_', 'IVA', '$'), '', .), starts_with('delta_')) %>%
  mutate(condition='IVA') %>% select('ID', 'condition', 'celltype', 'interaction', contains('delta_'))

RV_scores <- results_pathway_scores %>% 
  filter(!interaction %in% c('none') & condition %in% c('NI_rv', 'RV'))
RV_scores$condition <- gsub('NI_rv', 'NI', RV_scores$condition)
RV_scores <- RV_scores %>% group_by(condition, celltype, interaction) %>% 
  pivot_wider(names_from=condition, values_from=-c(ID, celltype, condition, interaction)) %>%
  ungroup() %>% mutate(across(ends_with('RV'), ~ . - get(sub(paste0('_', 'RV', '$'), '_NI', cur_column())),
                .names = 'delta_{.col}')) %>% rename_with(~ sub(paste0('_', 'RV', '$'), '', .), starts_with('delta_')) %>%
  mutate(condition='RV') %>% select('ID', 'condition', 'celltype', 'interaction', contains('delta_'))

full_scores <- rbind(IVA_scores, RV_scores) %>% left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'condition', 'celltype'))

# asthma
full_scores %>% filter(interaction=='asthma_alb') %>% 
  mutate(Recorded_Diagnosis=if_else(Recorded_Diagnosis=='Recorded_Asthma_Diagnosis', 'asthmatic', 'non_asthmatic')) %>%
  select(condition, celltype, Recorded_Diagnosis, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=condition, y=score, fill=Recorded_Diagnosis)) + geom_boxplot(alpha=0.5) +
  stat_compare_means(method='t.test', label='p.signif') + ylab('delta pathway score (Inf-NI)') + xlab(NULL) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_asthma_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='asthma_alb') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  mutate(Recorded_Diagnosis=if_else(Recorded_Diagnosis=='Recorded_Asthma_Diagnosis', 'asthmatic', 'non_asthmatic')) %>%
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score')
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=Recorded_Diagnosis)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=Recorded_Diagnosis), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_asthma_boxplots_new.png', height=5, width=10)

# income
full_scores %>% filter(interaction=='income') %>% 
  select(condition, celltype, income, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=condition, y=score, fill=income)) + geom_boxplot(alpha=0.5) +
  stat_compare_means(method='t.test', label='p.signif') + ylab('delta pathway score (Inf-NI)') + xlab(NULL) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_income_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='income') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score')
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=income)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=income), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_income_boxplots_new.png', height=5, width=10)

# ACT
full_scores %>% filter(interaction=='ACT') %>% 
  select(condition, celltype, ACT_score, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=ACT_score, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('ACT score') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_ACT_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='ACT') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(ACT_score, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_ACT_boxplots_new.png', height=5, width=10)

# ACE
full_scores %>% filter(interaction=='ACE') %>% 
  select(condition, celltype, ACE_result, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=ACE_result, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('ACE result') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_ACE_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='ACE') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(ACE_result, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_ACE_boxplots_new.png', height=5, width=10)

# resilience
full_scores %>% filter(interaction=='resilience') %>% 
  select(condition, celltype, Parent_Resilience_Score, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Parent_Resilience_Score, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('resilience') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_resilience_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='resilience') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Parent_Resilience_Score, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_resilience_boxplots_new.png', height=5, width=10)

# social support
full_scores %>% filter(interaction=='social_support') %>% 
  select(condition, celltype, Parents_Score_Avg, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Parents_Score_Avg, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('social support') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_social_support_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='social_support') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Parents_Score_Avg, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_social_support_boxplots_new.png', height=5, width=10)

# total racism
full_scores %>% filter(interaction=='total_racism') %>% 
  select(condition, celltype, Total_Racist_Events, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Total_Racist_Events, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('total racism') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_total_racism_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='total_racism') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Total_Racist_Events, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_total_racism_boxplots_new.png', height=5, width=10)

# year racism
full_scores %>% filter(interaction=='year_racism') %>% 
  select(condition, celltype, Year_Racist_events, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Year_Racist_events, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('year racism') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_year_racism_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='year_racism') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Year_Racist_events, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_year_racism_boxplots_new.png', height=5, width=10)

# life racism
full_scores %>% filter(interaction=='life_racism') %>% 
  select(condition, celltype, Life_Racist_events, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Life_Racist_events, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('life racism') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_life_racism_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='life_racism') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Life_Racist_events, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_life_racism_boxplots_new.png', height=5, width=10)

# stress racism
full_scores %>% filter(interaction=='stress_racism') %>% 
  select(condition, celltype, Racist_stress, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Racist_stress, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('stress racism') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_stress_racism_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='stress_racism') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Racist_stress, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_stress_racism_boxplots_new.png', height=5, width=10)

# kid24h racism
full_scores %>% filter(interaction=='kid_24h_racism') %>% 
  select(condition, celltype, Racism_child_24hr, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Racism_child_24hr, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('kid24h racism') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_kid_24h_racism_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='kid_24h_racism') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Racism_child_24hr, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_stress_kid_24h_racism_boxplots_new.png', height=5, width=10)

# kid discrimination
full_scores %>% filter(interaction=='kid_discrimination') %>% 
  select(condition, celltype, Experience_Discrimination_child, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=Experience_Discrimination_child, y=score, color=condition)) + geom_point() + ylab('delta pathway score (Inf-NI)') + xlab('kid discrimination') +
  geom_smooth(method='lm', se=TRUE) + stat_poly_eq(aes(label=paste(after_stat(p.value.label))), formula=y~x, parse=TRUE) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_kid_discrimination_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='kid_discrimination') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>%
  mutate(quartile=factor(ntile(Experience_Discrimination_child, 4), levels=1:4)) %>% filter(quartile==1 | quartile==4)
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=quartile)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=quartile), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_stress_kkid_discrimination_boxplots_new.png', height=5, width=10)

# infection at collection
full_scores %>% filter(interaction=='collection_infection') %>% 
  select(condition, celltype, infection_status, delta_INTERFERON_ALPHA_RESPONSE, delta_INTERFERON_GAMMA_RESPONSE) %>% 
  rename(IFNa=delta_INTERFERON_ALPHA_RESPONSE, IFNy=delta_INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score') %>% 
  ggplot(., aes(x=condition, y=score, fill=infection_status)) + geom_boxplot(alpha=0.5) +
  stat_compare_means(method='t.test', label='p.signif') + ylab('delta pathway score (Inf-NI)') + xlab(NULL) +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_delta_scores_collection_infection_boxplots_new.png', height=5, width=10)

infection_scores <- results_pathway_scores %>% filter(interaction=='collection_infection') %>%
  select(ID, condition, celltype, interaction, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE) %>%
  mutate(main_condition=recode(condition, 'NI_iva'='NI', 'NI_rv'='NI')) %>% 
  left_join(bulk_obj@meta.data, by=c('ID'='IDs', 'main_condition'='condition', 'celltype')) %>% 
  rename(IFNa=INTERFERON_ALPHA_RESPONSE, IFNy=INTERFERON_GAMMA_RESPONSE) %>% 
  pivot_longer(cols=c(IFNa, IFNy), names_to='pathway', values_to='score')
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
ggplot(infection_scores, aes(x=condition, y=score, fill=infection_status)) + geom_boxplot(alpha=0.5) + ylab('pathway score') +
  stat_compare_means(aes(group=infection_status), method='t.test', label='p.signif') +
  facet_grid(cols=vars(celltype), rows=vars(pathway)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSE_pathway_scores_collection_infection_boxplots_new.png', height=5, width=10)
