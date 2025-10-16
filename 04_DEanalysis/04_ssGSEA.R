library(GSVA)
library(Seurat)
library(msigdbr)
library(data.table)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(edgeR)
library(limma)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('none','asthma','income')

# get human hallmark gene sets
genes_hallmark <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)
names(genes_hallmark) <- gsub('HALLMARK_', '', names(genes_hallmark))

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# load pseudobulk seurat object
bulk_obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds') 
bulk_obj@meta.data$condition <- factor(bulk_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# merge metadata
mdata <- bulk_obj@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
bulk_obj@meta.data <- mdata

for (int in interactions){
  for (i in 1:length(conditions)){
    for (ctype in cells_seurat){
      print(c(int, conditions[i], ctype))
      
      # LOOK INTO INFECTION SCORES
      if (int == 'none'){
        # subset seurat object
        tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==conditions[i]))
        sub_mdata <- tmp@meta.data
        sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))

        # get samples/genes from voom expression dataframe (easier to filter stuff)
        v <- fread('../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_voom_expression.txt') %>% 
          column_to_rownames('Gene')
        v_samples <- colnames(v)
        v_genes <- rownames(v)
        rm(v)
        
        # subset metadata
        sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
        
        # extract and subset count matrix
        count <- tmp@assays$RNA$counts
        count <- count[rownames(count) %in% v_genes,]
        count <- count[,colnames(count) %in% v_samples]
        count <- DGEList(counts=count) %>% calcNormFactors()
        
        # create design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt, data=sub_mdata)
        
        # fit linear model
        voom_obj <- voom(count, design, plot=F)
        l_model <- lmFit(voom_obj, design=design) %>% eBayes()
        
        # extract residuals and add the intercept
        residuals <- residuals.MArrayLM(l_model, voom_obj)
        intercept <- l_model$coefficients[,'(Intercept)']
        corrected_expression <- residuals + intercept
        
        # compute scores
        tmp_scores <- gsva(ssgseaParam(corrected_expression, genes_hallmark))
        
        # reformat output
        tmp_scores <- tmp_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
          separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID) %>% filter(n_distinct(condition)>1)
        tmp_scores$condition <- gsub('NI', 'NI_'%&%tolower(conditions[i]), tmp_scores$condition)
        
        # group infection scores
        if (exists('infection_scores')){
          infection_scores <- rbind(infection_scores, tmp_scores) 
        } else {infection_scores <-  tmp_scores}


      # LOOK INTO INFECTION:ASTHMA SCORES
      } else if (int == 'asthma'){
        # subset seurat object
        tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==conditions[i]))
        sub_mdata <- tmp@meta.data
        sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
        sub_mdata$albuterol <- factor(sub_mdata$albuterol, levels=c('No', 'Yes'))

        # get samples/genes from voom expression dataframe (easier to filter stuff)
        v <- fread('../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_asthma_alb_voom_expression.txt') %>% 
          column_to_rownames('Gene')
        v_samples <- colnames(v)
        v_genes <- rownames(v)
        rm(v)
        
        # subset metadata
        sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
        
        # extract and subset count matrix
        count <- tmp@assays$RNA$counts
        count <- count[rownames(count) %in% v_genes,]
        count <- count[,colnames(count) %in% v_samples]
        count <- DGEList(counts=count) %>% calcNormFactors()

        # create design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol, data=sub_mdata)
        
        # fit linear model
        voom_obj <- voom(count, design, plot=F)
        l_model <- lmFit(voom_obj, design=design) %>% eBayes()
        
        # extract residuals and add the intercept
        residuals <- residuals.MArrayLM(l_model, voom_obj)
        intercept <- l_model$coefficients[,'(Intercept)']
        corrected_expression <- residuals + intercept
        
        # compute scores
        tmp_scores <- gsva(ssgseaParam(corrected_expression, genes_hallmark))
        
        # compute paired deltas
        tmp_scores <- tmp_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
          separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID) %>% filter(n_distinct(condition)>1) %>% 
          ungroup() %>% pivot_wider(names_from=condition, values_from=-c(ID, celltype, condition)) %>%
          mutate(across(ends_with(conditions[i]), ~ . - get(sub(paste0('_', conditions[i], '$'), '_NI', cur_column())),
            .names = 'delta_{.col}')) %>% rename_with(~ sub(paste0('_', conditions[i], '$'), '', .), starts_with('delta_')) %>%
          mutate(condition=conditions[i]) %>% select('ID', 'condition', 'celltype', contains('delta_'))
        
        # group asthma scores
        if (exists('asthma_scores')){
          asthma_scores <- rbind(asthma_scores, tmp_scores)
        } else {asthma_scores <-  tmp_scores}
        
        
        # LOOK INTO INFECTION:INCOME SCORES
      } else {
        # subset seurat object
        tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==conditions[i]))
        sub_mdata <- tmp@meta.data
        sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))

        # get samples/genes from voom expression dataframe (easier to filter stuff)
        v <- fread('../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_income_voom_expression.txt') %>% 
          column_to_rownames('Gene')
        v_samples <- colnames(v)
        v_genes <- rownames(v)
        rm(v)
        
        # subset metadata
        sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
        
        # extract and subset count matrix
        count <- tmp@assays$RNA$counts
        count <- count[rownames(count) %in% v_genes,]
        count <- count[,colnames(count) %in% v_samples]
        count <- DGEList(counts=count) %>% calcNormFactors()
        
        # create design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt, data=sub_mdata)
        
        # fit linear model
        voom_obj <- voom(count, design, plot=F)
        l_model <- lmFit(voom_obj, design=design) %>% eBayes()
        
        # extract residuals and add the intercept
        residuals <- residuals.MArrayLM(l_model, voom_obj)
        intercept <- l_model$coefficients[,'(Intercept)']
        corrected_expression <- residuals + intercept
        
        # compute scores
        tmp_scores <- gsva(ssgseaParam(corrected_expression, genes_hallmark))
        
        # compute paired deltas
        tmp_scores <- tmp_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
          separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID) %>% filter(n_distinct(condition)>1) %>% 
          ungroup() %>% pivot_wider(names_from=condition, values_from=-c(ID, celltype, condition)) %>%
          mutate(across(ends_with(conditions[i]), ~ . - get(sub(paste0('_', conditions[i], '$'), '_NI', cur_column())),
                        .names = 'delta_{.col}')) %>% rename_with(~ sub(paste0('_', conditions[i], '$'), '', .), starts_with('delta_')) %>%
          mutate(condition=conditions[i]) %>% select('ID', 'condition', 'celltype', contains('delta_'))
        
        # group asthma scores
        if (exists('income_scores')){
          income_scores <- rbind(income_scores, tmp_scores)
        } else {income_scores <-  tmp_scores}
      }
    }
  }
}

# visualize scores for a given pathway

## INFECTION SCORES
infection_scores$condition <- factor(infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
infection_scores %>% select(condition, celltype, INTERFERON_ALPHA_RESPONSE) %>%
  ggplot(., aes(x=condition, y=INTERFERON_ALPHA_RESPONSE)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format', comparisons=list(c('NI_iva', 'IVA'), c('NI_rv', 'RV'))) + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_infection_boxplots.pdf', height=4, width=9)

infection_scores %>% select(condition, celltype, INTERFERON_GAMMA_RESPONSE) %>%
  ggplot(., aes(x=condition, y=INTERFERON_GAMMA_RESPONSE)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format', comparisons=list(c('NI_iva', 'IVA'), c('NI_rv', 'RV'))) + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_infection_boxplots.pdf', height=4, width=9)

## INFECTION:ASTHMA SCORES
#join asthma status
asthma_scores_w_mdata <- left_join(asthma_scores, mdata, by=c('ID'='IDs', 'condition', 'celltype')) %>%
  select(-c(batch, n, avg_mt, age, gender, income, albuterol))
asthma_scores_w_mdata$asthma <- factor(asthma_scores_w_mdata$asthma, levels=c('No', 'Yes'))

asthma_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_ALPHA_RESPONSE, asthma) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_ALPHA_RESPONSE, fill=asthma)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_asthma_boxplots.pdf', height=4, width=9)

asthma_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_GAMMA_RESPONSE, asthma) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_GAMMA_RESPONSE, fill=asthma)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_asthma_boxplots.pdf', height=4, width=9)

## INFECTION:INCOME SCORES
#join income status
income_scores_w_mdata <- left_join(income_scores, mdata, by=c('ID'='IDs', 'condition', 'celltype')) %>%
  select(-c(batch, n, avg_mt, age, gender, asthma, albuterol))
income_scores_w_mdata$income <- ifelse(income_scores_w_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Lower', 'Higher')
income_scores_w_mdata$income <- factor(income_scores_w_mdata$income, levels=c('Lower','Higher'))

income_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_ALPHA_RESPONSE, income) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_ALPHA_RESPONSE, fill=income)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_income_boxplots.pdf', height=4, width=9)

income_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_GAMMA_RESPONSE, income) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_GAMMA_RESPONSE, fill=income)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_income_boxplots.pdf', height=4, width=9)
