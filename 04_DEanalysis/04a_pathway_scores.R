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
cells_seurat <- c('B','CD4-T','CD8-T','Mono','NK')
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
        
        # compute scores per hallmark pathway
        for (p in seq(length(genes_hallmark))){
          # subset genes
          subset_expression <- corrected_expression %>% as.data.frame() %>% 
            filter(rownames(.) %in% genes_hallmark[[p]]) %>% as.matrix()
          
          # scale and compute average pathway expression per individual
          subset_expression <- subset_expression %>% t() %>% scale() %>% t() %>% colMeans(na.rm=TRUE) %>% 
            as.data.frame() %>% rownames_to_column('ID')
          colnames(subset_expression)[2] <- names(genes_hallmark)[p]
          
          # group scores
          if (exists('infection_scores')){
            infection_scores <- inner_join(infection_scores, subset_expression, by=c('ID')) 
          } else {infection_scores <-  subset_expression}
        }
      
        # reformat output
        infection_scores <- infection_scores %>% separate(col=ID, into=c('ID', 'condition', 'celltype'), sep='_') %>%
          group_by(ID) %>% mutate(is_unique = n() == 1) %>% ungroup() %>% filter(is_unique==FALSE) %>% select(-is_unique) 
        infection_scores$condition <- gsub('NI', 'NI_'%&%tolower(conditions[i]), infection_scores$condition)
        
        # group infection scores
        if (exists('full_infection_scores')){
          full_infection_scores <- rbind(full_infection_scores, infection_scores) 
        } else {full_infection_scores <- infection_scores}
        rm(infection_scores)

        
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
        
        # compute scores per hallmark pathway
        for (p in seq(length(genes_hallmark))){
          # subset genes
          subset_expression <- corrected_expression %>% as.data.frame() %>% 
            filter(rownames(.) %in% genes_hallmark[[p]]) %>% as.matrix()
          
          # scale and compute average pathway expression per individual
          subset_expression <- subset_expression %>% t() %>% scale() %>% t() %>% colMeans(na.rm=TRUE) %>% 
            as.data.frame() %>% rownames_to_column('ID')
          colnames(subset_expression)[2] <- names(genes_hallmark)[p]
          
          # group scores
          if (exists('asthma_scores')){
            asthma_scores <- inner_join(asthma_scores, subset_expression, by=c('ID')) 
          } else {asthma_scores <-  subset_expression}
        }
        
        # compute paired deltas
        asthma_scores <- asthma_scores %>% separate(col=ID, into=c('ID', 'condition', 'celltype'), sep='_') %>% group_by(ID) %>% 
          filter(n_distinct(condition)>1) %>% ungroup() %>% pivot_wider(names_from=condition, values_from=-c(ID, celltype, condition)) %>%
          mutate(across(ends_with(conditions[i]), ~ . - get(sub(paste0('_', conditions[i], '$'), '_NI', cur_column())),
            .names = 'delta_{.col}')) %>% rename_with(~ sub(paste0('_', conditions[i], '$'), '', .), starts_with('delta_')) %>%
          mutate(condition=conditions[i]) %>% select('ID', 'condition', 'celltype', contains('delta_'))

        # group asthma scores
        if (exists('full_asthma_scores')){
          full_asthma_scores <- rbind(full_asthma_scores, asthma_scores) 
        } else {full_asthma_scores <- asthma_scores}
        rm(asthma_scores)
        
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
        
        # compute scores per hallmark pathway
        for (p in seq(length(genes_hallmark))){
          # subset genes
          subset_expression <- corrected_expression %>% as.data.frame() %>% 
            filter(rownames(.) %in% genes_hallmark[[p]]) %>% as.matrix()
          
          # scale and compute average pathway expression per individual
          subset_expression <- subset_expression %>% t() %>% scale() %>% t() %>% colMeans(na.rm=TRUE) %>% 
            as.data.frame() %>% rownames_to_column('ID')
          colnames(subset_expression)[2] <- names(genes_hallmark)[p]
          
          # group scores
          if (exists('income_scores')){
            income_scores <- inner_join(income_scores, subset_expression, by=c('ID')) 
          } else {income_scores <-  subset_expression}
        }
        
        # compute paired deltas
        income_scores <- income_scores %>% separate(col=ID, into=c('ID', 'condition', 'celltype'), sep='_') %>% group_by(ID) %>% 
          filter(n_distinct(condition)>1) %>% ungroup() %>% pivot_wider(names_from=condition, values_from=-c(ID, celltype, condition)) %>%
          mutate(across(ends_with(conditions[i]), ~ . - get(sub(paste0('_', conditions[i], '$'), '_NI', cur_column())),
                        .names = 'delta_{.col}')) %>% rename_with(~ sub(paste0('_', conditions[i], '$'), '', .), starts_with('delta_')) %>%
          mutate(condition=conditions[i]) %>% select('ID', 'condition', 'celltype', contains('delta_'))
        
        # group income scores
        if (exists('full_income_scores')){
          full_income_scores <- rbind(full_income_scores, income_scores)
        } else {full_income_scores <-  income_scores}
        rm(income_scores)
      }
    }
  }
}
# save files
fwrite(full_infection_scores, 'pathway_scores_infection.txt', col.names=TRUE, sep=' ')
fwrite(full_asthma_scores, 'pathway_scores_asthma.txt', col.names=TRUE, sep=' ')
fwrite(full_income_scores, 'pathway_scores_income.txt', col.names=TRUE, sep=' ')

# visualize scores for a given pathway
## INFECTION SCORES
full_infection_scores$condition <- factor(full_infection_scores$condition, levels=c('NI_iva', 'IVA', 'NI_rv', 'RV'))
full_infection_scores %>% select(condition, celltype, INTERFERON_ALPHA_RESPONSE) %>%
  ggplot(., aes(x=condition, y=INTERFERON_ALPHA_RESPONSE)) + geom_boxplot() +
  stat_compare_means(aes(group=ID), method='t.test', label='p.format', comparisons=list(c('NI_iva', 'IVA'), c('NI_rv', 'RV')), paired=T) + 
  facet_wrap(~celltype, nrow=1) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_infection_boxplots.pdf', height=3, width=10)
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_infection_boxplots.png', height=3, width=10)

full_infection_scores %>% select(condition, celltype, INTERFERON_GAMMA_RESPONSE) %>%
  ggplot(., aes(x=condition, y=INTERFERON_GAMMA_RESPONSE)) + geom_boxplot() +
  stat_compare_means(aes(group=ID), method='t.test', label='p.format', comparisons=list(c('NI_iva', 'IVA'), c('NI_rv', 'RV')), paired=T) + 
  facet_wrap(~celltype, nrow=1) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_infection_boxplots.pdf', height=3, width=10)
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_infection_boxplots.png', height=3, width=10)

## INFECTION:ASTHMA SCORES
#join asthma status
full_asthma_scores_w_mdata <- left_join(full_asthma_scores, mdata, by=c('ID'='IDs', 'condition', 'celltype')) %>%
  select(-c(batch, n, avg_mt, age, gender, income, albuterol))
full_asthma_scores_w_mdata$asthma <- factor(full_asthma_scores_w_mdata$asthma, levels=c('No', 'Yes'))

full_asthma_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_ALPHA_RESPONSE,asthma) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_ALPHA_RESPONSE, fill=asthma)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype, nrow=1) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_asthma_boxplots.pdf', height=4, width=10)

full_asthma_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_GAMMA_RESPONSE, asthma) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_GAMMA_RESPONSE, fill=asthma)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype, nrow=1) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_asthma_boxplots.pdf', height=4, width=10)

a <- full_asthma_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_ALPHA_RESPONSE, 
                                           delta_INTERFERON_GAMMA_RESPONSE, asthma) %>%
  filter(celltype=='CD8-T', condition=='RV') %>% pivot_longer(cols=contains('delta'))
a$name <- gsub('delta_INTERFERON_ALPHA_RESPONSE', 'delta_IFNa', a$name)
a$name <- gsub('delta_INTERFERON_GAMMA_RESPONSE', 'delta_IFNy', a$name)
ggplot(a, aes(x=name, y=value, fill=asthma)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype, nrow=1) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSEs_CD8-T_RV_ssGSEA_asthma_boxplots.pdf', height=4, width=5)
ggsave('INTERFERON_RESPONSEs_CD8-T_RV_ssGSEA_asthma_boxplots.png', height=4, width=5)


## INFECTION:INCOME SCORES
#join income status
full_income_scores_w_mdata <- left_join(full_income_scores, mdata, by=c('ID'='IDs', 'condition', 'celltype')) %>%
  select(-c(batch, n, avg_mt, age, gender, asthma, albuterol))
full_income_scores_w_mdata$income <- ifelse(full_income_scores_w_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Lower', 'Higher')
full_income_scores_w_mdata$income <- factor(full_income_scores_w_mdata$income, levels=c('Lower','Higher'))

full_income_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_ALPHA_RESPONSE, income) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_ALPHA_RESPONSE, fill=income)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_ALPHA_RESPONSE_ssGSEA_income_boxplots.pdf', height=4, width=9)

full_income_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_GAMMA_RESPONSE, income) %>%
  ggplot(., aes(x=condition, y=delta_INTERFERON_GAMMA_RESPONSE, fill=income)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_GAMMA_RESPONSE_ssGSEA_income_boxplots.pdf', height=4, width=9)

a <- full_income_scores_w_mdata %>% select(condition, celltype, delta_INTERFERON_ALPHA_RESPONSE, 
                                           delta_INTERFERON_GAMMA_RESPONSE, income) %>%
  filter(celltype %in% c('CD8-T','CD4-T','Mono') & condition=='IVA' |
         celltype %in% c('CD4-T','Mono') & condition=='RV') %>% pivot_longer(cols=contains('delta'))
a$name <- gsub('delta_INTERFERON_ALPHA_RESPONSE', 'delta_IFNa', a$name)
a$name <- gsub('delta_INTERFERON_GAMMA_RESPONSE', 'delta_IFNy', a$name)

ggplot(a, aes(x=name, y=value, fill=income)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_grid(cols=vars(celltype), rows=vars(condition)) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSEs_ssGSEA_income_grid.boxplots.pdf', height=4, width=7)
ggsave('INTERFERON_RESPONSEs_ssGSEA_income_grid.boxplots.png', height=4, width=7)

ggplot(a, aes(x=name, y=value, fill=income)) + geom_boxplot() +
  stat_compare_means(method='t.test', label='p.format') + 
  facet_wrap(~celltype+condition, nrow=1) + theme_bw() + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave('INTERFERON_RESPONSEs_ssGSEA_income_wrap.boxplots.pdf', height=4, width=9)
ggsave('INTERFERON_RESPONSEs_ssGSEA_income_wrap.boxplots.png', height=4, width=9)
