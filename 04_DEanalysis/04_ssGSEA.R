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

### important things to load
# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# get human hallmark gene sets
ifn_genes <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)

# filter for IFN modules 
ifn_genes <- ifn_genes[grep('INTERFERON', names(ifn_genes), value=T)]
names(ifn_genes) <- c('IFNa','IFNy')

# load pseudobulk seurat object
bulk_obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds') 
bulk_obj@meta.data$condition <- factor(bulk_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# merge metadata
mdata <- bulk_obj@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
bulk_obj@meta.data <- mdata
###

for (cond in c('RV')){
  for (ctype in c('T-CD8')){
    
    # subset seurat object
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$albuterol <- factor(sub_mdata$albuterol, levels=c('No', 'Yes'))
    sub_mdata$asthma <- factor(sub_mdata$asthma, levels=c('No', 'Yes'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_asthma_alb_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # subset metadata
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
    
    # extract and subset count matrix
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% v_samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # first try: scores with just log transformed counts, don't adjust for anything yet
    v <- voom(count, plot=F)
    first_ssgsea_scores <- gsva(ssgseaParam(v$E, ifn_genes))

    # second try: scores with voom-adjusted counts adjusted for covariates using removeBatchEffect()
    v <- voom(count, plot=F)
    sub_mdata$gender_num <- as.numeric(sub_mdata$gender)
    sub_mdata$albuterol_num <- as.numeric(sub_mdata$albuterol) 
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition+asthma+condition:asthma, data=sub_mdata))
    second_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # third try: same thing as before, but now voom takes a design matrix to compute the weights
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition+asthma+condition:asthma, data=sub_mdata))
    third_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))  
    
    # fourth try: fit a linear model on voom-adjusted reads, extract residuals and add the intercept back
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    intercept <- fit$coefficients[,'(Intercept)']
    corrected_expression <- residuals + intercept
    fourth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # compute paired deltas
    first_ssgsea_scores <- first_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='1')
    
    second_ssgsea_scores <- second_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='2')
    
    third_ssgsea_scores <- third_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='3')
    
    fourth_ssgsea_scores <- fourth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='4')
    
    # concatenate
    ssgsea_scores <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores) %>% 
      mutate(condition=cond)
  }
}

# merge w metadata to retrieve asthma status
ssgsea_scores.wasthma <- ssgsea_scores %>% inner_join(sample_m, by=c('ID')) %>% 
  select(ID, celltype, condition, deltaIFNa, deltaIFNy, asthma, method) %>% 
  pivot_longer(c(deltaIFNa, deltaIFNy), names_to='IFN', values_to='score') %>% 
  inner_join(mdata, by=c('ID'='IDs', 'celltype', 'condition', 'asthma')) %>% 
  select(ID, celltype, condition, method, IFN, score, asthma)

# read asthma manual scores
manual_asthma_scores <- fread('IFNscores_asthma_manual.txt') %>% rename(ID=donor) %>% mutate(method='manual') %>% 
  select(ID, celltype, condition, method, IFN, score, asthma)
manual_asthma_scores$IFN <- gsub('IFNa', 'deltaIFNa', manual_asthma_scores$IFN)
manual_asthma_scores$IFN <- gsub('IFNy', 'deltaIFNy', manual_asthma_scores$IFN)
manual_asthma_scores <- manual_asthma_scores %>% filter(celltype=='T-CD8', condition=='RV')
joint_asthma <- rbind(ssgsea_scores.wasthma, manual_asthma_scores) 
fwrite(joint_asthma, 'joint_asthma_scores.txt', sep=' ')
joint_asthma$method <- factor(joint_asthma$method, levels=c('1', '2', '3', '4', 'manual'))

# boxplot of scores split by asthma status
ggplot(joint_asthma, aes(x=condition, y=score, fill=asthma)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ggtitle('RV - T-CD8')
ggsave('IFNscores_asthma_boxplots.pdf', height=4, width=9)

# correlate scores
joint_asthma <- joint_asthma %>% pivot_wider(names_from=method, values_from=score)
for (i in (6:ncol(joint_asthma))) {
  for (j in (6:ncol(joint_asthma))) {
    
    corr <- joint_asthma %>% group_by(IFN) %>% do({
      test <- cor.test(as.matrix(.data[,i]), as.matrix(.data[,j]), method='spearman') 
      test <- broom::tidy(test) %>% mutate(method1=colnames(joint_asthma)[i], method2=colnames(joint_asthma)[j])
    })
    
    if (exists('final.cor')){
      final.cor <- rbind(final.cor, corr)
    } else (final.cor <- corr)
  }
}
final.cor$method1 <- factor(final.cor$method1, levels=c('1', '2', '3', '4', 'manual'))
final.cor$method2 <- factor(final.cor$method2, levels=c('1', '2', '3', '4', 'manual'))
ggplot(final.cor, aes(x=method1, y=method2, fill=estimate)) + geom_tile() + facet_wrap(~IFN) +
  theme_bw()
ggsave('IFNscores_asthma_corr_bwtnmethods_heatmap.pdf', height=3, width=6)
rm(ssgsea_scores)

#### INCOME ####
for (cond in c('IVA')){
  for (ctype in c('Mono', 'T-CD4')){
    
    # subset seurat object
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$income <- na_if(sub_mdata$income, '')
    sub_mdata$income <- ifelse(sub_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Lower', 'Higher')
    sub_mdata$income <- factor(sub_mdata$income, levels=c('Lower','Higher'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_income_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # subset metadata
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
    
    # extract and subset count matrix
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% v_samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # first try: scores with just log transformed counts, don't adjust for anything yet
    v <- voom(count, plot=F)
    first_ssgsea_scores <- gsva(ssgseaParam(v$E, ifn_genes))
    
    # second try: scores with voom-adjusted counts adjusted for covariates using removeBatchEffect()
    v <- voom(count, plot=F)
    sub_mdata$gender_num <- as.numeric(sub_mdata$gender)
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition+income+condition:income, data=sub_mdata))
    second_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # third try: same thing as before, but now voom takes a design matrix to compute the weights
    design <- model.matrix(~batch+age+gender+n+avg_mt, data=sub_mdata)
    v <- voom(count, design, plot=F)
    sub_mdata$gender_num <- as.numeric(sub_mdata$gender)
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition+income+condition:income, data=sub_mdata))
    third_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))  

    # fourth try: fit a linear model on voom-adjusted reads, extract residuals and add the intercept back
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    intercept <- fit$coefficients[,'(Intercept)']
    corrected_expression <- residuals + intercept
    fourth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # compute paired deltas
    first_ssgsea_scores <- first_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
              deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='1')
    
    second_ssgsea_scores <- second_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
              deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='2')
    
    third_ssgsea_scores <- third_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
              deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='3')
    
    fourth_ssgsea_scores <- fourth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      reframe(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
              deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% mutate(method='4')
    
    # concatenate
    if (exists('ssgsea_scores')){
    tmp <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores) %>% mutate(condition=cond)
    ssgsea_scores <- rbind(ssgsea_scores, tmp)
    } else {
      ssgsea_scores <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores) %>% mutate(condition=cond)
    }
  }
}
# merge w metadata to retrieve income status
ssgsea_scores.wincome <- ssgsea_scores %>% inner_join(sample_m, by=c('ID')) %>% 
  select(ID, celltype, condition, deltaIFNa, deltaIFNy, income, method) %>% 
  pivot_longer(c(deltaIFNa, deltaIFNy), names_to='IFN', values_to='score') %>% 
  inner_join(mdata, by=c('ID'='IDs', 'celltype', 'condition', 'income')) %>% 
  select(ID, celltype, condition, method, IFN, score, income)
ssgsea_scores.wincome$income <- na_if(ssgsea_scores.wincome$income, '')
ssgsea_scores.wincome$income <- ifelse(ssgsea_scores.wincome$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                                       'Lower', 'Higher')
ssgsea_scores.wincome$income <- factor(ssgsea_scores.wincome$income, levels=c('Lower','Higher'))

# read income manual scores
manual_income_scores <- fread('IFNscores_income_manual.txt') %>% rename(ID=donor) %>% mutate(method='manual') %>% 
  select(ID, celltype, condition, method, IFN, score, income) %>% filter(celltype %in% c('T-CD4', 'Mono'), condition=='IVA')
manual_income_scores$IFN <- gsub('IFNa', 'deltaIFNa', manual_income_scores$IFN)
manual_income_scores$IFN <- gsub('IFNy', 'deltaIFNy', manual_income_scores$IFN)
joint_income <- rbind(ssgsea_scores.wincome, manual_income_scores) 
fwrite(joint_income, 'joint_income_scores.txt', sep=' ')

# boxplot of scores split by income status
joint_income %>% filter(celltype=='Mono') %>% ggplot(., aes(x=condition, y=score, fill=income)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = income), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ggtitle('IVA - Mono')
ggsave('IFNscores_income_boxplots_mono.pdf', height=4, width=9)

joint_income %>% filter(celltype=='T-CD4') %>% ggplot(., aes(x=condition, y=score, fill=income)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = income), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ggtitle('IVA - T-CD4')
ggsave('IFNscores_income_boxplots_tcd4.pdf', height=4, width=9)

# correlate scores
joint_income <- joint_income %>% pivot_wider(names_from=method, values_from=score)
for (i in (6:ncol(joint_income))) {
  for (j in (6:ncol(joint_income))) {
    
    corr <- joint_income %>% group_by(celltype, IFN) %>% do({
      test <- cor.test(as.matrix(.data[,i]), as.matrix(.data[,j]), method='spearman') 
      test <- broom::tidy(test) %>% mutate(method1=colnames(joint_income)[i], method2=colnames(joint_income)[j])
    })
    
    if (exists('final.cor.in')){
      final.cor.in <- rbind(final.cor.in, corr)
    } else (final.cor.in <- corr)
  }
}
final.cor.in$method1 <- factor(final.cor.in$method1, levels=c('1', '2', '3', '4', 'manual'))
final.cor.in$method2 <- factor(final.cor.in$method2, levels=c('1', '2', '3', '4', 'manual'))
ggplot(final.cor.in, aes(x=method1, y=method2, fill=estimate)) + geom_tile() + facet_grid(cols=vars(IFN), rows=vars(celltype)) +
  theme_bw()
ggsave('IFNscores_income_corr_bwtnmethods_heatmap.pdf', height=5, width=7)


### official ssGSEA replace everything but this eventually
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
