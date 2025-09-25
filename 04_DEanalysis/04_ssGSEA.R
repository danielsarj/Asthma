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
    design <- model.matrix(~1, data=count$samples)
    v <- voom(count, design, plot=F)
    v <- v$E
    first_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))

    # second try: scores with design matrix containing all terms but asthma
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    second_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # third try: same design, but now using removeBatchEffect function
    sub_mdata$condition_num <- as.numeric(sub_mdata$condition)
    sub_mdata$gender_num <- as.numeric(sub_mdata$gender)
    sub_mdata$albuterol_num <- as.numeric(sub_mdata$albuterol) 
    sub_mdata$asthma_num <- as.numeric(sub_mdata$asthma)
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num','condition_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~asthma_num+condition_num:asthma_num, data=sub_mdata))
    third_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # fourth try: pull from the second try, but add the residuals back
    v <- voom(count, design, plot=F)
    v <- v$E
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    fourth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # fifth try: include asthma, but not the interaction, in the design
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition+asthma, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    fifth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # fifth try: same design, but now using removeBatchEffect function
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num','condition_num','asthma_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition_num:asthma_num, data=sub_mdata))
    sixth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # seventh try: include everything including the interaction in the design
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition*asthma, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    seventh_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # eighth try: same design, but now using removeBatchEffect function
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num','condition_num','asthma_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition_num:asthma_num, data=sub_mdata))
    eighth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # ninth try: a simpler design model
    ninth <- model.matrix(~0+batch, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    ninth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # tenth try: doing whatever it is when adding back batch effects from the design above?
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    tenth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # eleventh try: doing whatever it is when adding back effects from the design above, but now with everything but asthma
    design <- model.matrix(~0+age+gender+n+avg_mt+albuterol+condition, data=sub_mdata)
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    eleventh_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # twelfth try: doing whatever it is when adding back effects from the design above, but now with everything
    design <- model.matrix(~0+age+gender+n+avg_mt+albuterol+condition+asthma, data=sub_mdata)
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    twelfth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # thirteenth: the simplest model ever
    v <- voom(count, plot=F)
    v <- v$E
    thirteenth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # compute paired deltas
    first_ssgsea_scores <- first_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='1')
    
    second_ssgsea_scores <- second_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='2')
    
    third_ssgsea_scores <- third_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='3')
    
    fourth_ssgsea_scores <- fourth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='4')
    
    fifth_ssgsea_scores <- fifth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='5')
    
    sixth_ssgsea_scores <- sixth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='6')
    
    seventh_ssgsea_scores <- seventh_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='7')
    
    eighth_ssgsea_scores <- eighth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='8')
    
    ninth_ssgsea_scores <- ninth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='9')
    
    tenth_ssgsea_scores <- tenth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='10')
    
    eleventh_ssgsea_scores <- eleventh_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='11')
    
    twelfth_ssgsea_scores <- twelfth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='12')
    
    thirteenth_ssgsea_scores <- thirteenth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='13')
    
    # concatenate
    ssgsea_scores <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores, 
                           fifth_ssgsea_scores, sixth_ssgsea_scores, seventh_ssgsea_scores, eighth_ssgsea_scores,
                           ninth_ssgsea_scores, tenth_ssgsea_scores, eleventh_ssgsea_scores, twelfth_ssgsea_scores,
                           thirteenth_ssgsea_scores) %>% mutate(condition=cond)
  }
}

# merge w metadata to retrieve asthma status
ssgsea_scores.wasthma <- ssgsea_scores %>% inner_join(sample_m, by=c('ID')) %>% 
  select(ID, celltype, condition, deltaIFNa, deltaIFNy, asthma, method) %>% 
  pivot_longer(c(deltaIFNa, deltaIFNy), names_to='IFN', values_to='score') %>% 
  inner_join(mdata, by=c('ID'='IDs', 'celltype', 'condition', 'asthma'))
ssgsea_scores.wasthma$method <- factor(ssgsea_scores.wasthma$method, levels=c('1', '2', '3', '4', '5', '6', '7','8', '9',
                                                                              '10', '11', '12', '13'))

ggplot(ssgsea_scores.wasthma, aes(x=condition, y=score, fill=asthma)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ggtitle('RV - T-CD8')
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
    design <- model.matrix(~1, data=count$samples)
    v <- voom(count, design, plot=F)
    v <- v$E
    first_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # second try: scores with design matrix containing all terms but income
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    second_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # third try: same design, but now using removeBatchEffect function
    sub_mdata$condition_num <- as.numeric(sub_mdata$condition)
    sub_mdata$gender_num <- as.numeric(sub_mdata$gender)
    sub_mdata$income_num <- as.numeric(sub_mdata$income) 
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','condition_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~income_num+condition_num:income_num, data=sub_mdata))
    third_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # fourth try: pull from the second try, but add the residuals back
    v <- voom(count, design, plot=F)
    v <- v$E
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    fourth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # fifth try: include asthma, but not the interaction, in the design
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition+income, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    fifth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # sixth try: same design, but now using removeBatchEffect function
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','condition_num','income_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition_num:income_num, data=sub_mdata))
    sixth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # seventh try: include everything including the interaction in the design
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition*income, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    seventh_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # eighth try: same design, but now using removeBatchEffect function
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','condition_num','income_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition_num:income_num, data=sub_mdata))
    eighth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # ninth try: a simpler design model
    design <- model.matrix(~0+batch, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    ninth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # tenth try: doing whatever it is when adding back batch effects from the design above?
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    tenth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # eleventh try: doing whatever it is when adding back effects from the design above, but now with everything but income
    design <- model.matrix(~0+age+gender+n+avg_mt+condition, data=sub_mdata)
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    eleventh_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # twelfth try: doing whatever it is when adding back effects from the design above, but now with everything
    design <- model.matrix(~0+age+gender+n+avg_mt+condition+income, data=sub_mdata)
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    twelfth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # thirteenth: the simplest model ever
    v <- voom(count, plot=F)
    v <- v$E
    thirteenth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # compute paired deltas
    first_ssgsea_scores <- first_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='1')
    
    second_ssgsea_scores <- second_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='2')
    
    third_ssgsea_scores <- third_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='3')
    
    fourth_ssgsea_scores <- fourth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='4')
    
    fifth_ssgsea_scores <- fifth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='5')
    
    sixth_ssgsea_scores <- sixth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='6')
    
    seventh_ssgsea_scores <- seventh_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='7')
    
    eighth_ssgsea_scores <- eighth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='8')
    
    ninth_ssgsea_scores <- ninth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='9')
    
    tenth_ssgsea_scores <- tenth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='10')
    
    eleventh_ssgsea_scores <- eleventh_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='11')
    
    twelfth_ssgsea_scores <- twelfth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='12')
    
    thirteenth_ssgsea_scores <- thirteenth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='13')
    
    # concatenate
    if (exists('ssgsea_scores')){
    tmp <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores, 
                           fifth_ssgsea_scores, sixth_ssgsea_scores, seventh_ssgsea_scores, eighth_ssgsea_scores,
                           ninth_ssgsea_scores, tenth_ssgsea_scores, eleventh_ssgsea_scores, twelfth_ssgsea_scores,
                 thirteenth_ssgsea_scores) %>% mutate(condition=cond)
    ssgsea_scores <- rbind(ssgsea_scores, tmp)
    } else {
      ssgsea_scores <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores, 
                             fifth_ssgsea_scores, sixth_ssgsea_scores, seventh_ssgsea_scores, eighth_ssgsea_scores,
                             ninth_ssgsea_scores, tenth_ssgsea_scores, eleventh_ssgsea_scores, twelfth_ssgsea_scores,
                             thirteenth_ssgsea_scores) %>% mutate(condition=cond)
    }
  }
}

# merge w metadata to retrieve income status
ssgsea_scores.wincome <- ssgsea_scores %>% inner_join(sample_m, by=c('ID')) %>% 
  select(ID, celltype, condition, deltaIFNa, deltaIFNy, income, method) %>% 
  pivot_longer(c(deltaIFNa, deltaIFNy), names_to='IFN', values_to='score') %>% 
  inner_join(mdata, by=c('ID'='IDs', 'celltype', 'condition', 'income'))
ssgsea_scores.wincome$method <- factor(ssgsea_scores.wincome$method, levels=c('1', '2', '3', '4', '5', '6', '7','8', '9',
                                                                              '10', '11', '12', '13'))
ssgsea_scores.wincome$income <- na_if(ssgsea_scores.wincome$income, '')
ssgsea_scores.wincome$income <- ifelse(ssgsea_scores.wincome$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Lower', 'Higher')
ssgsea_scores.wincome$income <- factor(ssgsea_scores.wincome$income, levels=c('Lower','Higher'))


ssgsea_scores.wincome %>% filter(celltype=='Mono') %>% ggplot(., aes(x=condition, y=score, fill=income)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = income), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ggtitle('IVA - Mono')

ssgsea_scores.wincome %>% filter(celltype=='T-CD4') %>% ggplot(., aes(x=condition, y=score, fill=income)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = income), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + ggtitle('IVA - T-CD4')

## compare ssGSEA and manual scores
# read asthma manual scores
ssgsea_scores.wasthma <- ssgsea_scores.wasthma %>% select(ID, celltype, condition, method, IFN, score)
manual_asthma_scores <- fread('IFNscores_asthma_manual.txt') %>% rename(ID=donor) %>% mutate(method='manual') %>% 
  select(ID, celltype, condition, method, IFN, score)
manual_asthma_scores$IFN <- gsub('IFNa', 'deltaIFNa', manual_asthma_scores$IFN)
manual_asthma_scores$IFN <- gsub('IFNy', 'deltaIFNy', manual_asthma_scores$IFN)
manual_asthma_scores <- manual_asthma_scores %>% filter(celltype=='T-CD8', condition=='RV')
joint_asthma <- rbind(ssgsea_scores.wasthma, manual_asthma_scores) %>% pivot_wider(names_from=method, values_from=score)
fwrite(joint_asthma, 'joint_asthma_scores.txt', sep=' ')

for (i in (5:ncol(joint_asthma))) {
  for (j in (5:ncol(joint_asthma))) {
    
    corr <- joint_asthma %>% group_by(IFN) %>% do({
      test <- cor.test(as.matrix(joint_asthma[,i]), as.matrix(joint_asthma[,j]), method='spearman') 
      test <- broom::tidy(test) %>% mutate(method1=colnames(joint_asthma)[i], method2=colnames(joint_asthma)[j])
    })
    
    if (exists('final.cor')){
      final.cor <- rbind(final.cor, corr)
    } else (final.cor <- corr)
  }
}
final.cor$method1 <- factor(final.cor$method1, levels=c('1', '2', '3', '4', '5', '6', '7','8', '9',
                                                        '10', '11', '12', '13', 'manual'))
final.cor$method2 <- factor(final.cor$method2, levels=c('1', '2', '3', '4', '5', '6', '7','8', '9',
                                                        '10', '11', '12', '13', 'manual'))
ggplot(final.cor, aes(x=method1, y=method2, fill=estimate)) + geom_tile() + facet_wrap(~IFN) +
  theme_bw()
ggsave('IFNscores_asthma_corr_bwtnmethods_heatmap.pdf', height=5, width=10)

## income
ssgsea_scores.wincome <- ssgsea_scores.wincome %>% select(ID, celltype, condition, method, IFN, score)
manual_income_scores <- fread('IFNscores_income_manual.txt') %>% rename(ID=donor) %>% mutate(method='manual') %>% 
  select(ID, celltype, condition, method, IFN, score) %>% filter(celltype %in% c('T-CD4', 'Mono'), condition=='IVA')
manual_income_scores$IFN <- gsub('IFNa', 'deltaIFNa', manual_income_scores$IFN)
manual_income_scores$IFN <- gsub('IFNy', 'deltaIFNy', manual_income_scores$IFN)
joint_income <- rbind(ssgsea_scores.wincome, manual_income_scores) %>% pivot_wider(names_from=method, values_from=score)
fwrite(joint_income, 'joint_income_scores.txt', sep=' ')

for (i in (5:ncol(joint_income))) {
  for (j in (5:ncol(joint_income))) {
    
    corr <- joint_income %>% group_by(IFN, celltype) %>% do({
      test <- cor.test(as.matrix(joint_income[,i]), as.matrix(joint_income[,j]), method='spearman') 
      test <- broom::tidy(test) %>% mutate(method1=colnames(joint_income)[i], method2=colnames(joint_income)[j])
    })
    
    if (exists('final.cor.in')){
      final.cor.in <- rbind(final.cor.in, corr)
    } else (final.cor.in <- corr)
  }
}
final.cor.in$method1 <- factor(final.cor.in$method1, levels=c('1', '2', '3', '4', '5', '6', '7','8', '9',
                                                             '10', '11', '12', '13', 'manual'))
final.cor.in$method2 <- factor(final.cor.in$method2, levels=c('1', '2', '3', '4', '5', '6', '7','8', '9',
                                                              '10', '11', '12', '13', 'manual'))
ggplot(final.cor.in, aes(x=method1, y=method2, fill=estimate)) + geom_tile() + facet_grid(cols=vars(IFN), rows=vars(celltype)) +
  theme_bw()
ggsave('IFNscores_income_corr_bwtnmethods_heatmap.pdf', height=10, width=10)
