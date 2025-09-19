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
conditions <- c('RV')
celltypes <- c('T-CD8')

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

for (cond in conditions){
  for (ctype in celltypes){
    
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
    
    # fourth try: include asthma, but not the interaction, in the design
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition+asthma, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    fourth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # fifth try: same design, but now using removeBatchEffect function
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num','condition_num','asthma_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition_num:asthma_num, data=sub_mdata))
    fifth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # sixth try: include everything including the interaction in the design
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition*asthma, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    sixth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # seventh try: same design, but now using removeBatchEffect function
    v <- removeBatchEffect(v,
                           covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender_num','albuterol_num','condition_num','asthma_num')]),
                           batch = sub_mdata$batch,
                           design = model.matrix(~condition_num:asthma_num, data=sub_mdata))
    seventh_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    # eighth try: a simpler design model
    design <- model.matrix(~0+batch, data=sub_mdata)
    v <- voom(count, design, plot=F)
    v <- v$E
    eighth_ssgsea_scores <- gsva(ssgseaParam(v, ifn_genes))
    
    
    # ninth try: doing whatever it is when adding back batch effects from the design above?
    v <- voom(count, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    ninth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    # tenth try: doing whatever it is when adding back effects from the design above, but now with everything but asthma
    design <- model.matrix(~0+age+gender+n+avg_mt+albuterol+condition, data=sub_mdata)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    residuals <- residuals.MArrayLM(fit, v)
    avg_batch_effect <- rowMeans(fit$coefficients)
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    tenth_ssgsea_scores <- gsva(ssgseaParam(corrected_expression, ifn_genes))
    
    
    
    
    
    
    # compute paired deltas
    first_ssgsea_scores <- first_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='first')
    
    second_ssgsea_scores <- second_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='second')
    
    third_ssgsea_scores <- third_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='third')
    
    fourth_ssgsea_scores <- fourth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='fourth')
    
    fifth_ssgsea_scores <- fifth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='fifth')
    
    sixth_ssgsea_scores <- sixth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='sixth')
    
    seventh_ssgsea_scores <- seventh_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='seventh')
    
    eighth_ssgsea_scores <- eighth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='eighth')
    
    ninth_ssgsea_scores <- ninth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='ninth')
    
    tenth_ssgsea_scores <- tenth_ssgsea_scores %>% t() %>% as.data.frame() %>% rownames_to_column('temp') %>% 
      separate(temp, c('ID', 'condition', 'celltype'), '_') %>% group_by(ID, celltype) %>%
      summarise(deltaIFNa = IFNa[which(condition==cond)] - IFNa[which(condition=='NI')],
                deltaIFNy = IFNy[which(condition==cond)] - IFNy[which(condition=='NI')]) %>% ungroup() %>%
      mutate(method='tenth')
    
    # concatenate
    ssgsea_scores <- rbind(first_ssgsea_scores, second_ssgsea_scores, third_ssgsea_scores, fourth_ssgsea_scores, 
                           fifth_ssgsea_scores, sixth_ssgsea_scores, seventh_ssgsea_scores, eighth_ssgsea_scores,
                           ninth_ssgsea_scores, tenth_ssgsea_scores) %>% mutate(condition=cond)
  }
}

# merge w metadata to retrieve asthma status
ssgsea_scores.wasthma <- ssgsea_scores %>% inner_join(sample_m, by=c('ID')) %>% 
  select(ID, celltype, condition, deltaIFNa, deltaIFNy, asthma, method) %>% 
  pivot_longer(c(deltaIFNa, deltaIFNy), names_to='IFN', values_to='score') %>% 
  inner_join(mdata, by=c('ID'='IDs', 'celltype', 'condition', 'asthma'))
ssgsea_scores.wasthma$method <- factor(ssgsea_scores.wasthma$method, levels=c('first', 'second', 'third',
                                                                              'fourth', 'fifth', 'sixth', 
                                                                              'seventh', 'eighth', 'ninth',
                                                                              'tenth'))

ggplot(ssgsea_scores.wasthma, aes(x=condition, y=score, fill=asthma)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(method), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))



results <- ssgsea_scores.wasthma %>%
  group_by(IFN, method) %>%
  summarise(
    t_pval = t.test(score ~ asthma)$p.value,
    mean_asthma = mean(score[asthma == "Yes"], na.rm = TRUE),
    mean_control = mean(score[asthma == "No"], na.rm = TRUE),
    .groups = "drop"
  )


lm_results <- ssgsea_scores.wasthma %>%
  group_by(celltype, IFN, method) %>%
  do({
    fit <- lm(score ~ asthma + age + gender + batch + n + avg_mt + albuterol, data = .)
    broom::tidy(fit)   # gives coef, p-value, etc.
  }) %>%
  ungroup()






