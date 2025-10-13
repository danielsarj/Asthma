library(tidyverse)
library(data.table)
library(Seurat)
library(limma)
library(msigdbr)
'%&%' <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')

### GET SAMPLE METADATA ###
# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# load pseudobulk seurat object
bulk_obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds') 
bulk_obj@meta.data$condition <- factor(bulk_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# merge metadata
mdata <- bulk_obj@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
bulk_obj@meta.data <- mdata

# get human hallmark gene sets
ifn_genes <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)
ifn_genes <- ifn_genes[grep('INTERFERON', names(ifn_genes), value=T)]
names(ifn_genes) <- c('IFNa','IFNy')

### LOAD RESULTS ###
# load DE results
DE_results <- fread('NI_IVAxRV_integrated_limma_results.txt') %>% group_by(condition, interaction, celltype) %>% 
  arrange(desc(t), .by_group=T) %>% mutate(rank=row_number()) %>% ungroup()

# load GSEA results
GSEA_results <- fread('NI_IVAxRV_integrated_descGSEAresults.txt')

# get expression levels
for (int in c('none','asthma','income')){
  
  if (int == 'none'){
    for (cond in c('IVA', 'RV')){
      if (cond == 'IVA'){
        ctype <- 'T-CD8'
      } else { ctype <- 'Mono'}
      
      # subset seurat object
      tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
      sub_mdata <- tmp@meta.data
      sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    
      # get samples/genes from voom expression dataframe (easier to filter stuff)
      v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_voom_expression.txt') %>% 
        column_to_rownames('Gene')
      v_samples <- colnames(v)
      v_genes <- rownames(v)

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
      
      # average expression per infection status
      exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
        separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
        select(Gene, IDs, condition, celltype, value) %>% group_by(Gene, condition, celltype) %>% 
        summarise(mean_value=mean(value)) %>% ungroup() %>% group_by(celltype, condition) %>% 
        arrange(desc(mean_value), .by_group=T) %>% mutate(rank=row_number()) %>% ungroup()

      # see how rank changes
      inf_y <-  exp_long %>% filter(condition==cond) %>% select(Gene, rank) %>%
        mutate(type=cond)
      inf_n <-  exp_long %>% filter(condition=='NI') %>% select(Gene, rank) %>%
        mutate(type='NI')
      rank_genes <- rbind(inf_y, inf_n)
      
      (rank_genes %>% filter(Gene %in% ifn_genes[[1]]) %>% ggplot(., aes(x=rank, fill=type, color=type)) + 
          geom_density(aes(y=after_stat(count)), alpha=.4) + geom_rug() + theme_bw() + ggtitle('IFNa - ' %&% cond %&% ' - ' %&% ctype)) + 
        (rank_genes %>% filter(Gene %in% ifn_genes[[2]]) %>% ggplot(., aes(x=rank, fill=type, color=type)) + 
           geom_density(aes(y=after_stat(count)), alpha=.4) + geom_rug() + theme_bw() + ggtitle('IFNy - ' %&% cond %&% ' - ' %&% ctype))
      ggsave('IFN_gene_ranks_'%&%cond%&%'_infection_'%&%ctype%&%'_densityplots.pdf', height=4, width=9)
    }
  } else if (int == 'asthma'){
    cond <- 'RV'
    ctype <- 'T-CD8'
        
    # subset seurat object
     tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
     sub_mdata <- tmp@meta.data
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$albuterol <- factor(sub_mdata$albuterol, levels=c('No', 'Yes'))
      
     # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_asthma_alb_voom_expression.txt') %>% 
      column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)

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
      
    # average expression per asthma status
    exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
      separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
      inner_join(mdata, by=c('IDs', 'condition', 'celltype')) %>% 
      select(Gene, IDs, condition, celltype, value, asthma) %>% group_by(Gene, condition, celltype, asthma) %>% 
      summarise(mean_value=mean(value)) %>%
      ungroup() %>% group_by(Gene, celltype, asthma) %>% 
      summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>%
      ungroup() %>% group_by(celltype, asthma) %>% 
      arrange(desc(delta_value), .by_group=T) %>% mutate(rank=row_number()) %>% ungroup()
      
    # see how rank changes
    asthma_y <-  exp_long %>% filter(asthma=='Yes') %>% select(Gene, rank) %>%
      mutate(type='asthmaYes')
    asthma_n <-  exp_long %>% filter(asthma=='No') %>% select(Gene, rank) %>%
      mutate(type='asthmaNo')
    rank_genes <- rbind(asthma_y, asthma_n)
        
    (rank_genes %>% filter(Gene %in% ifn_genes[[1]]) %>% ggplot(., aes(x=rank, fill=type, color=type)) + 
      geom_density(aes(y=after_stat(count)), alpha=.4) + geom_rug() + theme_bw() + ggtitle('IFNa - ' %&% cond %&% ' - ' %&% ctype)) + 
      (rank_genes %>% filter(Gene %in% ifn_genes[[2]]) %>% ggplot(., aes(x=rank, fill=type, color=type)) + 
      geom_density(aes(y=after_stat(count)), alpha=.4) + geom_rug() + theme_bw() + ggtitle('IFNy - ' %&% cond %&% ' - ' %&% ctype))
    ggsave('IFN_gene_ranks_'%&%cond%&%'_asthma_'%&%ctype%&%'_densityplots.pdf', height=4, width=9)
  } else {
    for (ctype in c('Mono', 'T-CD4')){
      cond <- 'IVA'
      
      # subset seurat object
      tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
      sub_mdata <- tmp@meta.data
      sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
      sub_mdata$income <- ifelse(sub_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                                 'Lower', 'Higher')
      sub_mdata$income <- factor(sub_mdata$income, levels=c('Lower','Higher'))

      # get samples/genes from voom expression dataframe (easier to filter stuff)
      v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_income_voom_expression.txt') %>% 
        column_to_rownames('Gene')
      v_samples <- colnames(v)
      v_genes <- rownames(v)
      
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
      
      # average expression per asthma status
      exp_long <- corrected_expression %>% as.data.frame() %>% rownames_to_column('Gene') %>% reshape2::melt() %>%
        separate(variable, into=c('IDs', 'condition', 'celltype'), sep='_') %>% 
        inner_join(sub_mdata, by=c('IDs', 'condition', 'celltype')) %>% 
        select(Gene, IDs, condition, celltype, value, income) %>% group_by(Gene, condition, celltype, income) %>% 
        summarise(mean_value=mean(value)) %>%
        ungroup() %>% group_by(Gene, celltype, income) %>% 
        summarise(delta_value=.data$mean_value[condition==cond]-.data$mean_value[condition=='NI']) %>%
        ungroup() %>% group_by(celltype, income) %>% 
        arrange(desc(delta_value), .by_group=T) %>% mutate(rank=row_number()) %>% ungroup()
      
      # see how rank changes
      income_l <-  exp_long %>% filter(income=='Lower') %>% select(Gene, rank) %>%
        mutate(type='incomeLower')
      income_h <-  exp_long %>% filter(income=='Higher') %>% select(Gene, rank) %>%
        mutate(type='incomeHigher')
      rank_genes <- rbind(income_l, income_h)
      
      (rank_genes %>% filter(Gene %in% ifn_genes[[1]]) %>% ggplot(., aes(x=rank, fill=type, color=type)) + 
          geom_density(aes(y=after_stat(count)), alpha=.4) + geom_rug() + theme_bw() + ggtitle('IFNa - ' %&% cond %&% ' - ' %&% ctype)) + 
        (rank_genes %>% filter(Gene %in% ifn_genes[[2]]) %>% ggplot(., aes(x=rank, fill=type, color=type)) + 
           geom_density(aes(y=after_stat(count)), alpha=.4) + geom_rug() + theme_bw() + ggtitle('IFNy - ' %&% cond %&% ' - ' %&% ctype))
      ggsave('IFN_gene_ranks_'%&%cond%&%'_income_'%&%ctype%&%'_densityplots.pdf', height=4, width=9)
    }
  }
}
      
      