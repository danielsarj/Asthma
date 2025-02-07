library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')

# load sample metadata
sample_m <- fread('sample_metadata.txt')

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt')

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
  annotations$hgnc_symbol!='' &
  !grepl('^MT-', annotations$hgnc_symbol)]

for (i in 1:length(conditions)){
  print(c(conditions[i]))
  objs <- readRDS('NI_'%&%conditions[i]%&%'_pseudobulks.rds')
  
  # merge metadata
  mdata <- objs@meta.data
  mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
  objs@meta.data <- mdata
  
  # plots about the metadata
  summ <- mdata %>% group_by(condition, predicted.celltype.l1, asthma) %>% summarise(n=n())
  summ %>% filter(predicted.celltype.l1 %in% c('other', 'other-T')==FALSE) %>% 
    ggplot(.) + geom_col(aes(x=predicted.celltype.l1, y=n, fill=asthma), position='dodge') + 
    scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + facet_wrap(~condition)
  ggsave('NI_'%&%conditions[i]%&%'_SampleSizeByAsthmaStatus.pdf', height=3, width=7)
  
  summ <- mdata %>% group_by(condition, predicted.celltype.l1, income) %>% summarise(n=n())
  summ %>% filter(predicted.celltype.l1 %in% c('other', 'other-T')==FALSE) %>% 
    ggplot(.) + geom_col(aes(x=predicted.celltype.l1, y=n, fill=income), position='dodge') + 
    scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + facet_wrap(~condition)
  ggsave('NI_'%&%conditions[i]%&%'_SampleSizeByIncomeStatus.pdf', height=3, width=7)
  
  summ <- mdata %>% group_by(condition, predicted.celltype.l1, asthma, income) %>% summarise(n=n())
  summ %>% filter(predicted.celltype.l1 %in% c('other', 'other-T')==FALSE) %>% 
    ggplot(.) + geom_col(aes(x=predicted.celltype.l1, y=n, fill=condition), position='dodge') + 
    scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + 
    facet_grid(cols=vars(income), rows=vars(asthma))
  ggsave('NI_'%&%conditions[i]%&%'_SampleSizeByIncomeStatusANDAsthmaStatus.pdf', height=6, width=15)
  
  # celltype specific DE
  for (ctype in c('B','CD4-T','CD8-T','DC','Mono','NK')){
    print(ctype)
    # subset pseudobulk object
    tmp <- subset(objs, predicted.celltype.l1==ctype)
    
    # extract metadata and count matrices
    mdata <- tmp@meta.data
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$IDs <- as.factor(mdata$IDs)
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- factor(mdata$income, levels=c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', 
                                                  '$50,000-$69,999', '$70,000-$89,999'))
    
    for (interaction_term in c('asthma', 'income')){
      if (interaction_term=='asthma'){
        # filter count matrix (only keep protein coding genes)
        count <- tmp@assays$RNA$counts
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        
        # transform count into dge object
        count <- DGEList(counts=count)
        count <- calcNormFactors(count)
        
        # define design matrix
        design <- model.matrix(~0+age+condition+condition:asthma, data=mdata)
        
        # voom
        voom <- voom(count, design, plot=T)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        ni_results <- topTable(fit, coef='conditionNI:asthmaYes', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition='NI')
        results <- topTable(fit, coef='condition'%&%conditions[i]%&%':asthmaYes', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition=conditions[i])
        results <- rbind(ni_results, results)
        
        fwrite(results, '../DEanalysis/NI_'%&%conditions[i]%&%'_asthma_limma_'%&%ctype%&%'_results.txt',
               sep=' ', col.names=T, row.names=T)
      } else {
        # filter count matrix (only keep protein coding genes)
        count <- tmp@assays$RNA$counts
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        
        # transform count into dge object
        count <- DGEList(counts=count)
        count <- calcNormFactors(count)
        
        # define design matrix
        design <- model.matrix(~0+age+condition+condition:income, data=mdata)
        
        # voom
        voom <- voom(count, design, plot=T)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        tmp1_results <- topTable(fit, coef='conditionNI:income$10,000-$29,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition='NI', income='10,000-$29,999')
        tmp2_results <- topTable(fit, coef='conditionNI:income$30,000-$49,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition='NI', income='30,000-$49,999')
        tmp3_results <- topTable(fit, coef='conditionNI:income$50,000-$69,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition='NI', income='50,000-$69,999')
        tmp4_results <- topTable(fit, coef='conditionNI:income$70,000-$89,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition='NI', income='70,000-$89,999')
        tmp5_results <- topTable(fit, coef='condition'%&%conditions[i]%&%':income$10,000-$29,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition=conditions[i], income='10,000-$29,999')
        tmp6_results <- topTable(fit, coef='condition'%&%conditions[i]%&%':income$30,000-$49,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition=conditions[i], income='30,000-$49,999')
        tmp7_results <- topTable(fit, coef='condition'%&%conditions[i]%&%':income$50,000-$69,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition=conditions[i], income='50,000-$69,999')
        tmp8_results <- topTable(fit, coef='condition'%&%conditions[i]%&%':income$70,000-$89,999', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition=conditions[i], income='70,000-$89,999')
        
        results <- rbind(tmp1_results, tmp2_results, tmp3_results, tmp4_results,
                         tmp5_results, tmp6_results, tmp7_results, tmp8_results)
        
        fwrite(results, '../DEanalysis/NI_'%&%conditions[i]%&%'_income_limma_'%&%ctype%&%'_results.txt',
               sep=' ', col.names=T, row.names=T)
      }
    }
  }
}