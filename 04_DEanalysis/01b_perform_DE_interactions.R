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
sample_m <- fread('../sample_metadata.txt')

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt')

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# load seurat object
objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')

# merge metadata
mdata <- objs@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
objs@meta.data <- mdata

# plots about the metadata
summ <- mdata %>% group_by(condition, celltype, asthma) %>% summarise(n=n())
summ %>% ggplot(.) + geom_col(aes(x=celltype, y=n, fill=asthma), position='dodge') + 
  scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + facet_wrap(~condition)
ggsave('SampleSizeByAsthmaStatus.pdf', height=4, width=8)

summ <- mdata %>% group_by(condition, celltype, income) %>% summarise(n=n())
summ %>% drop_na() %>% ggplot(.) + geom_col(aes(x=celltype, y=n, fill=income), position='dodge') + 
  scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + facet_wrap(~condition)
ggsave('SampleSizeByIncomeStatus.pdf', height=4, width=8)

summ <- mdata %>% group_by(condition, celltype, asthma, income) %>% summarise(n=n())
summ %>% drop_na() %>% ggplot(.) + geom_col(aes(x=celltype, y=n, fill=condition), position='dodge') + 
  scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + 
  facet_grid(cols=vars(income), rows=vars(asthma))
ggsave('SampleSizeByIncomeStatusANDAsthmaStatus.pdf', height=6, width=12)

summ <- mdata %>% group_by(condition, celltype, asthma, albuterol) %>% summarise(n=n()) %>%
  ungroup() %>% mutate(across('albuterol', ~na_if(., '')))
summ %>% drop_na() %>% ggplot(.) + geom_col(aes(x=celltype, y=n, fill=condition), position='dodge') + 
  scale_y_continuous(breaks=seq(0, max(summ$n), by=1)) + theme_bw() + 
  facet_grid(cols=vars(albuterol), rows=vars(asthma))
ggsave('SampleSizeByAlbuterolStatusANDAsthmaStatus.pdf', height=6, width=7)

# condition specific DE
for (i in 1:length(conditions)){
  print(c(conditions[i]))
  
  # celltype specific DE
  for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
    print(ctype)
    
    # extract metadata for subsetting
    meta_df <- objs@meta.data
    filtered_meta <- meta_df %>% filter(celltype==ctype, condition %in% c(conditions[i], 'NI'))
    
    # subset bulk object
    matching_cells <- rownames(filtered_meta)
    tmp <- subset(objs, cells=matching_cells)
    rm(meta_df, filtered_meta, matching_cells)
    
    # extract metadata
    mdata <- tmp@meta.data
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    mdata$IDs <- as.factor(mdata$IDs)
    mdata$albuterol <- na_if(mdata$albuterol, '')
    mdata$albuterol <- factor(mdata$albuterol, levels=c('No', 'Yes'))
    mdata$asthma <- factor(mdata$asthma, levels=c('No', 'Yes'))
    mdata$income <- factor(mdata$income, levels=c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', 
                                                  '$50,000-$69,999', '$70,000-$89,999')) %>% as.numeric()
    no_NA_income <- mdata %>% filter(!is.na(income)) %>% rownames(.)
    no_NA_albuterol <- mdata %>% filter(!is.na(albuterol)) %>% rownames(.)
    
    for (interaction_term in c('asthma', 'income')){
      print(interaction_term)
      
      if (interaction_term=='asthma'){
        
        # if interaction is asthma, run DE depending on albuterol intake as well
        for (alb in c('no','yes')){
          
          # dont account for albuterol intake
          if (alb=='no'){
            # filter count matrix (only keep protein coding genes)
            count <- tmp@assays$RNA$counts
            count <- count[rownames(count) %in% annotations,]
            zero_var_genes <- apply(count, 1, var) == 0
            count <- count[!zero_var_genes, ]
            
            # transform count into dge object
            count <- DGEList(counts=count)
            count <- calcNormFactors(count)
            
            # define design matrix
            design <- model.matrix(~batch+age+gender+n+avg_mt+condition+condition:asthma, data=mdata)
            
            # voom
            voom <- voom(count, design, plot=F)
            
            # fit linear model 
            fit <- eBayes(lmFit(voom, design))
            
            # get results
            results <- topTable(fit, coef='condition'%&%conditions[i]%&%':asthmaYes', number=Inf, adjust='BH') %>% 
              rownames_to_column('gene') %>% mutate(condition=conditions[i])
            fwrite(results, 'NI_'%&%conditions[i]%&%'_asthma_limma_'%&%ctype%&%'_results.txt',
                   sep=' ', col.names=T, na='NA')
            
            # adjust for albuterol intake
          } else {
            # filter count matrix (only keep protein coding genes)
            count <- tmp@assays$RNA$counts
            count <- count[,colnames(count) %in% no_NA_albuterol]
            count <- count[rownames(count) %in% annotations,]
            zero_var_genes <- apply(count, 1, var) == 0
            count <- count[!zero_var_genes, ]
            
            # transform count into dge object
            count <- DGEList(counts=count)
            count <- calcNormFactors(count)
            
            # define design matrix
            design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition+condition:asthma, data=mdata)
            
            # voom
            voom <- voom(count, design, plot=F)
            
            # fit linear model 
            fit <- eBayes(lmFit(voom, design))
            
            # get results
            results <- topTable(fit, coef='condition'%&%conditions[i]%&%':asthmaYes', number=Inf, adjust='BH') %>% 
              rownames_to_column('gene') %>% mutate(condition=conditions[i])
            fwrite(results, 'NI_'%&%conditions[i]%&%'_asthma_alb_limma_'%&%ctype%&%'_results.txt',
                   sep=' ', col.names=T, na='NA')
          }
        }
        
        
      } else {
        # filter count matrix (only keep protein coding genes)
        count <- tmp@assays$RNA$counts
        count <- count[,colnames(count) %in% no_NA_income]
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        
        # transform count into dge object
        count <- DGEList(counts=count)
        count <- calcNormFactors(count)
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+condition+condition:income, data=mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        results <- topTable(fit, coef='condition'%&%conditions[i]%&%':income', number=Inf, adjust='BH') %>% 
          rownames_to_column('gene') %>% mutate(condition=conditions[i])
        fwrite(results, 'NI_'%&%conditions[i]%&%'_income_limma_'%&%ctype%&%'_results.txt',
               sep=' ', col.names=T, na='NA')
      }
    }
  }
}
