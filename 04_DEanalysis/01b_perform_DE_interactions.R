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

# define minimum average logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','T-CD4','T-CD8','Mono','NK',
                                'B','T-CD4','T-CD8','Mono','NK'),
                         threshold=c(4.9,1.9,1,3.4,5.6,
                                 3.5,3.6,3.1,3.4,5.6),
                         condition=c(rep('IVA',5),rep('RV',5)))

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
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- factor(mdata$income, levels=c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', 
                                  '$50,000-$69,999', '$70,000-$89,999')) %>% as.numeric()
    no_NA_income <- mdata %>% filter(!is.na(income)) %>% rownames(.)
    no_NA_albuterol <- mdata %>% filter(!is.na(albuterol)) %>% rownames(.)
    
    for (interaction_term in c('asthma', 'income')){
      print(interaction_term)
      
      if (interaction_term=='asthma'){
        # remove non protein coding genes from count matrix and genes with variance == 0
        count <- tmp@assays$RNA$counts
        count <- count[,colnames(count) %in% no_NA_albuterol]
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        count <- DGEList(counts=count)
        
        # remove lowly expressed genes based on logCPM threshold
        logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
          pull(threshold)
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMeans() %>% as.data.frame() %>% filter(.>=logcpm_threshold) %>%
          rownames_to_column() %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        asthma_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_albuterol)
        design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition*asthma, data=asthma_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef='condition'%&%conditions[i]%&%':asthmaYes', number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])
        
        # now do permutations where i shuffle asthma status in metadata
        for (j in (1:10)){
          
          # shuffle asthma status preserving infection condition
          permuted_mdata <- asthma_mdata
          for (ind in unique(permuted_mdata$IDs)){
            # flip coin to decide if asthma status will be reversed or not
            if (runif(1)<0.5){
              ix <- permuted_mdata$IDs == ind
              if (permuted_mdata$asthma[ix][1]=='Yes'){
                permuted_mdata$asthma[ix] <- 'No'
              } else {
                permuted_mdata$asthma[ix] <- 'Yes'
              }
            }
          }
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition*asthma, data=permuted_mdata)
          
          # voom
          voom <- voom(count, design, plot=F)
          
          # save voom-adjusted expression table
          exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
          fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_asthma_alb_voom_expression.txt', sep=' ')
          rm(exp)
          
          # fit linear model 
          fit <- eBayes(lmFit(voom, design))
          
          # save pvalues from permutation
          tmp_perm <- topTable(fit, coef=ncol(fit), number=Inf) %>% rownames_to_column('Gene') %>%
            select(Gene, P.Value)
          
          if (exists('compiled_perms')){
            compiled_perms <- inner_join(compiled_perms, tmp_perm, by='Gene')
          } else {compiled_perms <- tmp_perm}
        }
        
        # reorder compiled_perms df so gene order matches OG results
        compiled_perms <- compiled_perms[match(og_results$Gene, compiled_perms$Gene), ]
        compiled_perms <- compiled_perms %>% select(-Gene)
        
        # compute qvalues
        empP <- empPvals(stat=-log10(og_results$P.Value), stat0=-log10(as.matrix(compiled_perms[1:j])), pool=TRUE)
        og_results$qvals <- qvalue(empP)$qvalue
        
        # qqplot 
        pdf('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_asthma_alb_limma_results_qqplot.pdf', width=4, height=4)
        qqplot(x=-log10(compiled_perms[,1]), y=-log10(og_results$P.Value), main=conditions[i]%&%' '%&%ctype%&%' asthma', 
               xlab='-log10(permuted p-values)', ylab='-log10(true p-values)')
        abline(c(0,1), col='red')
        dev.off()
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_asthma_alb_limma_results_wqvals.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
        
      # now, income
      } else {
        # remove non protein coding genes from count matrix and genes with variance == 0
        count <- tmp@assays$RNA$counts
        count <- count[,colnames(count) %in% no_NA_income]
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        count <- DGEList(counts=count)
        
        # remove lowly expressed genes based on logCPM threshold
        logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
          pull(threshold)
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMeans() %>% as.data.frame() %>% filter(.>=logcpm_threshold) %>%
          rownames_to_column() %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        income_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_income)
        design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income, data=income_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_income_voom_expression.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef='condition'%&%conditions[i]%&%':income', number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])

        # now do permutations where i shuffle income status in metadata
        for (j in (1:10)){
          
          # copy metadata
          permuted_mdata <- income_mdata
          
          # shuffle income status preserving infection condition
          income_map <- permuted_mdata %>% select(IDs, income) %>% distinct()
          shuffled_income <- sample(income_map$income)
          income_map$income_perm <- shuffled_income
          
          # join back to permuted metadata
          permuted_mdata <- permuted_mdata %>% 
            left_join(income_map %>% select(IDs, income_perm), by=c('IDs'))
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income_perm, data=permuted_mdata)
          
          # voom
          voom <- voom(count, design, plot=F)
          
          # fit linear model 
          fit <- eBayes(lmFit(voom, design))
          
          # save pvalues from permutation
          tmp_perm <- topTable(fit, coef=ncol(fit), number=Inf) %>% rownames_to_column('Gene') %>%
            select(Gene, P.Value)
          
          if (exists('compiled_perms')){
            compiled_perms <- inner_join(compiled_perms, tmp_perm, by='Gene')
          } else {compiled_perms <- tmp_perm}
        }
        
        # reorder compiled_perms df so gene order matches OG results
        compiled_perms <- compiled_perms[match(og_results$Gene, compiled_perms$Gene), ]
        compiled_perms <- compiled_perms %>% select(-Gene)
        
        # compute qvalues
        empP <- empPvals(stat=-log10(og_results$P.Value), stat0=-log10(as.matrix(compiled_perms[1:j])), pool=TRUE)
        og_results$qvals <- qvalue(empP)$qvalue
        
        # qqplot 
        pdf('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_income_limma_results_qqplot.pdf', width=4, height=4)
        qqplot(x=-log10(compiled_perms[,1]), y=-log10(og_results$P.Value), main=conditions[i]%&%' '%&%ctype%&%' income', 
               xlab='-log10(permuted p-values)', ylab='-log10(true p-values)')
        abline(c(0,1), col='red')
        dev.off()
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_income_limma_results_wqvals.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
      }
    }
  }
}
