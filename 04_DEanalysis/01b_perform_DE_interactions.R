library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(qvalue)
library(matrixStats)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
celltype_vector <- c('B','CD4-T','CD8-T','Mono','NK')
interactions_vector <- c('asthma', 'income', 'ACT', 'ACE', 'resilience', 'social_support', 'total_racism',
  'year_racism', 'life_racism', 'racism_stress', 'racism_child_24hr', 'kid_discrimination', 'infection_at_collection')

# load gene annotation from ensembl
annotations <- fread('ensembl_genes.txt')

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# read seurat object
objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks_new.rds')

# remove batch 4
objs <- subset(objs, subset= batch!='B4')

# define minimum average logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','CD4-T','CD8-T','Mono','NK',
                                            'B','CD4-T','CD8-T','Mono','NK'),
                                 threshold=c(2.7,-0.5,1.0,3.9,2.8,
                                             4.0,-0.7,1.9,4.0,2.8),
                                 condition=c(rep('IVA',5),rep('RV',5)))

# condition specific DE
for (i in 1:length(conditions)){
  print(c(conditions[i]))
  
  # celltype specific DE
  for (ctype in celltype_vector){
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
    mdata$batch <- factor(mdata$batch, levels=c('B1','B2','B3'))
    mdata$income <- factor(mdata$income, levels=c('Low','High'))
    mdata$albuterol <- factor(mdata$albuterol, levels=c('No','Yes'))
    mdata$Recorded_Diagnosis <- factor(mdata$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
    mdata$infection_status <- factor(mdata$infection_status, levels=c('Negative', 'Positive'))
    
    no_NA_albuterol <- mdata %>% filter(!is.na(albuterol)) %>% rownames(.)
    no_NA_income <- mdata %>% filter(!is.na(income)) %>% rownames(.)
    no_NA_ACT <- mdata %>% filter(!is.na(ACT_score)) %>% rownames(.)
    no_NA_Parent_Resilience <- mdata %>% filter(!is.na(Parent_Resilience_Score)) %>% rownames(.)
    no_NA_Parent_support <- mdata %>% filter(!is.na(Parents_Score_Avg)) %>% rownames(.)
    no_NA_Racism_child_24hr <- mdata %>% filter(!is.na(Racism_child_24hr)) %>% rownames(.)
    no_NA_Discrimination_child <- mdata %>% filter(!is.na(Experience_Discrimination_child)) %>% rownames(.)
    no_NA_infection_status <- mdata %>% filter(!is.na(infection_status)) %>% rownames(.)
    
    for (interaction_term in interactions_vector){
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
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
        logCPM_pass$rowname <- rownames(count$counts)
        logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        asthma_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_albuterol)
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+albuterol+condition*Recorded_Diagnosis, data=asthma_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_asthma_alb_voom_expression_new.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])
        
        # now do permutations where i shuffle asthma status in metadata
        for (j in (1:10)){
          
          # shuffle asthma status preserving infection condition
          permuted_mdata <- asthma_mdata
          for (ind in unique(permuted_mdata$IDs)){
            # flip coin to decide if asthma status will be reversed or not
            if (runif(1)<0.5){
              ix <- permuted_mdata$IDs == ind
              if (permuted_mdata$Recorded_Diagnosis[ix][1]=='Recorded_Asthma_Diagnosis'){
                permuted_mdata$Recorded_Diagnosis[ix] <- 'No_Diagnosis'
              } else {
                permuted_mdata$Recorded_Diagnosis[ix] <- 'Recorded_Asthma_Diagnosis'
              }
            }
          }
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+prop+albuterol+condition*Recorded_Diagnosis, data=permuted_mdata)
          
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

        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_asthma_alb_limma_results_wqvals_new.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
        
      # now, income
      } else if (interaction_term=='income'){
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
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
        logCPM_pass$rowname <- rownames(count$counts)
        logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        income_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_income)
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*income, data=income_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_income_voom_expression_new.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])

        # now do permutations where i shuffle income status in metadata
        for (j in (1:10)){
          
          # shuffle income status preserving infection condition
          permuted_mdata <- income_mdata
          for (ind in unique(permuted_mdata$IDs)){
            # flip coin to decide if income status will be reversed or not
            if (runif(1)<0.5){
              ix <- permuted_mdata$IDs == ind
              if (permuted_mdata$income[ix][1]=='Low'){
                permuted_mdata$income[ix] <- 'High'
              } else {
                permuted_mdata$income[ix] <- 'Low'
              }
            }
          }
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*income, data=permuted_mdata)
          
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
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_income_limma_results_wqvals_new.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
        
        # now do ACT
      } else if (interaction_term=='ACT'){
        # remove non protein coding genes from count matrix and genes with variance == 0
        count <- tmp@assays$RNA$counts
        count <- count[,colnames(count) %in% no_NA_ACT]
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        count <- DGEList(counts=count)
        
        # remove lowly expressed genes based on logCPM threshold
        logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
          pull(threshold)
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
        logCPM_pass$rowname <- rownames(count$counts)
        logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        act_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_ACT)
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*ACT_score, data=act_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_ACT_voom_expression_new.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))

        # get results
        og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])
        
        # now do permutations where i shuffle ACT scores in metadata
        for (j in (1:10)){
          
          # shuffle ACT scores preserving infection condition
          permuted_mdata <- act_mdata
          shuffled_scores <- permuted_mdata %>% select(IDs, ACT_score) %>% distinct(IDs, .keep_all=TRUE) %>%
            mutate(ACT_score = sample(ACT_score)) 
          permuted_mdata <- permuted_mdata %>% select(-ACT_score) %>%   # remove old score
            left_join(shuffled_scores, by='IDs')
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*ACT_score, data=permuted_mdata)
          
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
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_ACT_limma_results_wqvals_new.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
      
        # now do ACE
      } else if (interaction_term=='ACE'){
        # remove non protein coding genes from count matrix and genes with variance == 0
        count <- tmp@assays$RNA$counts
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        count <- DGEList(counts=count)
        
        # remove lowly expressed genes based on logCPM threshold
        logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
          pull(threshold)
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
        logCPM_pass$rowname <- rownames(count$counts)
        logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        ace_mdata <- mdata
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*ACE_result, data=ace_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_ACE_voom_expression_new.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])
        
        # now do permutations where i shuffle ACE scores in metadata
        for (j in (1:10)){
          
          # shuffle ACE scores preserving infection condition
          permuted_mdata <- ace_mdata
          shuffled_scores <- permuted_mdata %>% select(IDs, ACE_result) %>% distinct(IDs, .keep_all=TRUE) %>%
            mutate(ACE_result = sample(ACE_result)) 
          permuted_mdata <- permuted_mdata %>% select(-ACE_result) %>%   # remove old score
            left_join(shuffled_scores, by='IDs')
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*ACE_result, data=permuted_mdata)
          
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
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_ACE_limma_results_wqvals_new.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
      
        # now do parent resilience score
      } else if (interaction_term=='resilience'){
        # remove non protein coding genes from count matrix and genes with variance == 0
        count <- tmp@assays$RNA$counts
        count <- count[,colnames(count) %in% no_NA_Parent_Resilience]
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        count <- DGEList(counts=count)
        
        # remove lowly expressed genes based on logCPM threshold
        logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
          pull(threshold)
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
        logCPM_pass$rowname <- rownames(count$counts)
        logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        resilience_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_Parent_Resilience)
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Parent_Resilience_Score, data=resilience_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_resilience_voom_expression_new.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])
        
        # now do permutations where i shuffle resilience scores in metadata
        for (j in (1:10)){
          
          # shuffle resilience scores preserving infection condition
          permuted_mdata <- resilience_mdata
          shuffled_scores <- permuted_mdata %>% select(IDs, Parent_Resilience_Score) %>% distinct(IDs, .keep_all=TRUE) %>%
            mutate(Parent_Resilience_Score = sample(Parent_Resilience_Score)) 
          permuted_mdata <- permuted_mdata %>% select(-Parent_Resilience_Score) %>%   # remove old score
            left_join(shuffled_scores, by='IDs')
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Parent_Resilience_Score, data=permuted_mdata)
          
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
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_resilience_limma_results_wqvals_new.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
      
        # now do parent social support
      } else if (interaction_term=='social_support'){
        # remove non protein coding genes from count matrix and genes with variance == 0
        count <- tmp@assays$RNA$counts
        count <- count[,colnames(count) %in% no_NA_Parent_support]
        count <- count[rownames(count) %in% annotations,]
        zero_var_genes <- apply(count, 1, var) == 0
        count <- count[!zero_var_genes, ]
        count <- DGEList(counts=count)
        
        # remove lowly expressed genes based on logCPM threshold
        logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
          pull(threshold)
        logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
        logCPM_pass$rowname <- rownames(count$counts)
        logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
        count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
        count <- calcNormFactors(count)
        
        # define design matrix
        social_support_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_Parent_support)
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Parents_Score_Avg, data=social_support_mdata)
        
        # voom
        voom <- voom(count, design, plot=F)
        
        # save voom-adjusted expression table
        exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
        fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_social_support_voom_expression_new.txt', sep=' ')
        rm(exp)
        
        # fit linear model 
        fit <- eBayes(lmFit(voom, design))
        
        # get results
        og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
          rownames_to_column('Gene') %>% mutate(condition=conditions[i])
        
        # now do permutations where i shuffle social support scores in metadata
        for (j in (1:10)){
          
          # shuffle social support scores preserving infection condition
          permuted_mdata <- social_support_mdata
          shuffled_scores <- permuted_mdata %>% select(IDs, Parents_Score_Avg) %>% distinct(IDs, .keep_all=TRUE) %>%
            mutate(Parents_Score_Avg = sample(Parents_Score_Avg)) 
          permuted_mdata <- permuted_mdata %>% select(-Parents_Score_Avg) %>%   # remove old score
            left_join(shuffled_scores, by='IDs')
          
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Parents_Score_Avg, data=permuted_mdata)
          
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
        
        # save result
        fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_social_support_limma_results_wqvals_new.txt',
               sep=' ', col.names=T, na='NA')
        rm(compiled_perms)
        
      # now do total racism
    } else if (interaction_term=='total_racism'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      total_racism_mdata <- mdata 
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Total_Racist_Events, data=total_racism_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_total_racism_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle total racism scores in metadata
      for (j in (1:10)){
        
        # shuffle total racism scores preserving infection condition
        permuted_mdata <- total_racism_mdata
        shuffled_scores <- permuted_mdata %>% select(IDs, Total_Racist_Events) %>% distinct(IDs, .keep_all=TRUE) %>%
          mutate(Total_Racist_Events = sample(Total_Racist_Events)) 
        permuted_mdata <- permuted_mdata %>% select(-Total_Racist_Events) %>%   # remove old score
          left_join(shuffled_scores, by='IDs')
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Total_Racist_Events, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_total_racism_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
    
      # now do year racism
    } else if (interaction_term=='year_racism'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      year_racism_mdata <- mdata 
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Year_Racist_events, data=year_racism_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_year_racism_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle year racism scores in metadata
      for (j in (1:10)){
        
        # shuffle year racism scores preserving infection condition
        permuted_mdata <- year_racism_mdata
        shuffled_scores <- permuted_mdata %>% select(IDs, Year_Racist_events) %>% distinct(IDs, .keep_all=TRUE) %>%
          mutate(Year_Racist_events = sample(Year_Racist_events)) 
        permuted_mdata <- permuted_mdata %>% select(-Year_Racist_events) %>%   # remove old score
          left_join(shuffled_scores, by='IDs')
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Year_Racist_events, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_year_racism_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
    
      # now do life racism
    } else if (interaction_term=='life_racism'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      life_racism_mdata <- mdata 
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Life_Racist_events, data=life_racism_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_life_racism_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle life racism scores in metadata
      for (j in (1:10)){
        
        # shuffle life racism scores preserving infection condition
        permuted_mdata <- life_racism_mdata
        shuffled_scores <- permuted_mdata %>% select(IDs, Life_Racist_events) %>% distinct(IDs, .keep_all=TRUE) %>%
          mutate(Life_Racist_events = sample(Life_Racist_events)) 
        permuted_mdata <- permuted_mdata %>% select(-Life_Racist_events) %>%   # remove old score
          left_join(shuffled_scores, by='IDs')
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Life_Racist_events, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_life_racism_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
   
      # now do stress racism
    } else if (interaction_term=='racism_stress'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      stress_racism_mdata <- mdata 
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Racist_stress, data=stress_racism_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_stress_racism_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle stress racism scores in metadata
      for (j in (1:10)){
        
        # shuffle stress racism scores preserving infection condition
        permuted_mdata <- stress_racism_mdata
        shuffled_scores <- permuted_mdata %>% select(IDs, Racist_stress) %>% distinct(IDs, .keep_all=TRUE) %>%
          mutate(Racist_stress = sample(Racist_stress)) 
        permuted_mdata <- permuted_mdata %>% select(-Racist_stress) %>%   # remove old score
          left_join(shuffled_scores, by='IDs')
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Racist_stress, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_stress_racism_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
    
      # now do past 24h kid racism
    } else if (interaction_term=='racism_child_24hr'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[,colnames(count) %in% no_NA_Racism_child_24hr]
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      kid24_racism_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_Racism_child_24hr)
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Racism_child_24hr, data=kid24_racism_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_kid_24h_racism_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle kid 24h racism scores in metadata
      for (j in (1:10)){
        
        # shuffle kid 24h racism scores preserving infection condition
        permuted_mdata <- kid24_racism_mdata
        shuffled_scores <- permuted_mdata %>% select(IDs, Racism_child_24hr) %>% distinct(IDs, .keep_all=TRUE) %>%
          mutate(Racism_child_24hr = sample(Racism_child_24hr)) 
        permuted_mdata <- permuted_mdata %>% select(-Racism_child_24hr) %>%   # remove old score
          left_join(shuffled_scores, by='IDs')
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Racism_child_24hr, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_kid_24h_racism_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
    
      # now do kid discrimination
    } else if (interaction_term=='kid_discrimination'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[,colnames(count) %in% no_NA_Discrimination_child]
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      kid_discrimination_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_Discrimination_child)
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Experience_Discrimination_child, data=kid_discrimination_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_kid_discrimination_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle kid discrimination scores in metadata
      for (j in (1:10)){
        
        # shuffle kid discrimination scores preserving infection condition
        permuted_mdata <- kid_discrimination_mdata
        shuffled_scores <- permuted_mdata %>% select(IDs, Experience_Discrimination_child) %>% distinct(IDs, .keep_all=TRUE) %>%
          mutate(Experience_Discrimination_child = sample(Experience_Discrimination_child)) 
        permuted_mdata <- permuted_mdata %>% select(-Experience_Discrimination_child) %>%   # remove old score
          left_join(shuffled_scores, by='IDs')
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*Experience_Discrimination_child, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_kid_discrimination_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
    
      # now do infection status at collection
    } else if (interaction_term=='infection_at_collection'){
      # remove non protein coding genes from count matrix and genes with variance == 0
      count <- tmp@assays$RNA$counts
      count <- count[,colnames(count) %in% no_NA_infection_status]
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
      count <- DGEList(counts=count)
      
      # remove lowly expressed genes based on logCPM threshold
      logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
        pull(threshold)
      logCPM_pass <- cpm(count, log=TRUE) %>% rowMedians() %>% as.data.frame() %>% rownames_to_column() 
      logCPM_pass$rowname <- rownames(count$counts)
      logCPM_pass <- logCPM_pass %>% filter(.>=logcpm_threshold) %>% pull(rowname)
      count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
      count <- calcNormFactors(count)
      
      # define design matrix
      infection_collection_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_infection_status)
      design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*infection_status, data=infection_collection_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # save voom-adjusted expression table
      exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
      fwrite(exp, '../scRNAanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_collection_infection_voom_expression_new.txt', sep=' ')
      rm(exp)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
      
      # get results
      og_results <- topTable(fit, coef=ncol(fit), number=Inf, adjust='BH') %>% 
        rownames_to_column('Gene') %>% mutate(condition=conditions[i])
      
      # now do permutations where i shuffle collection infection in metadata
      for (j in (1:10)){
        
        # shuffle collection infection labels preserving infection condition
        permuted_mdata <- infection_collection_mdata
        for (ind in unique(permuted_mdata$IDs)){
          # flip coin to decide if asthma status will be reversed or not
          if (runif(1)<0.5){
            ix <- permuted_mdata$IDs == ind
            if (permuted_mdata$infection_status[ix][1]=='Positive'){
              permuted_mdata$infection_status[ix] <- 'Negative'
            } else {
              permuted_mdata$infection_status[ix] <- 'Positive'
            }
          }
        }
        
        # define design matrix
        design <- model.matrix(~batch+age+gender+n+avg_mt+prop+condition*infection_status, data=permuted_mdata)
        
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
      
      # save result
      fwrite(og_results, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_collection_infection_limma_results_wqvals_new.txt',
             sep=' ', col.names=T, na='NA')
      rm(compiled_perms)
    }
    }
  }
}
