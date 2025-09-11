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
    mdata$income <- na_if(mdata$income, '')
    mdata$income <- factor(mdata$income, levels=c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', 
                                                  '$50,000-$69,999', '$70,000-$89,999'))
    no_NA_income <- mdata %>% filter(!is.na(income)) %>% rownames(.)

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
      
    # design 1: income as factors - dummy variables
    print('design 1')
    # define design matrix
    income_mdata <- mdata %>% filter(rownames(mdata) %in% no_NA_income)
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income, data=income_mdata)
        
    # voom
    voom <- voom(count, design, plot=T)
        
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # compute average residual standard deviation and median AIC across genes 
    d_1 <- mean(fit$sigma)
    
    AIC_all <- sapply(1:nrow(fit$coefficients), function(i) {
      sigma_i <- fit$sigma[i]        # residual standard deviation for gene i
      RSS <- sigma_i^2 * (ncol(voom$E) - ncol(fit$design))     # residual sum of squares
      ncol(voom$E) * log(RSS/ncol(voom$E)) + 2*ncol(fit$design)          # AIC formula
    })
    aic_1 <- median(AIC_all)
    
    # design 2: income as numeric values - one continuous scale
    print('design 2')
    income_mdata$income_num <- as.numeric(income_mdata$income)
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income_num, data=income_mdata)
    
    # voom
    voom <- voom(count, design, plot=T)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # compute average residual standard deviation and median AIC across genes 
    d_2 <- mean(fit$sigma)
    
    AIC_all <- sapply(1:nrow(fit$coefficients), function(i) {
      sigma_i <- fit$sigma[i]        # residual standard deviation for gene i
      RSS <- sigma_i^2 * (ncol(voom$E) - ncol(fit$design))     # residual sum of squares
      ncol(voom$E) * log(RSS/ncol(voom$E)) + 2*ncol(fit$design)          # AIC formula
    })
    aic_2 <- median(AIC_all)
    
    # design 3: income as binary variable (< $10,000 vs everything else)
    print('design 3')
    income_mdata$income_binary <- ifelse(income_mdata$income == '< $10,000', 'LOW', 'HIGH')
    income_mdata$income_binary <- factor(income_mdata$income_binary, levels=c('LOW','HIGH'))
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income_binary, data=income_mdata)
    
    # voom
    voom <- voom(count, design, plot=T)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # compute average residual standard deviation and median AIC across genes 
    d_3 <- mean(fit$sigma)
    
    AIC_all <- sapply(1:nrow(fit$coefficients), function(i) {
      sigma_i <- fit$sigma[i]        # residual standard deviation for gene i
      RSS <- sigma_i^2 * (ncol(voom$E) - ncol(fit$design))     # residual sum of squares
      ncol(voom$E) * log(RSS/ncol(voom$E)) + 2*ncol(fit$design)          # AIC formula
    })
    aic_3 <- median(AIC_all)
    
    # design 4: income as binary variable (< $10,000 and $10,000-$29,999 vs everything else)
    print('design 4')
    income_mdata$income_binary <- ifelse(income_mdata$income %in% c('< $10,000', '$10,000-$29,999'), 'LOW', 'HIGH')
    income_mdata$income_binary <- factor(income_mdata$income_binary, levels=c('LOW','HIGH'))
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income_binary, data=income_mdata)
    
    # voom
    voom <- voom(count, design, plot=T)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # compute average residual standard deviation and median AIC across genes 
    d_4 <- mean(fit$sigma)
    
    AIC_all <- sapply(1:nrow(fit$coefficients), function(i) {
      sigma_i <- fit$sigma[i]        # residual standard deviation for gene i
      RSS <- sigma_i^2 * (ncol(voom$E) - ncol(fit$design))     # residual sum of squares
      ncol(voom$E) * log(RSS/ncol(voom$E)) + 2*ncol(fit$design)          # AIC formula
    })
    aic_4 <- median(AIC_all)
    
    # design 5: income as binary variable (< $10,000, $10,000-$29,999, $30,000-$49,999 vs everything else)
    print('design 5')
    income_mdata$income_binary <- ifelse(income_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                                         'LOW', 'HIGH')
    income_mdata$income_binary <- factor(income_mdata$income_binary, levels=c('LOW','HIGH'))
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income_binary, data=income_mdata)
    
    # voom
    voom <- voom(count, design, plot=T)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # compute average residual standard deviation and median AIC across genes 
    d_5 <- mean(fit$sigma)
    
    AIC_all <- sapply(1:nrow(fit$coefficients), function(i) {
      sigma_i <- fit$sigma[i]        # residual standard deviation for gene i
      RSS <- sigma_i^2 * (ncol(voom$E) - ncol(fit$design))     # residual sum of squares
      ncol(voom$E) * log(RSS/ncol(voom$E)) + 2*ncol(fit$design)          # AIC formula
    })
    aic_5 <- median(AIC_all)
    
    # design 6: income as binary variable (< $10,000, $10,000-$29,999, $30,000-$49,999 vs everything else)
    print('design 6')
    income_mdata$income_binary <- ifelse(income_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', '$50,000-$69,999'),
                                         'LOW', 'HIGH')
    income_mdata$income_binary <- factor(income_mdata$income_binary, levels=c('LOW','HIGH'))
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition*income_binary, data=income_mdata)
    
    # voom
    voom <- voom(count, design, plot=T)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # compute average residual standard deviation and median AIC across genes 
    d_6 <- mean(fit$sigma)
    
    AIC_all <- sapply(1:nrow(fit$coefficients), function(i) {
      sigma_i <- fit$sigma[i]        # residual standard deviation for gene i
      RSS <- sigma_i^2 * (ncol(voom$E) - ncol(fit$design))     # residual sum of squares
      ncol(voom$E) * log(RSS/ncol(voom$E)) + 2*ncol(fit$design)          # AIC formula
    })
    aic_6 <- median(AIC_all)
    
    if (exists('mean_sigma_df')){
      mean_sigma_df <- rbind(mean_sigma_df, c(conditions[i], ctype, d_6, d_2, d_3, d_4, d_5, d_6, 
                                              aic_1, aic_2, aic_3, aic_4, aic_5, aic_6))
    } else {mean_sigma_df <- c(conditions[i], ctype, d_6, d_2, d_3, d_4, d_5, d_6, 
                               aic_1, aic_2, aic_3, aic_4, aic_5, aic_6)}
  }
}

mean_sigma_df <- as.data.frame(mean_sigma_df)
colnames(mean_sigma_df) <- c('condition', 'celltype', 'residualsd_1', 'residualsd_2', 'residualsd_3', 'residualsd_4', 
                             'residualsd_5', 'residualsd_6', 'aic_1', 'aic_2', 'aic_3', 'aic_4', 'aic_5', 'aic_6')
rownames(mean_sigma_df) <- NULL
mean_sigma_df[,3:14] <- lapply(mean_sigma_df[,3:14], as.numeric)

# change coefficients that were not estimated to NA
mean_sigma_df$residualsd_1[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='B')] <- NA
mean_sigma_df$residualsd_6[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='B')] <- NA
mean_sigma_df$residualsd_1[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='Mono')] <- NA
mean_sigma_df$residualsd_6[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='Mono')] <- NA
mean_sigma_df$residualsd_1[which(mean_sigma_df$condition=='IVA' & mean_sigma_df$celltype=='B')] <- NA
mean_sigma_df$residualsd_6[which(mean_sigma_df$condition=='IVA' & mean_sigma_df$celltype=='B')] <- NA

mean_sigma_df$aic_1[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='B')] <- NA
mean_sigma_df$aic_6[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='B')] <- NA
mean_sigma_df$aic_1[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='Mono')] <- NA
mean_sigma_df$aic_6[which(mean_sigma_df$condition=='RV' & mean_sigma_df$celltype=='Mono')] <- NA
mean_sigma_df$aic_1[which(mean_sigma_df$condition=='IVA' & mean_sigma_df$celltype=='B')] <- NA
mean_sigma_df$aic_6[which(mean_sigma_df$condition=='IVA' & mean_sigma_df$celltype=='B')] <- NA

# find best design (lowest mean sigma)
colMeans(mean_sigma_df[,3:14])
colMeans(mean_sigma_df[,3:8]) %>% which.min() # lowest residual sd
colMeans(mean_sigma_df[,9:14]) %>% which.min() # lowest aic

