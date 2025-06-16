library(Seurat)
library(SeuratData)
library(patchwork)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(data.table)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('RV', 'IVA')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
  annotations$hgnc_symbol!='' &
  !grepl('^MT-', annotations$hgnc_symbol)]

# load seurat object
objs <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')

# merge metadata
mdata <- objs@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
objs@meta.data <- mdata

for (i in 1:length(conditions)){
  print(c(conditions[i]))

  # celltype specific DE
  for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
    print(ctype)
    
    # subset pseudobulk object
    tmp <- subset(objs, celltype==ctype & (condition==conditions[i] | condition=='NI'))
    
    # extract metadata and count matrices
    mdata <- tmp@meta.data
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$IDs <- as.factor(mdata$IDs)
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))
    count <- tmp@assays$RNA$counts
    
    # filter count matrix (only keep protein coding genes)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% annotations,]
    zero_var_genes <- apply(count, 1, var) == 0
    count <- count[!zero_var_genes, ]
    
    # transform count into dge object
    count <- DGEList(counts=count)
    count <- calcNormFactors(count)
    
    # design matrix to remove batch effects
    design <- model.matrix(~0+batch, data=mdata)
    design <- design[,colSums(design!=0)>0]
    
    # apply mean variance weights and design matrix to phenotype data
    voom <- voom(count, design, plot=FALSE)
    fit <- lmFit(voom, design)
    fit <- eBayes(fit)
    
    # get residuals to regress out batch effect
    residuals <- residuals.MArrayLM(object=fit, voom)
    
    # calculate average batch effect for each gene
    avg_batch_effect <- rowMeans(fit$coefficients)
    
    # add the average batch effect back into residuals
    corrected_expression <- apply(residuals,2,function(x){x+avg_batch_effect})
    weights <- voom$weights
    colnames(weights) <- colnames(corrected_expression)
    rownames(weights) <- rownames(corrected_expression)
    
    # identify rows with no negative values and subset matrices
    keep_rows <- rowSums(corrected_expression<0)==0
    corrected_expression <- corrected_expression[keep_rows, ]
    weights <- weights[keep_rows, ]
    
    # save files
    fwrite(corrected_expression, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_corrected_expression.txt', 
           quote=F, sep=' ', col.names=T, row.names=T)
    fwrite(weights, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_weights.txt', quote=F,
           sep=' ', col.names=T, row.names=T)
    
    # model infection differential expression
    design <- model.matrix(~0+age+gender+n+condition, data=mdata)
    
    # fit linear model voom
    voom <- voom(corrected_expression, weights=weights, design, plot=F)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # get results
    results <- topTable(fit, coef=ncol(fit), number=Inf)
    fwrite(results, '../DEanalysis/NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_results_mod.txt',
           sep=' ', col.names=T, row.names=T)
  }
}
