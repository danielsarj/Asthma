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

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
  annotations$hgnc_symbol!='' &
  !grepl('^MT-', annotations$hgnc_symbol)]

# load seurat object
objs <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')

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
    count <- tmp@assays$RNA$counts
    
    # filter count matrix (only keep protein coding genes)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% annotations,]
    zero_var_genes <- apply(count, 1, var) == 0
    count <- count[!zero_var_genes, ]
    
    # transform count into dge object
    count <- DGEList(counts=count)
    count <- calcNormFactors(count)
    
    # define design matrix
    design <- model.matrix(~0+IDs+condition, data=mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # get results
    results <- topTable(fit, coef=ncol(fit), number=Inf)
    fwrite(results, '../DEanalysis/NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_results.txt',
           sep=' ', col.names=T, row.names=T)
  }
}
