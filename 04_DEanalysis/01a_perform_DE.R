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
    
    # define design matrix
    design <- model.matrix(~batch+age+gender+n+condition, data=mdata)
    
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
