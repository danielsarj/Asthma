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

for (i in 1:length(conditions)){
  print(c(conditions[i]))
  
  objs <- readRDS('NI_'%&%conditions[i]%&%'.integrated.w_celltype.rds')
  
  # clean and organize metadata
  meta.data <- objs@meta.data
  meta.data$condition <- factor(meta.data$condition, levels=c('NI', conditions[i]))
  meta.data$batch <- as.factor(meta.data$batch)
  meta.data$IDs <- as.factor(meta.data$IDs)
  meta.data$predicted.celltype.l1 <- gsub(' ', '-', meta.data$predicted.celltype.l1)
  meta.data$predicted.celltype.l2 <- gsub(' ', '-', meta.data$predicted.celltype.l2)
  objs@meta.data <- meta.data
  
  # compute number of cells per condition/IDs
  cells_summary <- meta.data %>% group_by(IDs,predicted.celltype.l1,condition) %>% summarise(n_cells=n())
  cells_summary$predicted.celltype.l1 <- gsub('_', '-', cells_summary$predicted.celltype.l1)
  
  # aggregate cell types per individual per condition (pseudobulk)
  bulk_objs <- AggregateExpression(objs, group.by=c('IDs','predicted.celltype.l1','condition'), 
                                   slot='counts', assays='RNA', return.seurat=T)
  
  # update pseudobulk metadata with number of cells
  meta.data <- bulk_objs@meta.data %>% inner_join(cells_summary)
  rownames(meta.data) <- meta.data$orig.ident
  bulk_objs@meta.data <- meta.data
  
  saveRDS(bulk_objs, file='../DEanalysis/NI_'%&%conditions[i]%&%'_pseudobulks.rds')

  # celltype specific DE
  for (ctype in c('B','CD4-T','CD8-T','Mono','NK')){
    print(ctype)
    # subset pseudobulk object
    tmp <- subset(bulk_objs, predicted.celltype.l1==ctype)
    
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
    voom <- voom(count, design, plot=T)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # get results
    results <- topTable(fit, coef=ncol(fit), number=Inf)
    fwrite(results, '../DEanalysis/NI_'%&%conditions[i]%&%'_limma_'%&%ctype%&%'_results.txt',
           sep=' ', col.names=T, row.names=T)
  }
}