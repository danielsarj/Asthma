library(Seurat)
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

# load seurat object
obj <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')

# load sample metadata and add to seurat's metadata
sample_m <- fread('../sample_metadata.txt')
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='ID'))
rownames(mdata) <- mdata$orig.ident
obj@meta.data <- mdata
rm(sample_m, mdata)

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt') %>% filter(gene_biotype=='protein_coding',
                                                                   hgnc_symbol!='', !grepl('^MT-', hgnc_symbol))
annotations$chromosome_name <- as.numeric(annotations$chromosome_name)
annotations <- annotations %>% drop_na()

# compute PCs per cell type
for (ct in unique(obj@meta.data$celltype)){
  print(ct)
  # subset seurat object
  subset_obj <- subset(obj, subset = celltype == ct)  
  
  # extract count matrix
  count_mat <- subset_obj@assays$RNA$counts
  
  # keep only protein-coding genes with variance > 0
  keep_genes <- rownames(count_mat) %in% annotations$hgnc_symbol
  count_mat <- count_mat[keep_genes, ]
  gene_vars <- apply(count_mat, 1, function(x) var(as.numeric(x)))
  count_mat <- count_mat[gene_vars > 0, ]
  rm(gene_vars, keep_genes)
  
  # convert to dense matrix for PC computing
  count_mat <- count_mat %>% as.matrix() %>% t()
  
  # compute expression PCs
  exp_pcs <- prcomp(count_mat, scale.=TRUE, center=TRUE)
  pve <- exp_pcs$sdev[1:30]^2 / sum(exp_pcs$sdev[1:30]^2) 
  exp_pcs_x <- exp_pcs$x[,1:30] %>% as.data.frame() %>% rownames_to_column('orig.ident')
  
  # add expression PCs to metadata
  pcs_mdata <- inner_join(subset_obj@meta.data, exp_pcs_x, by=('orig.ident'))
  pve <- pve %>% as.data.frame() %>% mutate(PC=seq(1:n()))
  colnames(pve)[1] <- 'PVE'
  
  # correlate PCs with covariates/confounders
  pcs_mdata$condition <- factor(pcs_mdata$condition, levels=c('NI','IVA','RV'))
  pcs_mdata$batch <- factor(pcs_mdata$batch, levels=c('B1','B2','B3','B4'))
  pcs_mdata$gender <- factor(pcs_mdata$gender, levels=c('Female', 'Male'))
  pcs_mdata$asthma <- factor(pcs_mdata$asthma, levels=c('No', 'Yes'))
  mmatrix <- model.matrix(~condition+batch+n+age+gender+avg_mt+asthma, data=pcs_mdata)[,-1]
  pcs_only <- pcs_mdata[,grep('^PC', colnames(pcs_mdata))]
  
  # create empty matrices
  cor_mat <- matrix(NA, nrow=ncol(mmatrix), ncol=ncol(pcs_only))
  pval_mat <- matrix(NA, nrow=ncol(mmatrix), ncol=ncol(pcs_only))
  rownames(cor_mat) <- rownames(pval_mat) <- colnames(mmatrix)
  colnames(cor_mat) <- colnames(pval_mat) <- colnames(pcs_only)
  
  # fill with Spearman rho and p-values
  for (i in seq_len(ncol(mmatrix))) {
    for (j in seq_len(ncol(pcs_only))) {
      test <- cor.test(mmatrix[,i], pcs_only[,j], method='spearman')
      cor_mat[i, j] <- test$estimate
      pval_mat[i, j] <- test$p.value
    }
  }
  pval_mat_adj <- matrix(p.adjust(pval_mat, method='fdr'), nrow=nrow(pval_mat),
                         dimnames=dimnames(pval_mat))
  cor_df <- as.data.frame(as.table(cor_mat)) %>% rename(Variable=Var1, PC=Var2, Spearman_rho=Freq)
  pval_df <- as.data.frame(as.table(pval_mat)) %>% rename(Variable=Var1, PC=Var2, P_value=Freq)
  cor_full <- left_join(cor_df, pval_df, by=c('Variable', 'PC'))
  cor_full$PC <- sub('^PC', '', cor_full$PC)
  cor_full$PC <- as.numeric(cor_full$PC)
  
  # add stars based on pval
  cor_full <- cor_full %>%
    mutate(sig = case_when(
      P_value <= 0.001 ~ '***',
      P_value <= 0.01  ~ '**',
      P_value <= 0.05  ~ '*',
      TRUE             ~ ''))
  
  # plot
  pc1pc2 <- ggplot(pcs_mdata, aes(x=PC1, y=PC2, color=condition, shape=asthma)) + geom_point() +
    theme_bw() + ggtitle(ct)
  pc2pc3 <- ggplot(pcs_mdata, aes(x=PC3, y=PC2, color=condition, shape=asthma)) + geom_point() +
    theme_bw() + ggtitle(ct)
  elbowp <- ggplot(pve, aes(x=PC, y=PVE)) + geom_line() + theme_bw() + ggtitle(ct) +
    scale_x_discrete(limits=as.character(1:30), breaks = c('1','5','10','15','20','25','30'))
  hmap <- ggplot(cor_full, aes(x=PC, y=Variable, fill=Spearman_rho)) + geom_tile() + theme_bw() +
    scale_fill_gradient2(low='blue',mid='white',high='red') + 
    scale_x_discrete(limits=as.character(1:30), breaks = c('1','5','10','15','20','25','30')) + 
  geom_text(aes(label=sig), color='black', size=2) + ggtitle(ct)

  (pc1pc2 + elbowp) / (pc2pc3 + hmap)
  ggsave(ct%&%'_expPCs_corr.pdf', height=6, width=12)
}
