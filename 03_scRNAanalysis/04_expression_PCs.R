library(Seurat)
library(tidyverse)
library(edgeR)
library(limma)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

# load seurat object
obj <- readRDS('NI_IVA_RV.integrated.pseudobulks_new.rds')

# load sample metadata (age, sex, asthma, income, albuterol) and add to seurat's metadata
sample_m <- fread('../sample_metadata.txt')
sample_m$income <- ifelse(sample_m$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'), 'Low', 'High')
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='ID'))
rownames(mdata) <- mdata$orig.ident
obj@meta.data <- mdata
rm(sample_m, mdata)

# load sample metadata (act, ace, resilience, social support, racism, albuterol) and add to seurat's metadata
sample_m <- fread('../Sample_test_scores_Araujo.tsv') %>% 
  mutate(ACE_result = coalesce(ACE_Percent_Self, ACE_Percent_Parent),
         ACT_conclusion = ifelse(Recorded_Diagnosis == 'No_Diagnosis', 'non_asthmatic', ACT_conclusion)) %>%
  select(Study_ID, Recorded_Diagnosis, ACT_score, ACT_conclusion, ACE_result, Parent_Resilience_Score, 
         Parents_Score_Avg, Parent_support_conclusion, Total_Racist_Events,
         Racist_stress, Racism_child_24hr, Experience_Discrimination_child)
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='Study_ID'))
rownames(mdata) <- mdata$orig.ident
obj@meta.data <- mdata
rm(sample_m, mdata)

# load sample metadata (infection at time of collectiong) and add to seurat's metadata
sample_m <- fread('../SEA_Metadata_Pathogen_Araujo.tsv') %>% 
  rename(infection_agent = Results, infection_status = Comment)
sample_m$infection_status <- gsub('infection', 'Positive', sample_m$infection_status)
mdata <- inner_join(obj@meta.data, sample_m, by=c('IDs'='Study.ID.'))
rownames(mdata) <- mdata$orig.ident
obj@meta.data <- mdata
rm(sample_m, mdata)

# update object with new metadata
saveRDS(obj, file='NI_IVA_RV.integrated.pseudobulks_new.rds')

# remove batch 4
obj <- subset(obj, subset= batch!='B4')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt') %>% filter(gene_biotype=='protein_coding',
                                                                   hgnc_symbol!='', !grepl('^MT-', hgnc_symbol))
annotations$chromosome_name <- as.numeric(annotations$chromosome_name)
annotations <- annotations %>% drop_na()

# compute PCs per cell type
for (ct in c(unique(obj@meta.data$celltype), 'PBMC')){
  print(ct)
  if (ct=='PBMC'){
    # create a PBMC-like pseudobulk
    subset_obj <- AggregateExpression(obj, group.by=c('IDs','condition'),
                                      slot='counts', assays='RNA', return.seurat=T)
    # pbmc-like metadata
    pbmc_m <- obj@meta.data %>% group_by(IDs, condition) %>% reframe(batch=batch, n=sum(n), avg_mt=mean(avg_mt),
                                 age=age, gender=gender, asthma=asthma, income=income, Recorded_Diagnosis=Recorded_Diagnosis,
                                 ACT_conclusion=ACT_conclusion, ACE_result=ACE_result, Parent_Resilience_Score=Parent_Resilience_Score, 
                                 Parents_Score_Avg=Parents_Score_Avg, Parent_support_conclusion=Parent_support_conclusion, 
                                 Total_Racist_Events=Total_Racist_Events, Racist_stress=Racist_stress, Racism_child_24hr=Racism_child_24hr,
                                 Experience_Discrimination_child=Experience_Discrimination_child, infection_status) %>% unique()

    pbmc_m <- inner_join(subset_obj@meta.data, pbmc_m, by=c('IDs', 'condition'))
    rownames(pbmc_m) <- pbmc_m$orig.ident
    subset_obj@meta.data <- pbmc_m
    rm(pbmc_m)
  } else {
    # subset seurat object
    subset_obj <- subset(obj, subset = celltype == ct)  
  }
  
  # extract count matrix
  count_mat <- subset_obj@assays$RNA$counts
  
  # keep only protein-coding genes with variance > 0
  keep_genes <- rownames(count_mat) %in% annotations$hgnc_symbol
  count_mat <- count_mat[keep_genes, ]
  gene_vars <- apply(count_mat, 1, function(x) var(as.numeric(x)))
  count_mat <- count_mat[gene_vars > 0, ]
  rm(gene_vars, keep_genes)
  
  # log transform
  dge <- DGEList(counts = count_mat)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design = NULL)

  # compute expression PCs
  exp_pcs <- prcomp(t(v$E), scale.=TRUE)
  pve <- exp_pcs$sdev[1:30]^2 / sum(exp_pcs$sdev[1:30]^2) 
  exp_pcs_x <- exp_pcs$x[,1:30] %>% as.data.frame() %>% rownames_to_column('orig.ident')
  
  # add expression PCs to metadata
  pcs_mdata <- inner_join(subset_obj@meta.data, exp_pcs_x, by=('orig.ident'))
  pve <- pve %>% as.data.frame() %>% mutate(PC=seq(1:n()))
  colnames(pve)[1] <- 'PVE'
  
  # correlate PCs with covariates/confounders
  pcs_mdata$condition <- factor(pcs_mdata$condition, levels=c('NI','IVA','RV'))
  pcs_mdata$batch <- factor(pcs_mdata$batch, levels=c('B1','B2','B3'))
  pcs_mdata$gender <- factor(pcs_mdata$gender, levels=c('Female', 'Male'))
  pcs_mdata$asthma <- factor(pcs_mdata$asthma, levels=c('No', 'Yes'))
  pcs_mdata$income <- factor(pcs_mdata$income, levels=c('Low','High'))
  pcs_mdata$Recorded_Diagnosis <- factor(pcs_mdata$Recorded_Diagnosis, levels=c('No_Diagnosis', 'Recorded_Asthma_Diagnosis'))
  pcs_mdata$ACT_conclusion <- factor(pcs_mdata$ACT_conclusion, levels=c('non_asthmatic', 'controlled', 'uncontrolled'))
  pcs_mdata$Parent_support_conclusion <- factor(pcs_mdata$Parent_support_conclusion, levels=c('low_support', 'moderate_support', 'high_support'))
  pcs_mdata$infection_status <- factor(pcs_mdata$infection_status, levels=c('Negative', 'Positive'))
  
  # take care of NAs
  num_vars <- names(pcs_mdata)[sapply(pcs_mdata, is.numeric)]
  pcs_mdata[num_vars] <- lapply(
    pcs_mdata[num_vars],
    function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
  factor_vars <- c('ACT_conclusion', 'Parent_support_conclusion', 'infection_status')
  pcs_mdata <- pcs_mdata %>% mutate(across(all_of(factor_vars), ~ .x %>% as.factor() %>%
                    fct_na_value_to_level(level = 'Missing') %>%
                    fct_relevel('Missing', after = 0)))
  
  if (ct=='PBMC'){
    mmatrix <- model.matrix(~condition+batch+n+age+gender+avg_mt+income+Recorded_Diagnosis+ACT_conclusion+
                              ACE_result+Parent_Resilience_Score+Parents_Score_Avg+Parent_support_conclusion+
                              Total_Racist_Events+Racist_stress+Racism_child_24hr+Experience_Discrimination_child+
                              infection_status, data=pcs_mdata)[,-1]
  } else {
    mmatrix <- model.matrix(~~condition+batch+n+prop+age+gender+avg_mt+income+Recorded_Diagnosis+ACT_conclusion+
                              ACE_result+Parent_Resilience_Score+Parents_Score_Avg+Parent_support_conclusion+
                              Total_Racist_Events+Racist_stress+Racism_child_24hr+Experience_Discrimination_child+
                              infection_status, data=pcs_mdata)[,-1]
  }
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
  pc1pc2 <- ggplot(pcs_mdata, aes(x=PC1, y=PC2, color=condition)) + geom_point() +
    theme_bw() + ggtitle(ct)
  pc1pc3 <- ggplot(pcs_mdata, aes(x=PC1, y=PC3, color=condition)) + geom_point() +
    theme_bw() + ggtitle(ct)
  elbowp <- ggplot(pve, aes(x=PC, y=PVE)) + geom_line() + theme_bw() + ggtitle(ct) +
    scale_x_discrete(limits=as.character(1:30), breaks = c('1','5','10','15','20','25','30'))
  hmap <- ggplot(cor_full, aes(x=PC, y=Variable, fill=Spearman_rho)) + geom_tile() + theme_bw() +
    scale_fill_gradient2(low='blue',mid='white',high='red') + 
    scale_x_discrete(limits=as.character(1:30), breaks = c('1','5','10','15','20','25','30')) + 
  geom_text(aes(label=sig), color='black', size=2) + ggtitle(ct)

  (pc1pc2 + pc1pc3) / (elbowp + hmap)
  ggsave(ct%&%'_expPCs_corr_new.pdf', height=6, width=12)
}
