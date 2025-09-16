library(Seurat)
library(msigdbr)
library(data.table)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(edgeR)
library(limma)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA', 'NI')
celltypes <- c('B','T-CD4','T-CD8','Mono','NK')

### important things to load
# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# get human hallmark gene sets
ifn_genes <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)

# filter for IFN modules 
ifn_genes <- ifn_genes[grep('INTERFERON', names(ifn_genes), value=T)]

# fgsea results
fgsea_full <- fread('NI_IVAxRV_integrated_descGSEAresults.txt')
###

## overall single-cell pattern of IFN score
# load sc seurat object
sc_obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.w_celltype.rds') %>% 
  SetIdent(value='celltype')
sc_obj@meta.data$condition <- factor(sc_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# add module score
sc_obj <- AddModuleScore(sc_obj, features=ifn_genes[1], name='IFNa_response')
sc_obj <- AddModuleScore(sc_obj, features=ifn_genes[2], name='IFNy_response')

# viz
FeaturePlot(sc_obj, features='IFNa_response1', label=TRUE, repel=TRUE, split.by='condition') /
FeaturePlot(sc_obj, features='IFNy_response1', label=TRUE, repel=TRUE, split.by='condition')
ggsave('UMAP_IFNresponse.pdf', height=6, width=10)

# viz
merged_mdata <- left_join(sc_obj@meta.data, sample_m, by=c('IDs'='ID'))
rownames(merged_mdata) <- rownames(sc_obj@meta.data)
sc_obj@meta.data <- merged_mdata
sc_obj$cond_asthma <- paste(sc_obj$condition, sc_obj$asthma, sep='_')
sc_obj@meta.data$cond_asthma <- factor(sc_obj@meta.data$cond_asthma, levels=c('NI_No', 'IVA_No', 'RV_No', 'NI_Yes', 'IVA_Yes', 'RV_Yes'))

FeaturePlot(sc_obj, features='IFNa_response1', label=TRUE, repel=TRUE, split.by='cond_asthma')  +
  patchwork::plot_layout(ncol=3, nrow=2)
ggsave('UMAP_IFNaresponse_byAsthma.pdf', height=6, width=10)
FeaturePlot(sc_obj, features='IFNy_response1', label=TRUE, repel=TRUE, split.by='cond_asthma')  +
  patchwork::plot_layout(ncol=3, nrow=2)
ggsave('UMAP_IFNyresponse_byAsthma.pdf', height=6, width=10)
rm(sc_obj, merged_mdata)

## pseudobulk level IFN scores
# load pseudobulk seurat object
bulk_obj <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds') 
bulk_obj@meta.data$condition <- factor(bulk_obj@meta.data$condition, levels=c('NI','IVA','RV'))

# merge metadata
mdata <- bulk_obj@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
bulk_obj@meta.data <- mdata

# paired-level, all pathway genes, no interaction yet
for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV', 'IVA')){
    print(cond)
    
    # subset
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
      
    # subset metadata
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
    sub_mdata <- sub_mdata %>% group_by(IDs) %>% mutate(is_unique = n() == 1) %>% ungroup() %>% 
      filter(is_unique==FALSE) %>% select(-is_unique)
    
    # extract and subset count matrix
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% sub_mdata$samples]
    count <- DGEList(counts=count) %>% calcNormFactors()

    # define design matrix (without condition term)
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=sub_mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # remove confounding effects
    sub_mdata$condition <- as.numeric(sub_mdata$condition)
    sub_mdata$gender <- as.numeric(sub_mdata$gender)
    adj_expr <- removeBatchEffect(voom$E,
                                  covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender')]),
                                  batch = sub_mdata$batch,
                                  design = model.matrix(~condition, data=sub_mdata))
    
    
    ## first, IFNa
    # subset to IFNa genes
    ifn_a <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[1]]) %>% 
      as.matrix() %>% as.data.frame() %>% rownames_to_column('Gene')

    # second, IFNy
    # subset to IFNy genes
    ifn_y <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[2]]) %>% 
      as.matrix() %>% as.data.frame() %>% rownames_to_column('Gene')

   for (id in unique(sub_mdata$IDs)){
     
      # ifna
      sub_ifn_a_scores <- ifn_a %>% select(Gene, contains(id)) %>% select(Gene, contains(cond), contains('NI'))
      sub_ifn_a_scores <- (sub_ifn_a_scores[,2] - sub_ifn_a_scores[,3]) %>% as.data.frame() %>% mutate(donor=id, IFNa_score=mean(.))
      sub_ifn_a_scores <- sub_ifn_a_scores %>% select(donor, IFNa_score) %>% unique()
      
      # ifna
      sub_ifn_y_scores <- ifn_y %>% select(Gene, contains(id)) %>% select(Gene, contains(cond), contains('NI'))
      sub_ifn_y_scores <- (sub_ifn_y_scores[,2] - sub_ifn_y_scores[,3]) %>% as.data.frame() %>% mutate(donor=id, IFNy_score=mean(.))
      sub_ifn_y_scores <- sub_ifn_y_scores %>% select(donor, IFNy_score) %>% unique()
      
      # join
      sub_ifn_scores <- inner_join(sub_ifn_a_scores, sub_ifn_y_scores, by=c('donor')) %>% mutate(celltype=ctype, condition=cond)
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, sub_ifn_scores)
      } else {paired.bulk.ifn.scores <- sub_ifn_scores}
    }
  }
}
# plot IFN scores
(ggplot(paired.bulk.ifn.scores, aes(x=celltype, y=IFNa_score, fill=condition)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
               alpha=0.4, size=1) + theme_bw()) /
  (ggplot(paired.bulk.ifn.scores, aes(x=celltype, y=IFNy_score, fill=condition)) + 
     geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
     geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
                alpha=0.4, size=1) + theme_bw())
ggsave('IFNscores_paired_adjusted_boxplots.pdf', height=6, width=10)
rm(paired.bulk.ifn.scores)

# paired-level, leading edge genes, no interaction yet
for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV', 'IVA')){
    print(cond)
    
    # subset
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # subset metadata
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples) %>% rownames_to_column('samples')
    sub_mdata <- sub_mdata %>% group_by(IDs) %>% mutate(is_unique = n() == 1) %>% ungroup() %>% 
      filter(is_unique==FALSE) %>% select(-is_unique)
    
    # extract and subset count matrix
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% sub_mdata$samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # define design matrix (without condition term)
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=sub_mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # remove confounding effects
    sub_mdata$condition <- as.numeric(sub_mdata$condition)
    sub_mdata$gender <- as.numeric(sub_mdata$gender)
    adj_expr <- removeBatchEffect(voom$E,
                                  covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender')]),
                                  batch = sub_mdata$batch,
                                  design = model.matrix(~condition, data=sub_mdata))
    
    # extract ifn leading edge genes
    ifn_a_edge <- fgsea_full %>% filter(celltype==ctype, condition==cond, interaction=='none', 
                                       pathway=='INTERFERON_ALPHA_RESPONSE') %>% pull(leadingEdge) %>%
      unlist() %>% strsplit('\\|') %>% unlist() %>% gsub('"', '', .) %>% trimws()
    ifn_y_edge <- fgsea_full %>% filter(celltype==ctype, condition==cond, interaction=='none', 
                                        pathway=='INTERFERON_GAMMA_RESPONSE') %>% pull(leadingEdge) %>%
      unlist() %>% strsplit('\\|') %>% unlist() %>% gsub('"', '', .) %>% trimws()
    
    ## first, IFNa
    # subset to IFNa genes
    ifn_a <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_a_edge) %>% 
      as.matrix() %>% as.data.frame() %>% rownames_to_column('Gene')
    
    # second, IFNy
    # subset to IFNy genes
    ifn_y <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_y_edge) %>% 
      as.matrix() %>% as.data.frame() %>% rownames_to_column('Gene')
    
    for (id in unique(sub_mdata$IDs)){
      
      # ifna
      sub_ifn_a_scores <- ifn_a %>% select(Gene, contains(id)) %>% select(Gene, contains(cond), contains('NI'))
      sub_ifn_a_scores <- (sub_ifn_a_scores[,2] - sub_ifn_a_scores[,3]) %>% as.data.frame() %>% mutate(donor=id, IFNa_score=mean(.))
      sub_ifn_a_scores <- sub_ifn_a_scores %>% select(donor, IFNa_score) %>% unique()
      
      # ifna
      sub_ifn_y_scores <- ifn_y %>% select(Gene, contains(id)) %>% select(Gene, contains(cond), contains('NI'))
      sub_ifn_y_scores <- (sub_ifn_y_scores[,2] - sub_ifn_y_scores[,3]) %>% as.data.frame() %>% mutate(donor=id, IFNy_score=mean(.))
      sub_ifn_y_scores <- sub_ifn_y_scores %>% select(donor, IFNy_score) %>% unique()
      
      # join
      sub_ifn_scores <- inner_join(sub_ifn_a_scores, sub_ifn_y_scores, by=c('donor')) %>% mutate(celltype=ctype, condition=cond)
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, sub_ifn_scores)
      } else {paired.bulk.ifn.scores <- sub_ifn_scores}
    }
  }
}
# plot IFN scores
(ggplot(paired.bulk.ifn.scores, aes(x=celltype, y=IFNa_score, fill=condition)) + 
    geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
    geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
               alpha=0.4, size=1) + theme_bw()) /
  (ggplot(paired.bulk.ifn.scores, aes(x=celltype, y=IFNy_score, fill=condition)) + 
     geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
     geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
                alpha=0.4, size=1) + theme_bw())
ggsave('IFNscores_paired_adjusted_leadingedge_boxplots.pdf', height=6, width=10)
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% mutate(interaction='none', leading_edge_score=TRUE)
full_df <- rbind(full_df, paired.bulk.ifn.scores)
rm(paired.bulk.ifn.scores)

# paired-level, all pathway genes, asthma vs infection interaction
for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)
    
    # subset
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$albuterol <- factor(sub_mdata$albuterol, levels=c('No', 'Yes'))
    sub_mdata$asthma <- factor(sub_mdata$asthma, levels=c('No', 'Yes'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_asthma_alb_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # extract and subset count matrix
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% v_samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # define design matrix (without asthma+interaction term)
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition, data=sub_mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # remove confounding effects
    sub_mdata$condition <- as.numeric(sub_mdata$condition)
    sub_mdata$gender <- as.numeric(sub_mdata$gender)
    sub_mdata$albuterol <- as.numeric(sub_mdata$albuterol) 
    sub_mdata$asthma <- as.numeric(sub_mdata$asthma)
    adj_expr <- removeBatchEffect(voom$E,
                                  covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender','albuterol','condition')]),
                                  batch = sub_mdata$batch,
                                  design = model.matrix(~asthma+condition:asthma, data=sub_mdata))
    
    ## first, IFNa
    # subset to IFNa genes
    ifn_a <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[1]]) %>% as.matrix()
    
    if (exists('ifn_sizes')){
      ifn_sizes <- rbind(ifn_sizes, c(ctype, cond, 'asthma', 'IFNa', nrow(ifn_a)))
    } else {ifn_sizes <- c(ctype, cond, 'asthma', 'IFNa', nrow(ifn_a))}
    
    # scale per gene and compute score
    ifn_a <- ifn_a %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_a)[2] <- 'IFNa_tmp_score'
    
    # second, IFNy
    # subset to IFNa genes
    ifn_y <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[2]]) %>% as.matrix()
    ifn_sizes <- rbind(ifn_sizes, c(ctype, cond, 'asthma', 'IFNy', nrow(ifn_y)))

    # scale per gene and compute score
    ifn_y <- ifn_y %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_y)[2] <- 'IFNy_tmp_score'
    
    # join and reshape dfs
    ifn_scores <- inner_join(ifn_a, ifn_y, by=c('ID')) %>% 
      separate(col=ID, into=c('donor', 'condition', 'celltype'), sep='_', remove=TRUE)
    
    # remove unique indvs
    ifn_scores <- ifn_scores %>% group_by(donor) %>% mutate(is_unique = n() == 1)
    ifn_scores <- ifn_scores %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(ifn_scores$donor)){
      # select indv
      sub_ifn_scores <- ifn_scores %>% filter(donor==id) 
      
      # confirm position of NI vs ifnection
      NI_row <- which(sub_ifn_scores$condition=='NI')
      Inf_row <- which(sub_ifn_scores$condition==cond)
      
      # compute ifn score
      sub_ifn_scores$IFNa_score <- sub_ifn_scores$IFNa_tmp_score[Inf_row] - sub_ifn_scores$IFNa_tmp_score[NI_row]
      sub_ifn_scores$IFNy_score <- sub_ifn_scores$IFNy_tmp_score[Inf_row] - sub_ifn_scores$IFNy_tmp_score[NI_row]
      
      # edit df
      sub_ifn_scores <- sub_ifn_scores %>% filter(condition==cond) %>% 
        select(donor, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, sub_ifn_scores)
      } else {paired.bulk.ifn.scores <- sub_ifn_scores}
    }
  }
}
# merge w metadata to retrieve asthma status
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% inner_join(sample_m, by=c('donor'='ID'))
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% select(celltype, condition, IFNa_score, IFNy_score, asthma) %>% 
  pivot_longer(cols=c(IFNa_score, IFNy_score), names_to='IFN', values_to='score')
paired.bulk.ifn.scores$IFN <- gsub('_score','',paired.bulk.ifn.scores$IFN)

# plot splitting by asthma status
ggplot(paired.bulk.ifn.scores, aes(x=condition, y=score, fill=asthma)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(celltype), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('IFNscores_asthma_paired_adjusted_noasthma_boxplots.pdf', height=4, width=10)

paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>%
  group_by(IFN, condition, celltype) %>%
  summarise(t_test_res = list(broom::tidy(t.test(score ~ asthma))),
    .groups = "drop") %>% tidyr::unnest(t_test_res) %>% 
  select(IFN, condition, celltype, p.value) %>% mutate(interaction='asthma', is_leading_edge=FALSE)
full_df <- paired.bulk.ifn.scores
rm(paired.bulk.ifn.scores)

# paired-level, leading edge genes, asthma vs infection interaction
for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)
    
    # subset
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$albuterol <- factor(sub_mdata$albuterol, levels=c('No', 'Yes'))
    sub_mdata$asthma <- factor(sub_mdata$asthma, levels=c('No', 'Yes'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_asthma_alb_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # extract and subset count matrix
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% v_samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # define design matrix (without asthma+interaction term)
    design <- model.matrix(~batch+age+gender+n+avg_mt+albuterol+condition, data=sub_mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # remove confounding effects
    sub_mdata$condition <- as.numeric(sub_mdata$condition)
    sub_mdata$gender <- as.numeric(sub_mdata$gender)
    sub_mdata$albuterol <- as.numeric(sub_mdata$albuterol) 
    sub_mdata$asthma <- as.numeric(sub_mdata$asthma)
    adj_expr <- removeBatchEffect(voom$E,
                                  covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender','albuterol','condition')]),
                                  batch = sub_mdata$batch,
                                  design = model.matrix(~asthma+condition:asthma, data=sub_mdata))
    
    # extract ifn leading edge genes
    ifn_a_edge <- fgsea_full %>% filter(celltype==ctype, condition==cond, interaction=='asthma', 
                                        pathway=='INTERFERON_ALPHA_RESPONSE') %>% pull(leadingEdge) %>%
      unlist() %>% strsplit('\\|') %>% unlist() %>% gsub('"', '', .) %>% trimws()
    ifn_y_edge <- fgsea_full %>% filter(celltype==ctype, condition==cond, interaction=='asthma', 
                                        pathway=='INTERFERON_GAMMA_RESPONSE') %>% pull(leadingEdge) %>%
      unlist() %>% strsplit('\\|') %>% unlist() %>% gsub('"', '', .) %>% trimws()
    
    ## first, IFNa
    # subset to IFNa genes
    ifn_a <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[1]]) %>% as.matrix()
    
    # scale per gene and compute score
    ifn_a <- ifn_a %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% as.data.frame() %>% 
      filter(rownames(.) %in% ifn_a_edge) %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_a)[2] <- 'IFNa_tmp_score'
    
    # second, IFNy
    # subset to IFNa genes
    ifn_y <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[2]]) %>% as.matrix()
    
    # scale per gene and compute score
    ifn_y <- ifn_y %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% as.data.frame() %>% 
      filter(rownames(.) %in% ifn_y_edge) %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_y)[2] <- 'IFNy_tmp_score'
    
    # join and reshape dfs
    ifn_scores <- inner_join(ifn_a, ifn_y, by=c('ID')) %>% 
      separate(col=ID, into=c('donor', 'condition', 'celltype'), sep='_', remove=TRUE)
    
    # remove unique indvs
    ifn_scores <- ifn_scores %>% group_by(donor) %>% mutate(is_unique = n() == 1)
    ifn_scores <- ifn_scores %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(ifn_scores$donor)){
      # select indv
      sub_ifn_scores <- ifn_scores %>% filter(donor==id) 
      
      # confirm position of NI vs ifnection
      NI_row <- which(sub_ifn_scores$condition=='NI')
      Inf_row <- which(sub_ifn_scores$condition==cond)
      
      # compute ifn score
      sub_ifn_scores$IFNa_score <- sub_ifn_scores$IFNa_tmp_score[Inf_row] - sub_ifn_scores$IFNa_tmp_score[NI_row]
      sub_ifn_scores$IFNy_score <- sub_ifn_scores$IFNy_tmp_score[Inf_row] - sub_ifn_scores$IFNy_tmp_score[NI_row]
      
      # edit df
      sub_ifn_scores <- sub_ifn_scores %>% filter(condition==cond) %>% 
        select(donor, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, sub_ifn_scores)
      } else {paired.bulk.ifn.scores <- sub_ifn_scores}
    }
  }
}
# merge w metadata to retrieve asthma status
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% inner_join(sample_m, by=c('donor'='ID'))
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% select(celltype, condition, IFNa_score, IFNy_score, asthma) %>% 
  pivot_longer(cols=c(IFNa_score, IFNy_score), names_to='IFN', values_to='score')
paired.bulk.ifn.scores$IFN <- gsub('_score','',paired.bulk.ifn.scores$IFN)

# plot splitting by asthma status
ggplot(paired.bulk.ifn.scores, aes(x=condition, y=score, fill=asthma)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(celltype), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = asthma), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('IFNscores_asthma_paired_adjusted_leadingedge_noasthma_boxplots.pdf', height=4, width=10)

paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>%
  group_by(IFN, condition, celltype) %>%
  summarise(t_test_res = list(broom::tidy(t.test(score ~ asthma))),
            .groups = "drop") %>% tidyr::unnest(t_test_res) %>% 
  select(IFN, condition, celltype, p.value) %>% mutate(interaction='asthma', is_leading_edge=TRUE)
full_df <- rbind(full_df, paired.bulk.ifn.scores)
rm(paired.bulk.ifn.scores)

# paired-level, all pathway genes, income vs infection interaction
for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)
    
    # subset
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$income <- ifelse(sub_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                           'Low', 'High')
    sub_mdata$income <- factor(sub_mdata$income, levels=c('Low','High'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_income_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # extract and subset count matrix
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% v_samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # define design matrix (without income+interaction term)
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=sub_mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # remove confounding effects
    sub_mdata$condition <- as.numeric(sub_mdata$condition)
    sub_mdata$gender <- as.numeric(sub_mdata$gender)
    sub_mdata$income <- as.numeric(sub_mdata$income)
    adj_expr <- removeBatchEffect(voom$E,
                                  covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender','condition')]),
                                  batch = sub_mdata$batch,
                                  design = model.matrix(~income+condition:income, data=sub_mdata))
    
    ## first, IFNa
    # subset to IFNa genes
    ifn_a <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[1]]) %>% as.matrix()
    ifn_sizes <- rbind(ifn_sizes, c(ctype, cond, 'income', 'IFNa', nrow(ifn_a)))

    # scale per gene and compute score
    ifn_a <- ifn_a %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_a)[2] <- 'IFNa_tmp_score'
    
    # second, IFNy
    # subset to IFNa genes
    ifn_y <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[2]]) %>% as.matrix()
    ifn_sizes <- rbind(ifn_sizes, c(ctype, cond, 'income', 'IFNy', nrow(ifn_y)))
    
    # scale per gene and compute score
    ifn_y <- ifn_y %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_y)[2] <- 'IFNy_tmp_score'
    
    # join and reshape dfs
    ifn_scores <- inner_join(ifn_a, ifn_y, by=c('ID')) %>% 
      separate(col=ID, into=c('donor', 'condition', 'celltype'), sep='_', remove=TRUE)
    
    # remove unique indvs
    ifn_scores <- ifn_scores %>% group_by(donor) %>% mutate(is_unique = n() == 1)
    ifn_scores <- ifn_scores %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(ifn_scores$donor)){
      # select indv
      sub_ifn_scores <- ifn_scores %>% filter(donor==id) 
      
      # confirm position of NI vs ifnection
      NI_row <- which(sub_ifn_scores$condition=='NI')
      Inf_row <- which(sub_ifn_scores$condition==cond)
      
      # compute ifn score
      sub_ifn_scores$IFNa_score <- sub_ifn_scores$IFNa_tmp_score[Inf_row] - sub_ifn_scores$IFNa_tmp_score[NI_row]
      sub_ifn_scores$IFNy_score <- sub_ifn_scores$IFNy_tmp_score[Inf_row] - sub_ifn_scores$IFNy_tmp_score[NI_row]
      
      # edit df
      sub_ifn_scores <- sub_ifn_scores %>% filter(condition==cond) %>% 
        select(donor, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, sub_ifn_scores)
      } else {paired.bulk.ifn.scores <- sub_ifn_scores}
    }
  }
}
# merge w metadata to retrieve income status
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% inner_join(sample_m, by=c('donor'='ID'))
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% select(celltype, condition, IFNa_score, IFNy_score, income) %>% 
  pivot_longer(cols=c(IFNa_score, IFNy_score), names_to='IFN', values_to='score')
paired.bulk.ifn.scores$IFN <- gsub('_score','',paired.bulk.ifn.scores$IFN)
paired.bulk.ifn.scores$income <- ifelse(paired.bulk.ifn.scores$income %in% 
                                          c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'), 'Low', 'High')
paired.bulk.ifn.scores$income <- factor(paired.bulk.ifn.scores$income, levels=c('Low','High'))

# plot splitting by income status
ggplot(paired.bulk.ifn.scores, aes(x=condition, y=score, fill=income)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(celltype), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = income), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('IFNscores_income_paired_adjusted_noincome_boxplots.pdf', height=4, width=10)

paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>%
  group_by(IFN, condition, celltype) %>%
  summarise(t_test_res = list(broom::tidy(t.test(score ~ income))),
            .groups = "drop") %>% tidyr::unnest(t_test_res) %>% 
  select(IFN, condition, celltype, p.value) %>% mutate(interaction='income', is_leading_edge=FALSE)
full_df <- rbind(full_df, paired.bulk.ifn.scores)
rm(paired.bulk.ifn.scores)

# paired-level, leading edge genes, income vs infection interaction
for (ctype in celltypes){
  print(ctype)
  for (cond in c('RV','IVA')){
    print(cond)
    
    # subset
    tmp <- subset(bulk_obj, celltype==ctype & (condition=='NI' | condition==cond))
    sub_mdata <- tmp@meta.data
    sub_mdata$gender <- factor(sub_mdata$gender, levels=c('Male','Female'))
    sub_mdata$condition <- factor(sub_mdata$condition, levels=c('NI', cond))
    sub_mdata$income <- ifelse(sub_mdata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'),
                               'Low', 'High')
    sub_mdata$income <- factor(sub_mdata$income, levels=c('Low','High'))
    
    # get samples/genes from voom expression dataframe (easier to filter stuff)
    v <- fread('../scRNAanalysis/NI_'%&%cond%&%'_'%&%ctype%&%'_income_voom_expression.txt') %>% column_to_rownames('Gene')
    v_samples <- colnames(v)
    v_genes <- rownames(v)
    
    # extract and subset count matrix
    sub_mdata <- sub_mdata %>% filter(rownames(.) %in% v_samples)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% v_genes,]
    count <- count[,colnames(count) %in% v_samples]
    count <- DGEList(counts=count) %>% calcNormFactors()
    
    # define design matrix (without income+interaction term)
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=sub_mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # remove confounding effects
    sub_mdata$condition <- as.numeric(sub_mdata$condition)
    sub_mdata$gender <- as.numeric(sub_mdata$gender)
    sub_mdata$income <- as.numeric(sub_mdata$income)
    adj_expr <- removeBatchEffect(voom$E,
                                  covariates = as.matrix(sub_mdata[, c('age','n','avg_mt','gender','condition')]),
                                  batch = sub_mdata$batch,
                                  design = model.matrix(~income+condition:income, data=sub_mdata))
    
    # extract ifn leading edge genes
    ifn_a_edge <- fgsea_full %>% filter(celltype==ctype, condition==cond, interaction=='income', 
                                        pathway=='INTERFERON_ALPHA_RESPONSE') %>% pull(leadingEdge) %>%
      unlist() %>% strsplit('\\|') %>% unlist() %>% gsub('"', '', .) %>% trimws()
    ifn_y_edge <- fgsea_full %>% filter(celltype==ctype, condition==cond, interaction=='income', 
                                        pathway=='INTERFERON_GAMMA_RESPONSE') %>% pull(leadingEdge) %>%
      unlist() %>% strsplit('\\|') %>% unlist() %>% gsub('"', '', .) %>% trimws()
    
    ## first, IFNa
    # subset to IFNa genes
    ifn_a <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[1]]) %>% as.matrix()
    
    # scale per gene and compute score
    ifn_a <- ifn_a %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% as.data.frame() %>% 
      filter(rownames(.) %in% ifn_a_edge) %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_a)[2] <- 'IFNa_tmp_score'
    
    # second, IFNy
    # subset to IFNa genes
    ifn_y <- adj_expr %>% as.data.frame() %>% filter(rownames(.) %in% ifn_genes[[2]]) %>% as.matrix()
    
    # scale per gene and compute score
    ifn_y <- ifn_y %>% t() %>% scale(center=TRUE, scale=TRUE) %>% t() %>% as.data.frame() %>% 
      filter(rownames(.) %in% ifn_y_edge) %>% colMeans(na.rm=TRUE) %>% 
      as.data.frame() %>% rownames_to_column('ID')
    colnames(ifn_y)[2] <- 'IFNy_tmp_score'
    
    # join and reshape dfs
    ifn_scores <- inner_join(ifn_a, ifn_y, by=c('ID')) %>% 
      separate(col=ID, into=c('donor', 'condition', 'celltype'), sep='_', remove=TRUE)
    
    # remove unique indvs
    ifn_scores <- ifn_scores %>% group_by(donor) %>% mutate(is_unique = n() == 1)
    ifn_scores <- ifn_scores %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    
    for (id in unique(ifn_scores$donor)){
      # select indv
      sub_ifn_scores <- ifn_scores %>% filter(donor==id) 
      
      # confirm position of NI vs ifnection
      NI_row <- which(sub_ifn_scores$condition=='NI')
      Inf_row <- which(sub_ifn_scores$condition==cond)
      
      # compute ifn score
      sub_ifn_scores$IFNa_score <- sub_ifn_scores$IFNa_tmp_score[Inf_row] - sub_ifn_scores$IFNa_tmp_score[NI_row]
      sub_ifn_scores$IFNy_score <- sub_ifn_scores$IFNy_tmp_score[Inf_row] - sub_ifn_scores$IFNy_tmp_score[NI_row]
      
      # edit df
      sub_ifn_scores <- sub_ifn_scores %>% filter(condition==cond) %>% 
        select(donor, condition, celltype, IFNa_score, IFNy_score) 
      
      if (exists('paired.bulk.ifn.scores')){
        paired.bulk.ifn.scores <- rbind(paired.bulk.ifn.scores, sub_ifn_scores)
      } else {paired.bulk.ifn.scores <- sub_ifn_scores}
    }
  }
}
# merge w metadata to retrieve income status
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% inner_join(sample_m, by=c('donor'='ID'))
paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>% select(celltype, condition, IFNa_score, IFNy_score, income) %>% 
  pivot_longer(cols=c(IFNa_score, IFNy_score), names_to='IFN', values_to='score')
paired.bulk.ifn.scores$IFN <- gsub('_score','',paired.bulk.ifn.scores$IFN)
paired.bulk.ifn.scores$income <- ifelse(paired.bulk.ifn.scores$income %in% 
                                          c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'), 'Low', 'High')
paired.bulk.ifn.scores$income <- factor(paired.bulk.ifn.scores$income, levels=c('Low','High'))

# plot splitting by income status
ggplot(paired.bulk.ifn.scores, aes(x=condition, y=score, fill=income)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8)) + 
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
             alpha=0.4, size=1) + facet_grid(cols=vars(celltype), rows=vars(IFN), scale='free') + theme_bw() +
  stat_compare_means(aes(group = income), method='t.test', label='p.format') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave('IFNscores_income_paired_adjusted_leadingedge_noincome_boxplots.pdf', height=4, width=10)

paired.bulk.ifn.scores <- paired.bulk.ifn.scores %>%
  group_by(IFN, condition, celltype) %>%
  summarise(t_test_res = list(broom::tidy(t.test(score ~ income))),
            .groups = "drop") %>% tidyr::unnest(t_test_res) %>% 
  select(IFN, condition, celltype, p.value) %>% mutate(interaction='income', is_leading_edge=TRUE)
full_df <- rbind(full_df, paired.bulk.ifn.scores)
rm(paired.bulk.ifn.scores)

### analyze leading edge sets
ifn_fgsea <- fgsea_full %>% filter(pathway %in% c('INTERFERON_ALPHA_RESPONSE', 'INTERFERON_GAMMA_RESPONSE')) %>% 
  select(pathway, celltype, condition, interaction, padj, NES, leadingEdge) %>% filter(interaction!='none')
ifn_fgsea$pathway <- gsub('_RESPONSE', '', ifn_fgsea$pathway)

# split leading edge column and append its size to full df
leading_edge_list <- lapply(ifn_fgsea$leadingEdge, function(x) {
  gsub("\"", "", unlist(strsplit(x, "\"\\|\"")))})
ifn_fgsea <- ifn_fgsea %>% mutate(leadingEdge_size=sapply(leading_edge_list, length), is_sig=ifelse(padj<0.05, 1, 0.5))
ggplot(ifn_fgsea) + geom_col(aes(x=celltype, y=leadingEdge_size, fill=interaction, alpha=is_sig), position='dodge') + 
  facet_grid(rows=vars(pathway), cols=vars(condition)) + theme_bw()
ggsave('LeadingEdgeSizes_IFNpathways_barplot.pdf', height=4, width=7)

# add number of genes to ifn fsega df
ifn_sizes <- as.data.frame(ifn_sizes)
rownames(ifn_sizes) <- NULL
colnames(ifn_sizes) <- c('celltype', 'condition', 'interaction', 'IFN', 'size')
ifn_sizes$IFN <- gsub('IFNa', 'INTERFERON_ALPHA', ifn_sizes$IFN)
ifn_sizes$IFN <- gsub('IFNy', 'INTERFERON_GAMMA', ifn_sizes$IFN)
ifn_sizes$size <- as.numeric(ifn_sizes$size)

ifn_fgsea <- left_join(ifn_fgsea, ifn_sizes, by=c('celltype', 'condition', 'interaction', 'pathway'='IFN'))
ifn_fgsea <- ifn_fgsea %>% mutate(prop_size=(leadingEdge_size/size))
ggplot(ifn_fgsea) + geom_col(aes(x=celltype, y=prop_size, fill=interaction, alpha=is_sig), position='dodge') + 
  facet_grid(rows=vars(pathway), cols=vars(condition)) + theme_bw()
ggsave('LeadingEdgeSizes_proportion_IFNpathways_barplot.pdf', height=4, width=7)

# add size to ful_df as well
full_df$IFN <- gsub('IFNa', 'INTERFERON_ALPHA', full_df$IFN)
full_df$IFN <- gsub('IFNy', 'INTERFERON_GAMMA', full_df$IFN)
full_df <- left_join(full_df, ifn_fgsea, by=c('IFN'='pathway', 'condition', 'celltype', 'interaction'))

# correlatin bwwn leadingEdge_size and scores pvalues
correlations_ifn <- full_df %>% filter(is_leading_edge==TRUE) %>%
  group_by(IFN, interaction) %>%
  summarise(
    cor_size=cor.test(leadingEdge_size, p.value, method='spearman')$estimate,
    pval_size=cor.test(leadingEdge_size, p.value, method='spearman')$p.value,
    cor_prop=cor.test(prop_size, p.value,  method='spearman')$estimate,
    pval_prop=cor.test(prop_size, p.value,  method='spearman')$p.value,
    .groups='drop')


de_results <- fread('NI_IVAxRV_integrated_limma_results.txt') %>% filter(interaction!='none')

for (i in 1:nrow(ifn_fgsea)){
  int <- ifn_fgsea$interaction[i]
  cond <- ifn_fgsea$condition[i]
  ctype <- ifn_fgsea$celltype[i]
  path <- ifn_fgsea$pathway[i]
  leading_edge <- lapply(ifn_fgsea$leadingEdge[i], function(x) {gsub("\"", "", unlist(strsplit(x, "\"\\|\"")))}) %>% unlist()
  
  sub_logfc_edge <- de_results %>% filter(condition==cond, interaction==int, celltype==ctype, Gene %in% leading_edge) %>%
    select(Gene, logFC, qvals, condition, interaction, celltype) %>% mutate(pathway=path)
  
  if (path=='INTERFERON_ALPHA'){
  sub_logfc <- de_results %>% filter(condition==cond, interaction==int, celltype==ctype, Gene %in% ifn_genes[[1]]) %>%
    select(Gene, logFC, qvals, condition, interaction, celltype) %>% mutate(pathway=path)
  } else {
    sub_logfc <- de_results %>% filter(condition==cond, interaction==int, celltype==ctype, Gene %in% ifn_genes[[2]]) %>%
      select(Gene, logFC, qvals, condition, interaction, celltype) %>% mutate(pathway=path)
  }
  
  if (exists('leadingEdge_logFC')){
    leadingEdge_logFC <- rbind(leadingEdge_logFC, sub_logfc_edge)
  } else {leadingEdge_logFC <- sub_logfc_edge}
  
  if (exists('all_logFC')){
    all_logFC <- rbind(all_logFC, sub_logfc)
  } else {all_logFC <- sub_logfc}
}

(ggplot(all_logFC, aes(x=celltype, y=logFC, fill=interaction)) + 
    facet_grid(rows=vars(pathway), cols=vars(condition)) + theme_bw() +
    geom_violin(position=position_dodge(width=0.8)) + 
    geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8), width=0.1, color='white') + theme_bw() ) /
(ggplot(leadingEdge_logFC, aes(x=celltype, y=logFC, fill=interaction)) + 
  facet_grid(rows=vars(pathway), cols=vars(condition)) + theme_bw() +
  geom_violin(position=position_dodge(width=0.8)) + 
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.8), width=0.1, color='white') + theme_bw()) 
ggsave('IFNpathways_interactions_logFCs_violinplot.pdf', height=8, width=10)

