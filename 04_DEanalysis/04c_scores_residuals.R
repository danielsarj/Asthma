library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(msigdbr)
library(ggpubr)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# load seurat object
objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')

# keep only protein coding and non-MT genes
annotations <- fread('ensembl_genes.txt')
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# merge metadata
mdata <- objs@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
objs@meta.data <- mdata

# get human hallmark gene sets
ifn_genes <- msigdbr(species='Homo sapiens', collection='H')  %>% 
  split(x=.$gene_symbol, f=.$gs_name)

# filter for IFN modules 
ifn_genes <- ifn_genes[grep('INTERFERON', names(ifn_genes), value=T)]

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
    mdata$income <- factor(mdata$income, levels=c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', 
                                                  '$50,000-$69,999', '$70,000-$89,999')) %>% as.numeric()
    no_NA_income <- mdata %>% filter(!is.na(income)) 
    no_NA_albuterol <- mdata %>% filter(!is.na(albuterol))
    
    for (interaction_term in c('asthma', 'income')){
      print(interaction_term)
      
      if (interaction_term=='asthma'){
          # filter count matrix (only keep protein coding genes)
          count <- tmp@assays$RNA$counts
          count <- count[,colnames(count) %in% rownames(no_NA_albuterol)]
          count <- count[rownames(count) %in% annotations,]
          zero_var_genes <- apply(count, 1, var) == 0
          count <- count[!zero_var_genes, ]
            
          # transform count into dge object
          count <- DGEList(counts=count)
          count <- calcNormFactors(count)
            
          # define design matrix
          design <- model.matrix(~batch+age+gender+n+albuterol, data=no_NA_albuterol)
          
          # voom
          voom <- voom(count, design, plot=F)
            
          # fit linear model 
          fit <- eBayes(lmFit(voom, design))
            
          # extract residuals
          residuals_mat <- residuals(fit, voom) %>% as.data.frame %>% rownames_to_column('Gene')
            
          # subset to IFN genes
          residuals_mat <- residuals_mat %>% filter(Gene %in% union(ifn_genes[[1]], ifn_genes[[2]])) %>% reshape2::melt()
            
          # compute scores
          ifna <- residuals_mat %>% filter(Gene %in% ifn_genes[[1]]) %>% group_by(variable) %>% 
            summarise(IFNa_score=mean(value))
          ifny<- residuals_mat %>% filter(Gene %in% ifn_genes[[2]]) %>% group_by(variable) %>% 
            summarise(IFNy_score=mean(value))
          IFN_scores <- full_join(ifna, ifny)
          rm(ifna, ifny)
            
          # compute paired scors
          IFN_scores <- IFN_scores %>% separate(variable, into=c('ID', 'condition', 'celltype'), sep='_')
            
          # remove unique indvs
          IFN_scores <- IFN_scores %>% group_by(ID) %>% mutate(is_unique = n() == 1)
          IFN_scores <- IFN_scores %>% 
            filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
            
          for (id in unique(IFN_scores$ID)){
            # select indv
            subset_score <- IFN_scores %>% filter(ID==id) 
              
            # confirm position of NI vs infection
            NI_row <- which(subset_score$condition=='NI')
            Inf_row <- which(subset_score$condition==conditions[i])
              
            # compute inf score
            subset_score$IFNa_score <- subset_score$IFNa_score[Inf_row] - subset_score$IFNa_score[NI_row]
            subset_score$IFNy_score <- subset_score$IFNy_score[Inf_row] - subset_score$IFNy_score[NI_row]
              
            # edit df
            subset_score <- subset_score %>% filter(condition==conditions[i]) %>% 
              select(ID, condition, celltype, IFNa_score, IFNy_score) 
              
            if (exists('paired.bulk.ifn.scores.alb.asthma')){
              paired.bulk.ifn.scores.alb.asthma <- rbind(paired.bulk.ifn.scores.alb.asthma, subset_score)
            } else {paired.bulk.ifn.scores.alb.asthma <- subset_score}
          }
        } else {
      # filter count matrix (only keep protein coding genes)
      count <- tmp@assays$RNA$counts
      count <- count[,colnames(count) %in% rownames(no_NA_income)]
      count <- count[rownames(count) %in% annotations,]
      zero_var_genes <- apply(count, 1, var) == 0
      count <- count[!zero_var_genes, ]
        
      # transform count into dge object
      count <- DGEList(counts=count)
      count <- calcNormFactors(count)
        
      # define design matrix
      design <- model.matrix(~batch+age+gender+n, data=no_NA_income)
        
      # voom
      voom <- voom(count, design, plot=F)
        
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
        
      # extract residuals
      residuals_mat <- residuals(fit, voom) %>% as.data.frame %>% rownames_to_column('Gene')
        
      # subset to IFN genes
      residuals_mat <- residuals_mat %>% filter(Gene %in% union(ifn_genes[[1]], ifn_genes[[2]])) %>% reshape2::melt()
        
      # compute scores
      ifna <- residuals_mat %>% filter(Gene %in% ifn_genes[[1]]) %>% group_by(variable) %>% 
        summarise(IFNa_score=mean(value))
      ifny<- residuals_mat %>% filter(Gene %in% ifn_genes[[2]]) %>% group_by(variable) %>% 
        summarise(IFNy_score=mean(value))
      IFN_scores <- full_join(ifna, ifny)
      rm(ifna, ifny)
        
      # compute paired scors
      IFN_scores <- IFN_scores %>% separate(variable, into=c('ID', 'condition', 'celltype'), sep='_')
        
      # remove unique indvs
      IFN_scores <- IFN_scores %>% group_by(ID) %>% mutate(is_unique = n() == 1)
      IFN_scores <- IFN_scores %>% 
        filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
        
      for (id in unique(IFN_scores$ID)){
        # select indv
        subset_score <- IFN_scores %>% filter(ID==id) 
          
        # confirm position of NI vs infection
        NI_row <- which(subset_score$condition=='NI')
        Inf_row <- which(subset_score$condition==conditions[i])
          
        # compute inf score
        subset_score$IFNa_score <- subset_score$IFNa_score[Inf_row] - subset_score$IFNa_score[NI_row]
        subset_score$IFNy_score <- subset_score$IFNy_score[Inf_row] - subset_score$IFNy_score[NI_row]
          
        # edit df
        subset_score <- subset_score %>% filter(condition==conditions[i]) %>% 
          select(ID, condition, celltype, IFNa_score, IFNy_score) 
          
        if (exists('paired.bulk.ifn.scores.income')){
          paired.bulk.ifn.scores.income <- rbind(paired.bulk.ifn.scores.income, subset_score)
        } else {paired.bulk.ifn.scores.income <- subset_score}
      }
        }
    }
  }
}    

# merge w metadata
paired.bulk.ifn.scores.alb.asthma <- paired.bulk.ifn.scores.alb.asthma %>% full_join(sample_m, by=c('ID'='ID'))    
paired.bulk.ifn.scores.alb.asthma <- paired.bulk.ifn.scores.alb.asthma %>% 
  select(celltype, condition, IFNa_score, IFNy_score, asthma, albuterol) %>%
  drop_na()    

# plot
paired.bulk.ifn.scores.alb.asthma %>% ggplot(., aes(x=condition, y=IFNa_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw()
ggsave('scGSEA_IFNascore_pairedID_adjusted.exp_asthma_boxplots.pdf', height=4, width=6)
paired.bulk.ifn.scores.alb.asthma %>% ggplot(., aes(x=condition, y=IFNy_score, fill=asthma)) + 
  geom_boxplot() + facet_wrap(~celltype, scale='free') + theme_bw()
ggsave('scGSEA_IFNyscore_pairedID_adjusted.exp_asthma_boxplots.pdf', height=4, width=6)
