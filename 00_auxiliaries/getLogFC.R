library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(ggpubr)
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

# read seurat object
objs <- readRDS('../scRNAanalysis/NI_IVA_RV.integrated.pseudobulks.rds')

# merge metadata
mdata <- objs@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
objs@meta.data <- mdata

for (i in 1:length(conditions)){
  print(c(conditions[i]))

  # celltype specific DE
  for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
    print(ctype)
    
    # extract metadata for subsetting
    meta_df <- objs@meta.data
    filtered_meta <- meta_df %>% filter(celltype==ctype, condition==conditions[i] | condition=='NI')
    
    # subset bulk object
    matching_cells <- rownames(filtered_meta)
    tmp <- subset(objs, cells=matching_cells)
    rm(meta_df, filtered_meta, matching_cells)
    
    # filter count matrix (only keep protein coding genes)
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% annotations,]
    zero_var_genes <- apply(count, 1, var) == 0
    count <- count[!zero_var_genes, ]
    
    # transform count into dge object and computting logCPM
    count <- DGEList(counts=count)
    count <- calcNormFactors(count)
    count <- cpm(count, log=TRUE)
    count <- count %>% as.data.frame()
    
    # compute the difference between logCPM for the same indvs. 
    # same as computing logFC
    for (id in unique(mdata$IDs)){
      # select indv
      subset_count <- count %>% select(contains(id%&%'_'))
      
      if (ncol(subset_count)==2){
        # confirm position of NI vs infection
        NI_col <- grep('NI', colnames(subset_count))
        Inf_col <- grep(conditions[i], colnames(subset_count))
      
        # compute logFC for every gene
        subset_count$logFC <- subset_count[,Inf_col] - subset_count[,NI_col]
      
        # add new columns
        subset_count <- subset_count %>% select(logFC) %>% 
          mutate(condition=conditions[i], celltype=ctype, ID=id) %>%
          rownames_to_column('gene')
      
        # compile results
        if (exists('compiled.logFC')){
          compiled.logFC <- rbind(compiled.logFC, subset_count)
        } else {compiled.logFC <- subset_count}
      }
    }
  }
}

# merge with metadata
full_df <- inner_join(compiled.logFC, sample_m, by=c('ID'='ID'))
full_df$income <- factor(full_df$income, levels=c('< $10,000', '$10,000-$29,999', '$30,000-$49,999', 
                                              '$50,000-$69,999', '$70,000-$89,999')) 

ggplot(full_df, aes(x=logFC, fill=asthma)) + geom_density(alpha=0.5, bw=0.5) + 
  facet_grid(cols=vars(celltype), rows=vars(condition)) 

full_df %>% filter(celltype=='B') %>% ggplot(., aes(x=logFC)) + geom_density(alpha=0.5) + 
  facet_grid(cols=vars(condition), rows=vars(asthma)) 

full_df %>% drop_na() %>% ggplot(., aes(x=logFC, fill=income)) + geom_density(alpha=0.5, bw=0.5) + 
  facet_grid(cols=vars(celltype), rows=vars(condition)) 

for (c in c('RV','IVA')){
  for (ct in c('B','T-CD4','T-CD8','Mono','NK')){
    full_subset <- full_df %>% filter(celltype==ct, condition==c)
    yes_asthma <- full_subset %>% filter(asthma=='Yes') %>% select(logFC) %>% 
      arrange(logFC) %>% rename(logFC_wAsthma = logFC)
    no_asthma <- full_subset %>% filter(asthma=='No') %>% select(logFC) %>% 
      arrange(logFC) %>% rename(logFC_woAsthma = logFC)
    
    pdf(c%&%'_'%&%ct%&%'_logFC_asthma_qqplot.pdf', width=4, height=4)
    qqplot(x=no_asthma$logFC_woAsthma, y=yes_asthma$logFC_wAsthma, 
           main=c%&%' '%&%ct%&%' logFC qqplot', xlab='logFC without asthma', 
           ylab='logFC with asthma')
    abline(c(0,1), col='red')
    dev.off()
  }
}









#### ASTHMA
# load asthma results
asthma_results <- fread('NI_IVAxRV_asthma_limma_results_avglogCPM.filtered.txt')  

# for a given gene, look for interactions with asthma
gene_of_interest <- 'DPYSL3'

# boxplots with t-test
full_df %>% filter(gene==gene_of_interest, condition=='RV', celltype=='Mono') %>% 
  ggplot(., aes(x=asthma, y=logFC)) +
  geom_boxplot() + facet_grid(cols=vars(celltype), rows=vars(condition)) + theme_bw() +
  stat_compare_means(comparisons=list(c('No','Yes')), method='t.test') +
  ylab(gene_of_interest %&% ' logFC')



#### INCOME
# load income results
income_results <- fread('NI_IVAxRV_income_limma_results_avglogCPM.filtered.txt')  

# for a given gene, look for interactions with asthma
gene_of_interest <- 'ETS2'

# boxplots with t-test
full_df %>% drop_na() %>% filter(gene==gene_of_interest) %>% ggplot(., aes(x=income, y=logFC)) +
  geom_boxplot() + facet_grid(rows=vars(celltype), cols=vars(condition)) + theme_bw() +
  ylab(gene_of_interest %&% ' logFC')