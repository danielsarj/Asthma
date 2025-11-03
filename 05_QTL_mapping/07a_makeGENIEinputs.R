library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/GENIE')
input_prefix <- c('IVA_B_7', 'NI_B_7', 'RV_B_16', 'IVA_Mono_1', 'NI_Mono_2', 'RV_Mono_1', 
                  'IVA_NK_2', 'NI_NK_3', 'RV_NK_4', 'IVA_T-CD4_3', 'NI_T-CD4_2', 'RV_T-CD4_10',
                  'IVA_T-CD8_1', 'NI_T-CD8_6', 'RV_T-CD8_7')

# load PLINK .fam file
plink_fam <- fread('filtered_merged.fam') %>% select(V1, V2) %>% rename(FID=V1, IID=V2) 

# load metadata file to create covariate and environment files
metadata <- fread('../../sample_metadata.txt') %>% select(ID, asthma, albuterol, income) %>% 
  mutate(FID=ID) %>% rename(IID=ID) %>% select(FID, IID, asthma, albuterol, income) 
metadata <- plink_fam %>% left_join(metadata, by=c('FID','IID'))
metadata$albuterol <- factor(metadata$albuterol, levels=c('No', 'Yes')) %>% as.numeric() %>% scale()
metadata$asthma <- factor(metadata$asthma, levels=c('No', 'Yes'))  %>% as.numeric() %>% scale()
metadata$income <- ifelse(metadata$income %in% c('< $10,000', '$10,000-$29,999', '$30,000-$49,999'), 'Low', 'High')
metadata$income <- factor(metadata$income, levels=c('Low','High'))  %>% as.numeric() %>% scale()
metadata_asthma <- metadata %>% select(FID, IID, asthma)
metadata_income <- metadata %>% select(FID, IID, income)
metadata_albuterol <- metadata %>% select(FID, IID, albuterol)

# read genetic PCs
gen_pcs <- fread('../PCAIR.eigenvec') %>% select(sample_id, V1, V2, V3, V4)
gen_pcs$sample_id <- gsub('SEA3', 'SEA-3', gen_pcs$sample_id)
gen_pcs <- plink_fam %>% left_join(gen_pcs, , by=c('IID'='sample_id'))
metadata_albuterol <- inner_join(metadata_albuterol, gen_pcs) %>% select(FID, IID, albuterol, V1, V2, V3, V4)
fwrite(metadata_asthma, 'inputs/ENV_asthma.txt', col.names=T, sep=' ', na='NA', quote=F)
fwrite(metadata_income, 'inputs/ENV_income.txt', col.names=T, sep=' ', na='NA', quote=F)
fwrite(metadata_albuterol, 'inputs/COV_albuterol.txt', col.names=T, sep=' ', na='NA', quote=F)
fwrite(gen_pcs, 'inputs/COV_PCsonly.txt', col.names=T, sep=' ', na='NA', quote=F)
rm(metadata, metadata_asthma, metadata_income, metadata_albuterol, gen_pcs)

# read gene annotation file
gene_local <- fread('../../DEanalysis/ensembl_genes.txt') %>% 
  filter(gene_biotype=='protein_coding' & chromosome_name!='MT') %>%
  select(hgnc_symbol, chromosome_name, start_position, end_position)
gene_local$chromosome_name <- as.numeric(gene_local$chromosome_name)
gene_local$hgnc_symbol <- gene_local$hgnc_symbol %>% na_if('')
gene_local <- gene_local %>% drop_na()

for (f in input_prefix){
  base <- str_remove(f, "_[^_]+$")
  print(base)
  
  # load gene expression matrix
  exp_matrix <- fread('../'%&%f%&%'PCs.txt')
  
  # subset gene annotation file and exp matrix
  sub_gene_anno <- gene_local %>% filter(hgnc_symbol %in% exp_matrix$GENES)
  fwrite(sub_gene_anno, 'inputs/'%&%base%&%'_gene_annotation.txt', col.names=T, sep='\t', quote=F)
  exp_matrix <- exp_matrix %>% filter(GENES %in% sub_gene_anno$hgnc_symbol)
  
  for (i in seq(1:nrow(exp_matrix))){
    # create one phenotype file per gene
    tmp_exp_matrix <- exp_matrix[i,] %>% t() %>% as.data.frame() %>% rownames_to_column() 
    colnames(tmp_exp_matrix) <- c('IID', tmp_exp_matrix[1,2])
    tmp_exp_matrix <- tmp_exp_matrix[-1,]
    
    # join with fam info
    tmp_exp_matrix <- plink_fam %>% left_join(tmp_exp_matrix, by=c('IID'))
    tmp_exp_matrix <- tmp_exp_matrix %>% select(FID, IID, everything())
    
    # save file
    fwrite(tmp_exp_matrix, 'inputs/'%&%base%&%'_'%&%colnames(tmp_exp_matrix)[3]%&%'.txt', col.names=T, sep=' ', na='NA', quote=F)
  }
}
