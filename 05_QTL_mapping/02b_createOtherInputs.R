library(tidyverse)
library(data.table)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

# create SNP location file
snp_local <- fread('../genotypes/all_genotypes.bim') %>% select(V2, V1, V4)
colnames(snp_local) <- c('snpid', 'chr', 'pos')
fwrite(snp_local, '../genotypes/snp_location.txt', col.names=T, sep='\t')

# create gene location file
gene_local <- fread('../DEanalysis/ensembl_genes.txt') %>% 
  filter(gene_biotype=='protein_coding' & chromosome_name!='MT') %>%
  select(hgnc_symbol, chromosome_name, start_position, end_position)
gene_local$chromosome_name <- as.numeric(gene_local$chromosome_name)
gene_local$hgnc_symbol <- gene_local$hgnc_symbol %>% na_if('')
gene_local <- gene_local %>% drop_na()
colnames(gene_local) <- c('geneid', 'chr', 'left', 'right')
fwrite(gene_local, 'gene_location.txt', col.names=T, sep='\t')

# create dosage file
dosage <- fread('../genotypes/filtered_genotypes_dosage.raw') %>% 
  select(-c(FID, PAT, MAT, SEX, PHENOTYPE)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% row_to_names(row_number=1)
tmp_colnames <- colnames(dosage)
tmp_colnames <- gsub('SEA3', 'SEA-3', tmp_colnames)
colnames(dosage) <- tmp_colnames
dosage$IID <- gsub('_A', '', dosage$IID)
dosage$IID <- gsub('_T', '', dosage$IID)
dosage$IID <- gsub('_C', '', dosage$IID)
dosage$IID <- gsub('_G', '', dosage$IID)
colnames(dosage)[1] <- c('snpid')
fwrite(dosage, '../genotypes/snp_dosage.txt', col.names=T, sep='\t')
