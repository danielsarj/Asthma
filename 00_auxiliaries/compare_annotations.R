library(Seurat)
library(SeuratData)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

# load seurat objects
old_obj <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')
old_obj <- old_obj@meta.data %>% select(orig.ident, IDs, condition, celltype, batch, n, avg_mt, prop)
new_obj <- readRDS('NI_IVA_RV.integrated.pseudobulks_new.rds')
new_obj <- new_obj@meta.data %>% select(orig.ident, IDs, condition, celltype, batch, n, avg_mt, prop)

# merge
full_obj <- inner_join(old_obj, new_obj, by=c('orig.ident', 'IDs', 'condition', 'celltype', 'batch'))

full_obj %>% group_by(condition) %>% summarise(old=sum(n.x), new=sum(n.y))

# absolute number of cells
ggplot(full_obj, aes(x=n.x, y=n.y)) + geom_point() + 
  facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(color='blue') + theme_bw() +
  xlab('Old annotation') + ylab('New annotation')

# proportions of cells
ggplot(full_obj, aes(x=prop.x, y=prop.y)) + geom_point() + 
  facet_grid(cols=vars(condition), rows=vars(celltype)) +
  geom_abline(color='blue') + theme_bw() +
  xlab('Old annotation') + ylab('New annotation')

