library(tidyverse)
library(data.table)
library(janitor)
library(ggpointdensity)
library(viridis)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')
conditions <- c('NI', 'IVA')
celltypes <- c('B', 'CD4-T', 'CD8-T', 'Mono', 'NK')

haley_b <- fread('haley_betas.txt') %>% separate(gene_SNP, 
                                                 into=c('gene', 
                                                        'Chromosome', 
                                                        'Position', 'Ref', 'Alt'), 
                                                 sep='_') %>% 
  mutate(snps=Chromosome%&%':'%&%Position) %>% select(-c(Chromosome, Position, Ref, Alt)) %>%
  select(gene, snps, everything()) %>% drop_na()

haley_cols <- colnames(haley_b)
haley_cols <- gsub('monocytes', 'Mono', haley_cols)
haley_cols <- gsub('flu', 'IVA', haley_cols)
haley_cols <- gsub('CD8T', 'CD8-T', haley_cols)
haley_cols <- gsub('CD4T', 'CD4-T', haley_cols)
colnames(haley_b) <- haley_cols

for (c in conditions){
  print(c)
  for (ct in celltypes){
    print(ct)
    daniel_b <- fread('matrixEQTL_results/'%&%c%&%'_'%&%ct%&%'_elbowPCs_cisQTL_sumstats.txt') %>% 
      select(gene, snps, beta)
    
    haley_subset <- haley_b %>% select(gene, snps, ct%&%'_'%&%c)
    haley_daniel_b <- inner_join(haley_subset, daniel_b, by=c('gene', 'snps')) %>%
      mutate(condition=c, celltype=ct)
    colnames(haley_daniel_b)[3:4] <- c('H_beta', 'D_beta')
    
    if (exists('compiled.betas')){
      compiled.betas <- rbind(compiled.betas, haley_daniel_b)
    } else {compiled.betas <- haley_daniel_b}
  }
}


# scatter plot with point density 
ggplot(compiled.betas) + geom_pointdensity(aes(D_beta, H_beta), show.legend=F) +
  stat_smooth(aes(D_beta, H_beta), method='lm', geom='smooth', formula=y~x) +
  geom_abline(slope=1, color='red') + theme_bw() + 
  facet_grid(cols=vars(celltype), rows=vars(condition)) + scale_color_viridis()
ggsave('eQTLmapping_DanielvsHaley_betas.pdf', height=5, width=10)

