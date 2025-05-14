library(tidyverse)
library(data.table)
library(janitor)
library(ggpointdensity)
library(viridis)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')
conditions <- c('NI', 'IVA')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')

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
haley_cols <- gsub('CD8T', 'T-CD8', haley_cols)
haley_cols <- gsub('CD4T', 'T-CD4', haley_cols)
colnames(haley_b) <- haley_cols
haley_snps <- haley_b$snps 

daniel_b <- fread('mashr/mashr_out_beta_df.txt') %>% separate(snps, into=c('snpss','alt'), sep='_') %>% 
  select(-c(alt, contains('RV'))) %>% rename(snps=snpss) 

for (c in conditions){
  print(c)
  for (ct in celltypes){
    print(ct)

    daniel_subset <- daniel_b %>% select(gene, snps, c%&%'_'%&%ct%&%'_beta')
    haley_subset <- haley_b %>% select(gene, snps, ct%&%'_'%&%c)
    
    haley_daniel_b <- inner_join(haley_subset, daniel_subset, by=c('gene', 'snps')) %>%
      mutate(condition=c, celltype=ct)
    
    colnames(haley_daniel_b)[3:4] <- c('H_beta', 'D_beta')
    
    if (exists('compiled.betas')){
      compiled.betas <- rbind(compiled.betas, haley_daniel_b)
    } else {compiled.betas <- haley_daniel_b}
  }
}


# scatter plot with point density 
ggplot(compiled.betas) + geom_pointdensity(aes(abs(D_beta), abs(H_beta)), show.legend=F) +
  stat_smooth(aes(abs(D_beta), abs(H_beta)), method='lm', geom='smooth', formula=y~x) +
  geom_abline(slope=1, color='red') + theme_bw() + 
  facet_grid(cols=vars(celltype), rows=vars(condition)) + scale_color_viridis()
ggsave('mashr/eQTLmapping_DanielvsHaley_betas_absbetas_postmash.pdf', height=5, width=10)
