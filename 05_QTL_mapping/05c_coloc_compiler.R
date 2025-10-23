library(tidyverse)
library(data.table)
library(ggrepel)
library(viridis)
"%&%" <- function(a,b) paste(a,b, sep="")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/colocalization')

# set file vectors
coloc_outputs <- list.files(pattern='_coloc_results.txt')

for (f in coloc_outputs){
  info <- str_split_1(f, pattern='_')
  
  # read coloc results and get the SNP-Gene pair with highest H4 posterior probabilities
  tmp <- fread(f) %>% group_by(gene) %>% slice_max(PP.H4.abf, with_ties=TRUE) %>%
    slice_max(SNP.PP.H4, with_ties=FALSE) %>% ungroup %>% mutate(condition=info[1], celltype=info[2], gwas=info[3])

  if (exists('coloc.compiled')){
    coloc.compiled <- rbind(coloc.compiled, tmp)
  } else {coloc.compiled <- tmp}
}

# reformat coloc.compiled to make a manhattan plot
coloc_plot <- coloc.compiled %>% separate(snp, into=c('chr', 'pos'), sep=':', remove=FALSE) %>%
  mutate(chr=factor(chr, levels=as.character(1:22)),
         pos=as.numeric(pos)) %>% drop_na(chr, pos) %>% arrange(chr, pos) %>% 
  group_by(chr) %>% mutate(chr_len=max(pos)) %>% ungroup() %>% 
  mutate(chr_start=cumsum(lag(chr_len, default=0)), pos_cum=pos+chr_start)

# compute chromosome center for labeling
axis_df <- coloc_plot %>% group_by(chr) %>% summarize(center=mean(pos_cum))

# find top SNP (max PP.H4.abf) per facet (condition × gwas)
top_hits <- coloc_plot %>% group_by(condition, gwas) %>%
  slice_max(PP.H4.abf, n=1, with_ties=FALSE) %>% ungroup()

# Plot
ggplot(coloc_plot, aes(x=pos_cum, y=PP.H4.abf, color=as.factor(as.numeric(chr) %% 2), shape=celltype)) +
  geom_point() + scale_color_viridis_d() +
  scale_x_continuous(label=axis_df$chr, breaks=axis_df$center) +
  labs(x='Chromosome', y='Posterior Probability of Colocalization (PP.H4)') +
  theme_bw() + theme(legend.position='none') + facet_grid(cols=vars(condition), rows=vars(gwas)) +
  geom_text_repel(data=top_hits, aes(label=paste0(gene,' (',celltype,')')),
    size=2.8, color='black', segment.color='gray60', box.padding=0.4, point.padding=0.3,
    max.overlaps=Inf)