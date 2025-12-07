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

# save compiled best coloc results
fwrite(coloc.compiled, 'best_coloc_results.txt', sep=' ')
coloc.compiled$celltype <- gsub('T-CD4', 'CD4-T', coloc.compiled$celltype)
coloc.compiled$celltype <- gsub('T-CD8', 'CD8-T', coloc.compiled$celltype)

# reformat coloc.compiled to make a manhattan plot
coloc_plot <- coloc.compiled %>% separate(snp, into=c('chr', 'pos'), sep=':', remove=FALSE) %>%
  mutate(chr=factor(chr, levels=as.character(1:22)),
         pos=as.numeric(pos)) %>% drop_na(chr, pos) %>% arrange(chr, pos) %>% 
  group_by(chr) %>% mutate(chr_len=max(pos)) %>% ungroup() %>% 
  mutate(chr_start=cumsum(lag(chr_len, default=0)), pos_cum=pos+chr_start)
coloc_plot$condition <- factor(coloc_plot$condition, levels=c('NI', 'IVA', 'RV'))

# compute chromosome center for labeling
axis_df <- coloc_plot %>% group_by(chr) %>% summarize(center=mean(pos_cum))

# find top SNP (max PP.H4.abf) per facet (condition × gwas × celltype)
top_hits <- coloc_plot %>% group_by(condition, celltype, gwas) %>%
  slice_max(PP.H4.abf, n=1, with_ties=FALSE) %>% ungroup()

# Plot
ggplot(coloc_plot, aes(x=pos_cum, y=PP.H4.abf, color=as.factor(as.numeric(chr) %% 2), shape=celltype)) +
  geom_point() + scale_x_continuous(label=axis_df$chr, breaks=axis_df$center) +
  labs(x='Chromosome', y='Posterior Probability of Colocalization (PP.H4)') +
  theme_bw() + theme(legend.position='none') + facet_grid(cols=vars(condition), rows=vars(gwas)) +
  geom_text_repel(data=top_hits, aes(label=paste0(gene,' (',celltype,')')),
    size=2.8, color='black', segment.color='gray60', box.padding=0.4, point.padding=0.3,
    max.overlaps=Inf)
ggsave('coloc_volcanoplot_allgenes.pdf', height=5, width=12)
ggsave('coloc_volcanoplot_allgenes.png', height=5, width=12)

# subset to mash-significant eGenes
mash_sig <- fread('../mashr/mashr_out_allstats_df.txt') %>% group_by(condition, celltype) %>%
  filter(lfsr<0.05) %>% ungroup() %>% separate(snps, into=c('snp', 'effectallele'), sep='_')
mash_sig$celltype <- gsub('T-CD4', 'CD4-T', mash_sig$celltype)
mash_sig$celltype <- gsub('T-CD8', 'CD8-T', mash_sig$celltype)
sub_coloc.compiled <- coloc.compiled %>% inner_join(mash_sig, by=c('gene', 'snp', 'condition', 'celltype'))
sub_coloc.compiled$condition <- factor(sub_coloc.compiled$condition, levels=c('NI', 'IVA', 'RV'))

# find top SNP (max PP.H4.abf) per condition × gwas × celltype
sub_top_hits <- sub_coloc.compiled %>% group_by(condition, celltype, gwas) %>%
  slice_max(PP.H4.abf, n=1, with_ties=FALSE) %>% ungroup()

ggplot(sub_coloc.compiled, aes(x=celltype, y=PP.H4.abf, fill=condition)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6, position=position_dodge(width=0.8)) +
  geom_point(aes(color=gwas), position=position_jitterdodge(jitter.width=0.15, dodge.width=0.8),
    alpha=0.7, size=2) + facet_wrap(~ gwas, scales='free_y') + scale_fill_brewer(palette='Set2') +
  labs(y='Posterior probability (PP.H4)', x='Cell type') + theme_bw() +
  ggrepel::geom_text_repel(data=sub_top_hits, aes(label=paste0(gene, ' (', condition, ')')),
    size=2.8, color='black', segment.color='gray60', box.padding=0.4, point.padding=0.3, max.overlaps=Inf)
ggsave('coloc_PPs_boxplot_mashgenes.pdf', height=4, width=10)
ggsave('coloc_PPs_boxplot_mashgenes.png', height=4, width=10)
