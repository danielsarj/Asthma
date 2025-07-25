library(tidyverse)
library(data.table)
library(fgsea)
library(msigdbr)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','T-CD4','T-CD8','Mono','NK')
interactions <- c('none','asthma', 'asthma_alb','income')

# retrieve human genes for the hallmark collection gene sets
human.path.list <- msigdbr(species='Homo sapiens', collection='H') %>% 
  split(x=.$gene_symbol, f=.$gs_name)

for (int in interactions){
  # load DE limma results
  if (int=='none'){
    DE_results <- fread('NI_IVAxRV_limma_results_mod.txt') %>% rename(gene=Gene)
  } else {
  DE_results <- fread('NI_IVAxRV_'%&%int%&%'_limma_results_mod.txt')
  }
  
  for (i in 1:length(conditions)){
    for (ctype in cells_seurat){
      
      # subset based on celltype and condition
      subset_DE_results <- DE_results %>% filter(celltype==ctype, condition==conditions[i]) 
      
      # format values for fGSEA
      subset_DE_results <- subset_DE_results %>% arrange(t) %>% select(gene, t)
      subset_DE_results <- setNames(subset_DE_results$t, subset_DE_results$gene)
      
      # run fGSEA
      fgseaRes <- fgsea(pathways=human.path.list,
                        stats=subset_DE_results, 
                        minSize=15, 
                        maxSize=1000)
      fgseaRes <- fgseaRes %>% mutate(celltype=ctype, condition=conditions[i], interaction=int) %>% 
        select('pathway','pval','padj','log2err','ES','NES','size','celltype','condition','interaction','leadingEdge')
      
      topPathwaysUp <- fgseaRes[NES>0][head(order(pval), n=5), pathway]
      topPathwaysDown <- fgseaRes[NES<0][head(order(pval), n=5), pathway]
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
      plotGseaTable(human.path.list[topPathways], subset_DE_results, fgseaRes, gseaParam=0.5)
      ggsave('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_'%&%int%&%'_topSigPathways_mod.pdf', height=6, width=10)
      
      # compile results
      if (exists('compiled.fgseaRes')){
        compiled.fgseaRes <- rbind(compiled.fgseaRes, fgseaRes)
      } else {
        compiled.fgseaRes <- fgseaRes
      }
    }
  }
}
fwrite(compiled.fgseaRes, 'NI_IVAxRV_GSEAresults_mod.txt', sep=' ', na='NA', col.names=T)

# plot number of enriched pathways per celltype/condition/interaction
summ <- compiled.fgseaRes %>% group_by(celltype, condition, interaction) %>% filter(padj<0.05) %>%
  summarise(n=n())
ggplot(summ) + geom_col(aes(x=celltype, y=n)) + theme_bw() + 
  facet_grid(cols=vars(interaction), rows=vars(condition))
ggsave('NI_IVAxRV_SigGeneSets_acrossInteractions_mod.pdf', height=4, width=8)

# bubble plots
for (int in interactions){
  compiled.fgseaRes %>% filter(padj<0.05, abs(NES)>=1.5, interaction==int) %>% 
    ggplot(., aes(x=NES, y=reorder(pathway, NES), size=size, color=-log10(padj))) +
    geom_point(alpha=0.8) + scale_size(range=c(1,10)) + scale_color_gradient(low='blue', high='red') +  
    labs(x='Normalized Enrichment Score (NES)', y='Pathway', size='Gene Count', color='-log10(FDR)') +
    theme_bw() + facet_grid(cols=vars(celltype), rows=vars(condition))
  ggsave('NI_IVAxRV_SigGeneSets_'%&%int%&%'_bubbleplot_mod.pdf', height=8, width=12)
}
