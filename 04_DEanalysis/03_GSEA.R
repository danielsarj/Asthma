library(tidyverse)
library(data.table)
library(fgsea)
library(msigdbr)
library(limma)
'%&%' <- function(a,b) paste(a,b, sep = '')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')
conditions <- c('RV', 'IVA')
cells_seurat <- c('B','CD4-T','CD8-T','Mono','NK')
interactions <- c('none','asthma','income')

# retrieve human genes for the hallmark collection gene sets
human.path.list <- msigdbr(species='Homo sapiens', collection='H') %>% 
  split(x=.$gene_symbol, f=.$gs_name)

# read results
DE_results <- fread('NI_IVAxRV_integrated_limma_results.txt')

for (int in interactions){
  for (i in 1:length(conditions)){
    for (ctype in cells_seurat){
      
      # GSEA  
      # subset based on interaction, celltype, and condition
      subset_DE_results <- DE_results %>% filter(celltype==ctype, condition==conditions[i], interaction==int) 
      
      # format values for fGSEA
      subset_DE_results <- subset_DE_results %>% arrange(desc(t)) %>% select(Gene, t)
      subset_DE_results <- setNames(subset_DE_results$t, subset_DE_results$Gene)

      # run fGSEA
      fgseaRes <- fgseaMultilevel(pathways=human.path.list,
                        stats=subset_DE_results, 
                        minSize=15, 
                        maxSize=1000,
                        nPermSimple=100000)
      fgseaRes <- fgseaRes %>% mutate(celltype=ctype, condition=conditions[i], interaction=int) %>% 
        select('pathway','pval','padj','log2err','ES','NES','size','celltype','condition','interaction','leadingEdge')
      
      topPathwaysUp <- fgseaRes[NES>0][head(order(pval), n=5), pathway]
      topPathwaysDown <- fgseaRes[NES<0][head(order(pval), n=5), pathway]
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
      plotGseaTable(human.path.list[topPathways], subset_DE_results, fgseaRes, gseaParam=0.5)
      ggsave('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_'%&%int%&%'_desc_topSigPathways.pdf', height=6, width=10)
      
      genes_in_pathway_a <- human.path.list[['HALLMARK_INTERFERON_ALPHA_RESPONSE']]
      stats_ranked_a <- sort(subset_DE_results, decreasing=TRUE)
      genes_in_pathway_y <- human.path.list[['HALLMARK_INTERFERON_GAMMA_RESPONSE']]
      stats_ranked_y <- sort(subset_DE_results, decreasing=TRUE)     
      plotEnrichment(genes_in_pathway_a, stats_ranked_a) + 
        ggtitle('Enrichment Plot: HALLMARK_INTERFERON_ALPHA_RESPONSE') + 
        plotEnrichment(genes_in_pathway_y, stats_ranked_y) + 
        ggtitle('Enrichment Plot: HALLMARK_INTERFERON_GAMMA_RESPONSE')
      ggsave('NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_'%&%int%&%'_INF_enrichmentplots.pdf', height=6, width=15)
      
      # CAMERA
      cam <- cameraPR(statistic=subset_DE_results, index=human.path.list) %>% rownames_to_column('pathway') %>%
        mutate(celltype=ctype, condition=conditions[i], interaction=int) 
      
      # compile results
      if (exists('compiled.fgseaRes')){
        compiled.fgseaRes <- rbind(compiled.fgseaRes, fgseaRes)
      } else {
        compiled.fgseaRes <- fgseaRes
      }
      
      # compile results
      if (exists('compiled.cameraRes')){
        compiled.cameraRes <- rbind(compiled.cameraRes, cam)
      } else {
        compiled.cameraRes <- cam
      }
    }
  }
}
compiled.fgseaRes$pathway <- gsub('HALLMARK_', '', compiled.fgseaRes$pathway)
fwrite(compiled.fgseaRes, 'NI_IVAxRV_integrated_descGSEAresults.txt', sep=' ', na='NA', col.names=T)
compiled.cameraRes$pathway <- gsub('HALLMARK_', '', compiled.cameraRes$pathway)
fwrite(compiled.cameraRes, 'NI_IVAxRV_integrated_descCAMERAresults.txt', sep=' ', na='NA', col.names=T)

# plot number of enriched pathways per celltype/condition/interaction
summ <- compiled.fgseaRes %>% group_by(celltype, condition, interaction) %>% filter(padj<0.05) %>%
  summarise(n=n())
summ$interaction <- factor(summ$interaction, levels=c('none','asthma','income'))
ggplot(summ) + geom_col(aes(x=celltype, y=n)) + theme_bw() + ylab('Number of enriched pathways') +
  facet_grid(cols=vars(interaction), rows=vars(condition))
ggsave('NI_IVAxRV_descSigGeneSets_acrossInteractions.pdf', height=4, width=8)

# filter CAMERA results by significance
compiled.cameraRes <- compiled.cameraRes %>% filter(FDR<0.05)

# filter GSEA results by CAMERA results
joint_gsea_camera <- inner_join(compiled.cameraRes, compiled.fgseaRes, by=c('pathway', 'condition', 'interaction', 'celltype'))
summ <- joint_gsea_camera %>% group_by(celltype, condition, interaction) %>% 
  filter(padj<0.05, (Direction=='Up' & NES>0) | (Direction=='Down' & NES<0)) %>% summarise(n=n())
summ$interaction <- factor(summ$interaction, levels=c('none','asthma','income'))
ggplot(summ) + geom_col(aes(x=celltype, y=n)) + theme_bw() + ylab('Number of enriched pathways') +
  facet_grid(cols=vars(interaction), rows=vars(condition))
ggsave('NI_IVAxRV_descSigGeneSets_acrossInteractions_afterCAMERAfiltering.pdf', height=4, width=8)
ggsave('NI_IVAxRV_descSigGeneSets_acrossInteractions_afterCAMERAfiltering.png', height=4, width=8)

# bubble plots for significantly enriched pathways 
for (int in interactions){
  joint_gsea_camera %>% filter(padj<0.05, interaction==int) %>% 
    ggplot(., aes(x=celltype, y=pathway, size=-log10(padj), color=NES)) +
    geom_point(alpha=0.8) + scale_size(range=c(1,10)) + scale_color_gradient(low='blue', high='red') +  
    labs(x='Cell type', y=NULL, size='-log10(padj)', color='NES') +
    theme_bw() + facet_wrap(~condition)

  ggsave('NI_IVAxRV_'%&%int%&%'_descSigGeneSets_bubbleplot.pdf', height=5, width=10)
  ggsave('NI_IVAxRV_'%&%int%&%'_descSigGeneSets_bubbleplot.png', height=5, width=10)
}
