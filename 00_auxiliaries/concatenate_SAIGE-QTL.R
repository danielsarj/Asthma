library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/Saige/step3/outputs/')
celltypes <- c('NK_NI')
pcs <- c('PCs1', 'PCs2', 'PCs3', 'PCs4', 'PCs5','PCs6','PCs7','PCs8','PCs9','PCs10', 
         'PCs11', 'PCs12', 'PCs13', 'PCs14','PCs15')
thresholds <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05)

final_df <- fread('compiled_step3_outputs.txt')
for (ct in celltypes){
  print(ct)
  for (pc in pcs){
    print(pc)
    # find all SAIGE-QTL step3 outputs
    saige_files <- list.files(ct%&%'/exp_'%&%pc, full.names=T)
  
    for (f in saige_files){
      # read file
      tmp <- fread(f) %>% mutate(celltype=ct, exp_pcs=pc)
    
      # combine into a sigle dataframe
      if (exists('final_df')){
      final_df <- rbind(final_df, tmp)
      } else {final_df <- tmp}
    }
  }
}
fwrite(final_df, 'compiled_step3_outputs.txt', sep=' ', col.names=T)
final_df$exp_pcs <- factor(final_df$exp_pcs, levels=pcs)

# read Haleys results
for (k in seq(1:15)){
  print(k)
  tmp <- fread('../../../HALEYs/matrixEQTL_results/NI_NK_adj_'%&%k%&%'PCs_cisQTL_sumstats.txt') %>% group_by(gene) %>% 
    slice_min(pvalue, with_ties=F) %>% ungroup() %>% mutate(exp_pcs='PCs'%&%k)
  
  # combine into a sigle dataframe
  if (exists('haleys')){
    haleys <- rbind(haleys, tmp)
  } else {haleys <- tmp}
}
fwrite(haleys, 'compiled_matrixeqtl_outputs.txt', sep=' ', col.names=T)
haleys_5e5 <- fread('compiled_matrixeqtl_outputs.txt')
haleys_5e5$exp_pcs <- factor(haleys_5e5$exp_pcs, levels=pcs)

for (k in seq(1:15)){
  print(k)
  tmp <- fread('../../../HALEYs/matrixEQTL_results/NI_NK_adj_'%&%k%&%'PCs_1e5_cisQTL_sumstats.txt') %>% group_by(gene) %>% 
    slice_min(pvalue, with_ties=F) %>% ungroup() %>% mutate(exp_pcs='PCs'%&%k)
  
  # combine into a sigle dataframe
  if (exists('haleys_1e5')){
    haleys_1e5 <- rbind(haleys_1e5, tmp)
  } else {haleys_1e5 <- tmp}
}
fwrite(haleys_1e5, 'compiled_matrixeqtl_1e5_outputs.txt', sep=' ', col.names=T)
haleys_1e5 <- fread('compiled_matrixeqtl_1e5_outputs.txt')
haleys_1e5$exp_pcs <- factor(haleys_1e5$exp_pcs, levels=pcs)

# find eGenes
saigeqtl_summary <- final_df %>%
  crossing(threshold = thresholds) %>%
  group_by(celltype, exp_pcs, threshold) %>%
  summarise(
    n_eGenes = sum(top_pval < threshold, na.rm = TRUE),
    all_genes = n(),
    prop = n_eGenes / all_genes,
    .groups = "drop"
  ) %>%
  mutate(method = "saigeqtl")
haleys_1e5_summary <- haleys_1e5 %>%
  crossing(threshold = thresholds) %>%
  group_by(celltype, exp_pcs, threshold) %>%
  summarise(
    n_eGenes = sum(pvalue < threshold, na.rm = TRUE),
    all_genes = n(),
    prop = n_eGenes / all_genes,
    .groups = "drop"
  ) %>%
  mutate(method = "matrixeqtl_1e5")

haleys_5e5_summary <- haleys_5e5 %>%
  crossing(threshold = thresholds) %>%
  group_by(celltype, exp_pcs, threshold) %>%
  summarise(
    n_eGenes = sum(pvalue < threshold, na.rm = TRUE),
    all_genes = n(),
    prop = n_eGenes / all_genes,
    .groups = "drop"
  ) %>%
  mutate(method = "matrixeqtl_5e5")

# merge dfs
summary_all <- rbind(saigeqtl_summary, haleys_1e5_summary, haleys_5e5_summary)
summary_all$celltype <- gsub('_NI', '', summary_all$celltype)

# plot 
ggplot(summary_all, aes(x=exp_pcs, y=n_eGenes, color=method, group=method)) + geom_point() + geom_line() + theme_bw() +
   facet_wrap(~threshold, scales='free')

# what if we subset for the same genes?
dfs <- list(haleys_1e5 = haleys_1e5, haleys_5e5 = haleys_5e5, saigeqtl= final_df)

# unique genes per exp_pcs for each dataframe
genes_by_pcs <- function(df) {
  df %>%
    distinct(exp_pcs, gene) %>%
    group_by(exp_pcs) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
}
d1 <- genes_by_pcs(haleys_1e5)
d2 <- genes_by_pcs(haleys_5e5)
d3 <- genes_by_pcs(final_df)

# find intersection
intersections <- reduce(list(d1, d2, d3), full_join, by = "exp_pcs") %>%
  rowwise() %>%
  mutate(intersect_genes = list(Reduce(intersect, list(genes.x, genes.y, genes)))) %>%
  select(exp_pcs, intersect_genes) %>%
  ungroup()

summarise_on_intersect <- function(df, p_col, method, intersections, thresholds) {
  df %>%
    # keep only intersecting genes
    inner_join(intersections %>% unnest(intersect_genes),
               by = c("exp_pcs", "gene" = "intersect_genes")) %>%
    # expand thresholds
    crossing(threshold = thresholds) %>%
    group_by(celltype, exp_pcs, threshold) %>%
    summarise(
      n_eGenes  = sum({{p_col}} < threshold, na.rm = TRUE),
      all_genes = n(),
      prop      = n_eGenes / all_genes,
      .groups   = "drop"
    ) %>%
    mutate(method = method)
}

haleys_1e5_summary_intersect <- summarise_on_intersect(
  haleys_1e5, pvalue, "matrixeqtl_1e5", intersections, thresholds
)

haleys_5e5_summary_intersect <- summarise_on_intersect(
  haleys_5e5, pvalue, "matrixeqtl_5e5", intersections, thresholds
)

saigeqtl_summary_intersect <- summarise_on_intersect(
  final_df, top_pval, "saigeqtl", intersections, thresholds
)

# combine all
all_summaries <- bind_rows(
  haleys_1e5_summary_intersect,
  haleys_5e5_summary_intersect,
  saigeqtl_summary_intersect
)

all_summaries$celltype <- gsub('_NI', '', all_summaries$celltype)

# plot 
ggplot(all_summaries, aes(x=exp_pcs, y=n_eGenes, color=method, group=method)) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(~threshold, scales='free')

# merge haleys and saige
haleys_5e5$snps <- gsub('_A', '', haleys_5e5$snps)
haleys_5e5$snps <- gsub('_T', '', haleys_5e5$snps)
haleys_5e5$snps <- gsub('_C', '', haleys_5e5$snps)
haleys_5e5$snps <- gsub('_G', '', haleys_5e5$snps)

merged <- inner_join(haleys_5e5, final_df, by=c('gene'='gene', 'snps'='top_MarkerID', 'exp_pcs'='exp_pcs'))
ggplot(merged, aes(x=-log10(pvalue), y=-log10(top_pval))) + geom_point() + theme_bw() +
  facet_wrap(~exp_pcs) + geom_abline(slope=1, color='red')
