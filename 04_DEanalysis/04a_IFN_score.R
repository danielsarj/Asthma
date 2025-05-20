library(Seurat)
library(SeuratData)
library(data.table)
library(msigdbr)
library(janitor)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# retrieve INF gamma and alpha genes for the hallmark collection gene sets
gene.list <- msigdbr(species='Homo sapiens', collection='H') %>% 
  split(x=.$gene_symbol, f=.$gs_name)
gene.list <- gene.list[grep('INTERFERON', names(gene.list), value=T)]

# load seurat object
objs <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')

for (cond in c('IVA','RV')){
  print(cond)
  # celltype specific expression
  for (ctype in c('B','T-CD4','T-CD8','Mono','NK')){
    print(ctype)
    
    # extract metadata
    meta_df <- objs@meta.data
    filtered_meta <- meta_df[meta_df$celltype==ctype & (meta_df$condition=='NI' | meta_df$condition==cond),]
    
    # remove unique IDs
    filtered_meta <- filtered_meta %>% group_by(IDs) %>% mutate(is_unique = n() == 1)
    filtered_meta <- filtered_meta %>% filter(is_unique==FALSE) %>% select(-is_unique) %>% as.data.frame()
    rownames(filtered_meta) <- filtered_meta$orig.ident
    
    # subset bulk object
    matching_cells <- rownames(filtered_meta)
    tmp <- subset(objs, cells=matching_cells)
    
    # extract count 
    count_df <- tmp@assays$RNA$counts %>% as.data.frame()
    
    # IFN alpha
    sub_count_df <- count_df[rownames(count_df) %in% gene.list[[1]],]
    sub_count_df <- sub_count_df[apply(sub_count_df, 1, sd) != 0, ]
    sub_count_df <- sub_count_df %>% t() %>% scale() %>% t() %>% as.data.frame()
    sub_count_df <- colMeans(sub_count_df) %>% as.data.frame()
    colnames(sub_count_df) <- names(gene.list)[1]
    sub_count_df <- sub_count_df %>% rownames_to_column('ID') %>% 
      separate(ID, into=c('ID','condition','celltype'), sep='_')
  
    # IFN gamma
    sub_count_df_2 <- count_df[rownames(count_df) %in% gene.list[[2]],]
    sub_count_df_2 <- sub_count_df_2[apply(sub_count_df_2, 1, sd) != 0, ]
    sub_count_df_2 <- sub_count_df_2 %>% t() %>% scale() %>% t() %>% as.data.frame()
    sub_count_df_2 <- colMeans(sub_count_df_2) %>% as.data.frame()
    colnames(sub_count_df_2) <- names(gene.list)[2]
    sub_count_df_2 <- sub_count_df_2 %>% rownames_to_column('ID') %>% 
      separate(ID, into=c('ID','condition','celltype'), sep='_')
  
    joint_df <- full_join(sub_count_df, sub_count_df_2, by=c('ID','condition','celltype'))
    
    for (ifn in c('ALPHA','GAMMA')){
      print(ifn)
  
      sub_compiled <- joint_df %>% select(ID, condition, celltype, contains(ifn))
      
      for (id in unique(sub_compiled$ID)){
        # select indv
        subset_score <- sub_compiled %>% filter(ID==id) %>% as.data.table() %>% melt(measure.vars=4)
        
        # confirm position of NI vs infection
        NI_row <- which(subset_score$condition=='NI')
        Inf_row <- which(subset_score$condition==cond)
          
        # compute inf score
        subset_score$score <- subset_score$value[Inf_row] - subset_score$value[NI_row]
          
        # add columns
        subset_score <- subset_score %>% filter(condition==cond) %>% select(ID, celltype, score, condition) %>% 
          mutate(interferon=ifn) 
          
        if (exists('ifn.scores')){
          ifn.scores <- rbind(ifn.scores, subset_score)
        } else {ifn.scores <- subset_score}
      }
    }
  }
}

# merge w metadata
ifn.scores <- ifn.scores %>% full_join(sample_m, by=c('ID'))
    
ifn.scores %>% drop_na() %>% ggplot(., aes(x=condition, y=score, fill=asthma)) + geom_boxplot() + 
  facet_grid(cols=vars(interferon), rows=vars(celltype), scale='free') + theme_bw() 
ggsave('IFNscore_asthma.pdf', height=7, width=6)

ifn.scores %>% mutate(across(where(is.character), ~na_if(.x, ''))) %>% drop_na() %>% ggplot(., aes(x=condition, y=score, fill=income)) + geom_boxplot() + 
  facet_grid(cols=vars(interferon), rows=vars(celltype), scale='free') + theme_bw() 
ggsave('IFNscore_income.pdf', height=7, width=6)
