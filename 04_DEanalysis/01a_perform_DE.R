library(Seurat)
library(SeuratData)
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/scRNAanalysis')
conditions <- c('RV', 'IVA')
celltypes <- c('B','CD4-T','CD8-T','Mono','NK')

# load gene annotation from ensembl
annotations <- fread('../DEanalysis/ensembl_genes.txt')

# load sample metadata
sample_m <- fread('../sample_metadata.txt')

# keep only protein coding and non-MT genes
annotations <- annotations$hgnc_symbol[
  annotations$gene_biotype=='protein_coding' &
    annotations$hgnc_symbol!='' &
    !grepl('^MT-', annotations$hgnc_symbol)]

# load seurat object
objs <- readRDS('NI_IVA_RV.integrated.pseudobulks.rds')

# merge metadata
mdata <- objs@meta.data
mdata <- inner_join(mdata, sample_m, by=c('IDs'='ID')) %>% column_to_rownames('orig.ident')
mdata$batch <- factor(mdata$batch, levels=c('B1','B2','B3','B4'))
objs@meta.data <- mdata

ggplot(mdata, aes(x=reorder(IDs, as.numeric(batch)), y=n, fill=batch)) + geom_col() + theme_bw() +
  facet_grid(cols=vars(factor(condition, levels=c('NI','IVA','RV'))), 
             rows=vars(celltype), scales='free') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab(NULL)
ggsave('Pseudobulksizes_byCondition.pdf', height=6, width=12)

ggplot(mdata, aes(x=reorder(IDs, as.numeric(batch)), y=n, fill=asthma)) + geom_col() + theme_bw() +
  facet_grid(cols=vars(factor(condition, levels=c('NI','IVA','RV'))), 
             rows=vars(celltype), scales='free') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab(NULL)
ggsave('../DEanalysis/Pseudobulksizes_byCondition_Asthmastatus.pdf', height=6, width=12)

# define minimum average logCPM thresholds
logCPMfilter_table <- data.frame(celltype=c('B','CD4-T','CD8-T','Mono','NK',
                                            'B','CD4-T','CD8-T','Mono','NK'),
                                 threshold=c(4.9,1.9,1,3.4,5.6,
                                             3.5,3.6,3.1,3.4,5.6),
                                 condition=c(rep('IVA',5),rep('RV',5)))

# condition-celltype specific DE
for (i in 1:length(conditions)){
  print(c(conditions[i]))
  for (ctype in celltypes){
    print(ctype)
    
    # subset pseudobulk object
    tmp <- subset(objs, celltype==ctype & (condition==conditions[i] | condition=='NI'))
    
    # extract metadata and count matrices
    mdata <- tmp@meta.data
    mdata$condition <- factor(mdata$condition, levels=c('NI', conditions[i]))
    mdata$gender <- factor(mdata$gender, levels=c('Male','Female'))

    # remove non protein coding genes from count matrix and genes with variance == 0
    count <- tmp@assays$RNA$counts
    count <- count[rownames(count) %in% annotations,]
    zero_var_genes <- apply(count, 1, var) == 0
    count <- count[!zero_var_genes, ]
    count <- DGEList(counts=count)
    
    # remove lowly expressed genes based on logCPM threshold
    logcpm_threshold <- logCPMfilter_table %>% filter(celltype==ctype, condition==conditions[i]) %>%
      pull(threshold)
    logCPM_pass <- cpm(count, log=TRUE) %>% rowMeans() %>% as.data.frame() %>% filter(.>=logcpm_threshold) %>%
      rownames_to_column() %>% pull(rowname)
    count <- count[logCPM_pass, , keep.lib.sizes=FALSE]
    count <- calcNormFactors(count)
    
    # define design matrix
    design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=mdata)
    
    # voom
    voom <- voom(count, design, plot=F)
    
    # save voom-adjusted expression table and weights
    exp <- voom$E %>% as.data.frame() %>% rownames_to_column('Gene')
    fwrite(exp, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_voom_expression.txt', sep=' ')
    wgts <- voom$weights %>% as.data.frame() %>% rownames_to_column('Gene')
    wgts$Gene <- exp$Gene
    colnames(wgts) <- colnames(exp)
    fwrite(wgts, 'NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_voom_weights.txt', sep=' ')
    rm(exp, wgts)
    
    # fit linear model 
    fit <- eBayes(lmFit(voom, design))
    
    # get results
    og_results <- topTable(fit, coef=ncol(fit), number=Inf) %>% 
      rownames_to_column('Gene') %>% mutate(condition=conditions[i])
    
    # now do permutations where i shuffle condition labels in metadata
    for (j in (1:10)){
      
      # shuffle condition labels preserving individual pairing
      permuted_mdata <- mdata
      for (ind in unique(permuted_mdata$IDs)){
        # flip coin to decide if condition label will be reversed or not
        if (runif(1)<0.5){
          ix <- permuted_mdata$IDs == ind
          permuted_mdata$condition[ix] <- rev(permuted_mdata$condition[ix])
        }
      }
    
      # define design matrix
      design <- model.matrix(~batch+age+gender+n+avg_mt+condition, data=permuted_mdata)
      
      # voom
      voom <- voom(count, design, plot=F)
      
      # fit linear model 
      fit <- eBayes(lmFit(voom, design))
    
      # save pvalues from permutation
      tmp_perm <- topTable(fit, coef=ncol(fit), number=Inf) %>% rownames_to_column('Gene') %>%
        select(Gene, P.Value)
    
      if (exists('compiled_perms')){
        compiled_perms <- inner_join(compiled_perms, tmp_perm, by='Gene')
      } else {compiled_perms <- tmp_perm}
    }
    
    # reorder compiled_perms df so gene order matches OG results
    compiled_perms <- compiled_perms[match(og_results$Gene, compiled_perms$Gene), ]
    compiled_perms <- compiled_perms %>% select(-Gene)
    
    # compute qvalues
    empP <- empPvals(stat=-log10(og_results$P.Value), stat0=-log10(as.matrix(compiled_perms[1:j])), pool=TRUE)
    og_results$qvals <- qvalue(empP)$qvalue
    
    # qqplot 
    pdf('../DEanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_limma_results_qqplot.pdf', width=4, height=4)
    qqplot(x=-log10(compiled_perms[,1]), y=-log10(og_results$P.Value), main=conditions[i]%&%' '%&%ctype, 
           xlab='-log10(permuted p-values)', ylab='-log10(true p-values)')
    abline(c(0,1), col='red')
    dev.off()
    
    # save result
    fwrite(og_results, '../DEanalysis/NI_'%&%conditions[i]%&%'_'%&%ctype%&%'_limma_results_wqvals.txt',
           sep=' ', col.names=T)
    rm(compiled_perms)
  }
}
