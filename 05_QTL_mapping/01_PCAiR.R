library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/genotypes')

# make GDS
snpgdsBED2GDS(bed.fn='filtered_genotypes.bed',
              bim.fn='filtered_genotypes.bim',
              fam.fn='filtered_genotypes.fam',
              out.gdsfn='filtered_genotypes.gds')
  
# kinship estimation
gds <- snpgdsOpen('filtered_genotypes.gds')
kinMat <- gds %>% snpgdsIBDKING()
samples <- kinMat$sample.id
kinMat <- kinMat$kinship
colnames(kinMat) <- samples
row.names(kinMat) <- samples

# LD prune
snpset <- snpgdsLDpruning(gds, method='corr', slide.max.bp=10e6,
                            ld.threshold=sqrt(0.3), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)

# run PCAiR
mypcair <- pcair(gds, kinobj=kinMat, divobj=kinMat, snp.include=pruned)
pdf('PCAIR_PC1_VS_PC2.pdf', height=4, width=4, onefile=F)
plot(mypcair)
dev.off()

# retrieve eigenvectors and eigenvalues
eigenvec <- mypcair$vectors %>% as.data.frame() %>% rownames_to_column(var='sample_id')
eigenval <- mypcair$values %>% as.data.frame()
fwrite(eigenvec, '../QTLmapping/PCAIR.eigenvec', col.names=T, row.names=F, sep='\t')
fwrite(eigenval, '../QTLmapping/PCAIR.eigenval')

# plot PCs and elbow plot
eigenval <- eigenval %>% mutate(perc=./sum(eigenval[1])*100)
ggplot(eigenval) + geom_line(aes(x=c(seq(1:nrow(eigenval))), y=perc)) + 
  scale_x_continuous(breaks=seq(1, nrow(eigenval), 1)) + theme_bw()
ggsave('../QTLmapping/PCAIR_elbowplot.pdf', height=4, width=7)

ggplot(eigenvec) + geom_point(aes(x=V1, y=V2)) + theme_bw() +
  xlab('PC1 ('%&%round(eigenval$perc[1], digits=2)%&%'%)') + 
  ylab('PC2 ('%&%round(eigenval$perc[2], digits=2)%&%'%)')
ggsave('../QTLmapping/PCAIR_elbowplot.pdf', height=4, width=4)