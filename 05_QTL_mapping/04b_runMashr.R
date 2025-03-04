library(tidyverse)
library(reshape2) 
library(data.table)
library(mashr)
library(pheatmap)
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

# read in random input dataframe and split into betas/SEs
random_df <- fread('mashr_in_random_df.txt') %>%
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))
random_beta_df <- random_df %>% select(contains('beta')) %>% as.matrix()
random_se_df <- random_df %>% select(contains('SE')) %>% as.matrix()
mash_random <- mash_set_data(random_beta_df, random_se_df)

# estimate correlation structure in random data
Vhat <- estimate_null_correlation_simple(mash_random)

# set up main random object again
mash_random <- mash_set_data(random_beta_df, random_se_df, V=Vhat)
rm(random_df, random_beta_df, random_se_df)

# read in strong input dataframe and split into betas/SEs
strong_df <- fread('mashr_in_strong_df.txt')  %>%
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))
strong_beta_df <- strong_df %>% select(contains('beta')) %>% as.matrix()
strong_se_df <- strong_df %>% select(contains('SE')) %>% as.matrix()
mash_strong <- mash_set_data(strong_beta_df, strong_se_df, V=Vhat)
strong_pairs <- strong_df %>% select(gene, snps)
rm(strong_df, strong_beta_df, strong_se_df)

# obtain data driven covariances
U.pca <- cov_pca(mash_strong, 15)
U.ed <- cov_ed(mash_strong, U.pca)
U.c <- cov_canonical(mash_random)

# create infection covariance matrices
U.NI <- matrix(nrow=15, ncol=15)
U.NI[,1:5] <- c(rep(1,5),rep(0,10))
U.NI[,6:15] <- c(rep(0,15))

U.RV <- matrix(nrow=15, ncol=15)
U.RV[,1:5] <- c(rep(0,15))
U.RV[,11:15] <- c(rep(0,15))
U.RV[,6:10] <- c(rep(0,5),rep(1,5),rep(0,5))

U.IVA <- matrix(nrow=15, ncol=15)
U.IVA[,1:10] <- c(rep(0,15))
U.IVA[,11:15] <- c(rep(0,10),rep(1,5))

U.infection <- list(U.NI, U.RV, U.IVA)
names(U.infection) <- c('NI','RV','IVA')

# fit mash models 
m <- mash(mash_random, Ulist=c(U.ed,U.c,U.infection), outputlevel=1)
m2 <- mash(mash_strong, g=get_fitted_g(m), fixg=TRUE)

# assess sharing of significant signals among each pair of conditions by posterior means
## lfsr = 0.1
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh=0.1, factor=0.5)
colnames(m.pairwise_PM) <- gsub('_beta', '', colnames(m.pairwise_PM))
rownames(m.pairwise_PM) <- gsub('_beta', '', rownames(m.pairwise_PM))
fwrite(m.pairwise_PM, 'mashr_correlation_sigresults_lfsr0.1.txt', quote=F, sep=' ', 
       row.names=T, col.names=T)
pdf('mashr_correlation_sigresults_lfsr0.1.pdf', height=5, width=6)
pheatmap(m.pairwise_PM, scale='none', clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean', clustering_method='complete', 
         angle_col=45)
dev.off()

## lfsr = 0.05
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh=0.05, factor=0.5)
colnames(m.pairwise_PM) <- gsub('_beta', '', colnames(m.pairwise_PM))
rownames(m.pairwise_PM) <- gsub('_beta', '', rownames(m.pairwise_PM))
fwrite(m.pairwise_PM, 'mashr_correlation_sigresults_lfsr0.05.txt', quote=F, sep=' ', 
       row.names=T, col.names=T)
pdf('mashr_correlation_sigresults_lfsr0.05.pdf', height=5, width=6)
pheatmap(m.pairwise_PM, scale='none', clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean', clustering_method='complete', 
         angle_col=45)
dev.off()

# get posterior summaries
p_lfsr <- get_lfsr(m2) # local false sign rate
p_lfsr <- cbind(strong_pairs, p_lfsr)
colnames(p_lfsr) <- gsub('_beta', '_lfsr', colnames(p_lfsr))
p_mean <- get_pm(m2) # new betas
p_mean <- cbind(strong_pairs, p_mean)
p_sd <- get_psd(m2) # standard deviation
p_sd <- cbind(strong_pairs, p_sd)
colnames(p_sd) <- gsub('_beta', '_SD', colnames(p_sd))

# save outputs
fwrite(p_lfsr, 'mashr_out_lfsr_df.txt', quote=F, sep=' ')
fwrite(p_mean, 'mashr_out_beta_df.txt', quote=F, sep=' ')
fwrite(p_sd, 'mashr_out_sd_df.txt', quote=F, sep=' ')