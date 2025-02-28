library(tidyverse)
library(data.table)
library(mashr)
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
U.pca <- cov_pca(mash_strong, ncol(mash_strong[['Bhat']]))
U.ed <- cov_ed(mash_strong, U.pca)
U.c <- cov_canonical(mash_random)

# fit mash models 
m <- mash(mash_random, Ulist=c(U.ed,U.c), outputlevel=1)
m2 <- mash(mash_strong, g=get_fitted_g(m), fixg=TRUE)

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
