library(data.table)
library(tidyverse)
library(patchwork)
'%&%' <- function(a,b) paste(a,b, sep='')
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping')

celltypes <- c('NK', 'CD4_T')
# NK has 10686 genes, 6366 cells
# CD4T has 10809 genes, 73323 cells

# SAIGEQTL
for (ct in celltypes){
  print(ct)
  ## step 1
  print('step1')
  step_files <- list.files('Saige/benchmarks/step1/', pattern='^'%&%ct%&%'_NI_')
  for (f in step_files){
    tmp <- fread('Saige/benchmarks/step1/'%&%f) %>% select(s, max_rss) %>% 
      mutate(step='step1', celltype=ct)
    if (exists('step1_compiled')){
      step1_compiled <- bind_rows(step1_compiled, tmp)
    } else {step1_compiled <- tmp}
  }
  ## step 2
  print('step2')
  step_files <- list.files('Saige/benchmarks/step2/', pattern='^'%&%ct%&%'_NI.*no_perm.tsv$')
  for (f in step_files){
    tmp <- fread('Saige/benchmarks/step2/'%&%f) %>% select(s, max_rss) %>% 
      mutate(step='step2', celltype=ct)
    if (exists('step2_compiled')){
      step2_compiled <- bind_rows(step2_compiled, tmp)
    } else {step2_compiled <- tmp}
  }
  ## step 3
  print('step3')
  step_files <- list.files('Saige/benchmarks/step3/', pattern='^'%&%ct%&%'_NI.*no_perm.tsv$')
  for (f in step_files){
    tmp <- fread('Saige/benchmarks/step3/'%&%f) %>% select(s, max_rss) %>% 
      mutate(step='step3', celltype=ct)
    if (exists('step3_compiled')){
      step3_compiled <- bind_rows(step3_compiled, tmp)
    } else {step3_compiled <- tmp}
  }
}

# MATRIXEQTL
for (ct in celltypes){
  print(ct)
  tmp <- fread('HALEYs/benchmark/'%&%ct%&%'.tsv') %>% select(s, max_rss) %>% 
    mutate(step='unique', celltype=ct)
  if (exists('unique_compiled')){
    unique_compiled <- bind_rows(unique_compiled, tmp)
  } else {unique_compiled <- tmp}
}

# compute summaries
compiled_saige <- bind_rows(step1_compiled, step2_compiled, step3_compiled)
compiled_saige_stats <- compiled_saige %>%
  group_by(celltype, step) %>%
  summarise(runtime=mean(s), max_ram=mean(max_rss)) %>%
  ungroup() %>%
  mutate(
    n_perms = if_else(step == 'step1', 1L, 11L),
    
    # total across all genes and perms (serial)
    runtime_allgenes_allperms = runtime * 10686 * n_perms,
    max_ram_allgenes_allperms = max_ram * 10686 * n_perms,
    
    # factoring in 1000 parallel jobs
    wallclock_allgenes_allperms = runtime_allgenes_allperms / 1000,
    peak_ram_parallel = max_ram * 1000,
    
    # human-readable
    wallclock_hours = wallclock_allgenes_allperms / 3600,
    total_cpu_hours = runtime_allgenes_allperms / 3600,
    peak_ram_parallel_GB = peak_ram_parallel / 1024
  ) %>% select(-n_perms)

unique_stats <- unique_compiled %>%
  mutate(
    wallclock_hours = s / 3600,
    peak_ram_GB = max_rss / 1024)

JOINT <- bind_rows(compiled_saige_stats, unique_stats) %>%
  mutate(method = if_else(step == 'unique', 'matrixeqtl', 'saigeqtl')) %>%
  select(celltype, step, method, wallclock_hours)

ggplot(JOINT, aes(x=method, y=wallclock_hours, fill=step)) + geom_col() + 
  facet_wrap(~celltype) + theme_bw()

