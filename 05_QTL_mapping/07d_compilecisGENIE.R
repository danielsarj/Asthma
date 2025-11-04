library(tidyverse)
library(stringr)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/QTLmapping/GENIE')
conditions <- c('NI', 'RV', 'IVA')
celltypes <- c('B', 'T-CD4', 'T-CD8', 'Mono', 'NK')
environments <- c('income', 'asthma')

# function to extract GENIE's results (such a badly formatted output. ew)
extract_values <- function(f) {
  lines <- readLines(f, warn = FALSE)
  
  # anchor the patterns so 'Total h2_g' doesn't match 'Total h2_gxe'
  h2g_line   <- lines[grepl('^\\s*Total\\s+h2_g\\s*:', lines)]
  h2gxe_line <- lines[grepl('^\\s*Total\\s+h2_gxe\\s*:', lines)]
  
  # use str_match to capture groups: value and SE (allow optional spaces and optional colon after SE)
  h2g_m <- if (length(h2g_line) >= 1) str_match(h2g_line[1],
                                                'Total\\s+h2_g\\s*:\\s*([-0-9.]+)\\s*SE\\s*:?\\s*([-0-9.]+)') else matrix(NA, nrow=1, ncol=3)
  h2gxe_m <- if (length(h2gxe_line) >= 1) str_match(h2gxe_line[1],
                                                    'Total\\s+h2_gxe\\s*:\\s*([-0-9.]+)\\s*SE\\s*:?\\s*([-0-9.]+)') else matrix(NA, nrow=1, ncol=3)
  
  h2g_val    <- as.numeric(h2g_m[1,2])
  h2g_se     <- as.numeric(h2g_m[1,3])
  h2gxe_val  <- as.numeric(h2gxe_m[1,2])
  h2gxe_se   <- as.numeric(h2gxe_m[1,3])
  
  # parse filename robustly: condition_celltype_gene..._environment.txt
  fname <- basename(f)
  fname_no_ext <- sub('\\.txt$', '', fname)
  parts <- strsplit(fname_no_ext, '_')[[1]]
  environment <- tail(parts, 1)
  condition <- parts[1]
  celltype <- if (length(parts) >= 2) parts[2] else NA
  if (length(parts) > 3) {
    gene <- paste(parts[3:(length(parts)-1)], collapse = '_')
  } else if (length(parts) == 3) {
    gene <- parts[3]
  } else {
    gene <- NA
  }
  
  tibble(
    condition = condition,
    celltype  = celltype,
    gene      = gene,
    h2g       = h2g_val,
    h2g_se    = h2g_se,
    h2gxe     = h2gxe_val,
    h2gxe_se  = h2gxe_se,
    environment = environment
  )
}


for (cond in conditions){
  for (ct in celltypes){
    for (env in environments){
      # list all files matching the pattern
      files <- list.files('outputs/'%&%cond%&%'_'%&%ct, pattern='_'%&%env%&%'\\.txt$', full.names=TRUE)
      
      # run over all files
      results <- files %>% map(~ safely(extract_values)(.x)) %>% map('result') %>% compact() %>% list_rbind()
      
      # compile results
      if (exists('compiled.results')){
        compiled.results <- rbind(compiled.results, results)
      } else {compiled.results <- results}
    }
  }
}
compiled.results <- compiled.results %>% drop_na()
compiled.results$h2g <- as.numeric(compiled.results$h2g)
compiled.results$h2gxe <- as.numeric(compiled.results$h2gxe)
fwrite(compiled.results, 'compiled.results.txt', sep=' ')

sig_results <- compiled.results %>% mutate(sig_h2g=h2g-(2*h2g_se)>0, sig_h2gxe=h2gxe-(2*h2gxe_se)>0) %>% filter(sig_h2g==TRUE & sig_h2gxe==TRUE)

sig_results %>% filter(h2g>0 & h2g <1 & h2gxe>0 & h2gxe<1) %>% 
  ggplot(.) + geom_point(aes(x=h2gxe, y=h2g)) + theme_bw() + facet_grid(cols=vars(celltype), rows=vars(condition)) +
  geom_abline(color='red')