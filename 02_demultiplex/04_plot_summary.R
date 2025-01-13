library(tidyverse)
library(data.table)
library(janitor)
library(forcats)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project/lbarreiro/USERS/daniel/asthma_project/alignment')

batches <- c('B1','B2','B3','B4')
conditions <- c('NI', 'RV', 'IVA')

for (b in 1:length(batches)){
  for (c in 1:length(conditions)){
    f <- fread(batches[b] %&%'-' %&% conditions[c] %&%'_GRCh38/demuxalot/demuxalot_summary.tsv') %>% 
      mutate(batch=batches[b], condition=conditions[c]) %>% clean_names()
    if (exists('summary_df')){
      summary_df <- rbind(summary_df, f)
    } else {summary_df <- f}
  }
}
summary_df$classification <- gsub('\\+', '', summary_df$classification) 
summary_df$classification <- gsub('\\_.*', '', summary_df$classification)

summary_df <- summary_df %>% group_by(classification, batch, condition) %>% 
  summarise('n_assigs'=sum(assignment_n))

fwrite(summary_df, 'demuxalot_assignments_summary.txt', sep=' ', col.names=T)

summary_df %>% filter(classification %in% c('doublet')==F) %>%
  ggplot(.) + geom_col(aes(x=fct_reorder(classification, batch), y=n_assigs, fill=batch)) + coord_flip() +
  facet_wrap(~condition) + theme_bw()

ggsave('demultiplexing_summary.pdf', height=6, width=8)

summary_df %>% 
  ggplot(.) + geom_col(aes(x=fct_reorder(classification, batch), y=n_assigs, fill=batch)) + coord_flip() +
  facet_wrap(~condition) + scale_y_continuous(breaks = seq(0,16000, by=2000)) + theme_bw()

ggsave('demultiplexing_summary_wdoublet.pdf', height=6, width=10)