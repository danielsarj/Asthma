library(tidyverse)
library(data.table)
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')

# load sample metadata
sample_m <- fread('sample_metadata.txt')

ggplot(sample_m) + geom_boxplot(aes(x=age, y=gender, fill=gender), show.legend =F) +
  geom_point(aes(x=age, y=gender), position=position_dodge(width=0.75)) + 
  theme_bw()

ggsave('../Gender.Age.boxplot.pdf', height=3, width=4)

sum <- sample_m %>% group_by(asthma) %>% summarise(n=n())
ggplot(sum) + geom_col(aes(x=asthma, y=n), show.legend=F, position='dodge') + 
  scale_y_continuous(breaks=seq(0, max(sum$n), by=1)) + theme_bw()

ggsave('../Asthma.status_barplot.pdf', height=5, width=4)