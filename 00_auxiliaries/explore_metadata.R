library(tidyverse)
library(data.table)
library(zipcodeR)
library(sf)
library(tigris)
setwd('/project/lbarreiro/USERS/daniel/asthma_project/DEanalysis')

# load sample metadata
sample_m <- fread('sample_metadata.txt')

ggplot(sample_m) + geom_boxplot(aes(x=age, y=gender, fill=gender), show.legend =F) +
  geom_point(aes(x=age, y=gender), position=position_dodge(width=0.75)) + 
  theme_bw()
ggsave('../gender.age.boxplot.pdf', height=3, width=4)

sum <- sample_m %>% group_by(asthma) %>% summarise(n=n())
ggplot(sum) + geom_col(aes(x=asthma, y=n), position='dodge') + 
  scale_y_continuous(breaks=seq(0, max(sum$n), by=1)) + theme_bw()

ggsave('../asthma.status_barplot.pdf', height=5, width=4)

sum <- sample_m %>% group_by(asthma, gender) %>% summarise(n=n())
ggplot(sum) + geom_col(aes(x=asthma, y=n, fill=gender), position='stack') + 
  scale_y_continuous(breaks=seq(0, sum(sum$n), by=1)) + theme_bw()

ggsave('../asthma.status.by_gender_barplot.pdf', height=5, width=4)


zips <- c(64132,64133,64128,66109,64133,64133,64130,64052,64124,64109,64130,66104,
          64128,64128,64128,64128,64130,64130,64109,64134,64063,66027,64110,64134,
          64130,64050,64050,64050,64155,66102,64052,64106,64106,64106,64133,64133,
          64127,64130,64118,64124,64128,64055,64153,64130,64130,64130)
zip_counts <- as.data.frame(table(zips)) %>%
  rename(zip = zips, freq = Freq) %>%
  mutate(zip = as.character(zip))

# get ZIP code polygons
options(tigris_use_cache=TRUE)
kc_zips <- kc_zips <- c(
  '64101','64102','64105','64106','64108','64109','64110','64111','64112','64113',
  '64114','64116','64117','64118','64119','64120','64123','64124','64125','64126',
  '64127','64128','64129','64130','64131','64132','64133','64134','64136','64137',
  '64138','64139','64145','64146','64147','64149','64151','64152','64153','64154',
  '64155','64156','64157','64158','64161','64163','64164','64165','64166','64167',
  '66101','66102','66103','66104','66105','66106','66109','66111','66112','66115',
  '66118','66160','64150')
zip_shapes <- zctas(cb=FALSE, starts_with=kc_zips)

# merge frequencies into shapefile
zip_shapes <- zip_shapes %>%
  left_join(zip_counts, by = c('ZCTA5CE20' = 'zip'))

ggplot(zip_shapes) +
  geom_sf(aes(fill = freq), color='black') +
  scale_fill_viridis_c(
    option = 'viridis', direction = -1,
    na.value = 'gray90',
    breaks = pretty(range(zip_shapes$freq, na.rm=TRUE)),
    labels = pretty(range(zip_shapes$freq, na.rm=TRUE))) +
  theme_minimal() + labs(fill = 'Frequency')

ggsave('../localization.kansas_city_map.pdf', height=8, width=8)
