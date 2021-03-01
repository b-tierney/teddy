#parse lmm output

library(tidyverse)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/lmm_cags/')

gada_db <- read.csv('diabimmune_association_outputs_gada.csv')
#gada_db=gada_db %>% filter(message=='') %>% mutate(by=p.adjust(p.value,method='BY'))
gada_db=gada_db %>% mutate(by=p.adjust(p.value,method='BY'))
gada_teddy <- read.csv('teddy_association_outputs_gada.csv')

iaa_db <- read.csv('diabimmune_association_outputs_iaa.csv')
#iaa_db=iaa_db  %>% filter(message=='') %>% mutate(by=p.adjust(p.value,method='BY'))
iaa_db=iaa_db %>% mutate(by=p.adjust(p.value,method='BY'))
iaa_teddy <- read.csv('teddy_association_outputs_iaa.csv')

ia2a_db <- read.csv('diabimmune_association_outputs_ia2a.csv')
ia2a_db=ia2a_db %>% mutate(by=p.adjust(p.value,method='BY'))
#ia2a_db=ia2a_db  %>% filter(message=='') %>% mutate(by=p.adjust(p.value,method='BY'))
ia2a_teddy <- read.csv('teddy_association_outputs_ia2a.csv')

###at time of sampling


db <- gada_db %>% filter(message=='' & round(by,2)<=.05 & term == 'dataset[, antibody]') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
teddy <- gada_teddy %>% filter(p.value<=.05 & term == 'dataset[, antibody]1') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% db$CAG)

print(teddy)

db <- iaa_db %>% filter(message=='' & round(by,2)<=.05 & term == 'dataset[, antibody]') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
teddy <- iaa_teddy %>% filter(p.value<=.05 & term == 'dataset[, antibody]1') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% db$CAG)

print(teddy)

db <- ia2a_db %>% filter(message=='' & round(by,2)<=.05 & term == 'dataset[, antibody]') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
teddy <- ia2a_teddy %>% filter(p.value<=.05 & term == 'dataset[, antibody]1') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% db$CAG)

print(teddy)


###ever

gada_db <- read.csv('diabimmune_association_outputs_gada_subset.csv')
gada_db=gada_db  %>% mutate(by=p.adjust(p.value,method='BY'))
gada_teddy <- read.csv('teddy_association_outputs_gada_subset.csv')

iaa_db <- read.csv('diabimmune_association_outputs_iaa_subset.csv')
iaa_db=iaa_db%>% mutate(by=p.adjust(p.value,method='BY'))
iaa_teddy <- read.csv('teddy_association_outputs_iaa_subset.csv')

ia2a_db <- read.csv('diabimmune_association_outputs_ia2a_subset.csv')
ia2a_db=ia2a_db   %>% mutate(by=p.adjust(p.value,method='BY'))
ia2a_teddy <- read.csv('teddy_association_outputs_ia2a_subset.csv')

db <- gada_db %>% filter(round(by,2)<=.05 & term == 'dataset[, antibody]') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
teddy <- gada_teddy %>% filter(p.value<=.05 & term == 'dataset[, antibody]1') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% db$CAG)

print(teddy)

db <- iaa_db %>% filter(round(by,2)<=.05 & term == 'dataset[, antibody]') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
teddy <- iaa_teddy %>% filter(p.value<=.05 & term == 'dataset[, antibody]1') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% db$CAG)

print(teddy)

db <- ia2a_db %>% filter(round(by,2)<=.05 & term == 'dataset[, antibody]') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
teddy <- ia2a_teddy %>% filter(p.value<=.05 & term == 'dataset[, antibody]1') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% db$CAG)

print(teddy)
