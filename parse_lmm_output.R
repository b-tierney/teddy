#parse lmm output

library(tidyverse)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/lmm_cags/')

full_diabimmune <- read.csv('diabimmune_association_outputs_full_adjusted.csv')
full_teddy <- read.csv('teddy_association_outputs_full_adjusted.csv')

full_diabimmune_samponly <- read.csv('diabimmune_association_outputs_at_sampling_only_adjusted.csv')
full_teddy_samponly  <- read.csv('teddy_association_outputs_full_at_sampling_only.csv')

presero_diabimmune <- read.csv('diabimmune_association_outputs_presero_adjusted.csv')
presero_teddy <- read.csv('teddy_association_outputs_presero_adjusted.csv')

###full validation check

full_diabimmune_significant_samp_only <- full_diabimmune_samponly %>% filter(round(by,2)<=.05 & term == 'seroconverted_at_sampling') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
full_teddy_significant_samp_only_validated <- full_teddy_samponly %>% filter(p.value<=.05 & term == 'seroconverted_at_sampling') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% full_diabimmune_significant_samp_only$CAG)

presero_diabimmune_significant_ever <- presero_diabimmune %>% filter(by<=.05 & term == 'seroconverted_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
presero_teddy_significant_ever_validated <- presero_teddy %>% filter(p.value<=.05 & term == 'seroconverted_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% presero_diabimmune_significant_ever$CAG)

###############################diabetes

full_diabimmune_samponly <- read.csv('diabimmune_association_outputs_t1d_sampling_full_adjusted.csv')
full_teddy_samponly  <- read.csv('teddy_association_outputs_t1d_sampling_full.csv')

presero_diabimmune <- read.csv('diabimmune_association_outputs_pret1d_adjusted.csv')
presero_teddy <- read.csv('teddy_association_outputs_pret1d.csv')


full_diabimmune_significant_samp_only <- full_diabimmune_samponly %>% filter(round(by,2)<=.05 & term == 'diabetes_at_sampling') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
full_teddy_significant_samp_only_validated <- full_teddy_samponly %>% filter(round(p.value,2)<=.05 & term == 'diabetes_at_sampling') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% full_diabimmune_significant_samp_only$CAG)


presero_diabimmune_significant_ever <- presero_diabimmune %>% filter(round(by,2)<=.05 & term == 't1d_ever1') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
presero_teddy_significant_ever_validated <- presero_teddy %>% filter(round(p.value,2)<=.05 & term == 't1d_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% presero_diabimmune_significant_ever$CAG)




















#######only using seroconverted, nonsingular terms

full_diabimmune_significant_ever <- full_diabimmune %>% filter(term == 'seroconverted_ever' & message=='') %>% mutate(by=p.adjust(p.value,method='BY')) %>% filter(by<=.05) %>% select(CAG) %>% mutate_if(is.factor,as.character) 
full_teddy_significant_ever_validated <- full_teddy %>% filter(p.value<=.05 & term == 'seroconverted_ever' & message=='') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% full_diabimmune_significant_ever$CAG)

full_diabimmune_significant_at_sampling <- full_diabimmune %>% filter(term == 'seroconverted_at_sampling' & message=='') %>% mutate(by=p.adjust(p.value,method='BY')) %>% filter(by<=.05) %>% select(CAG) %>% mutate_if(is.factor,as.character) 
full_teddy_significant_at_sampling_validated <- full_teddy %>% filter(p.value<=.05 & term == 'seroconverted_at_sampling' & message=='') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% full_diabimmune_significant_at_sampling$CAG)

presero_diabimmune_significant_ever <- presero_diabimmune %>% filter(term == 'seroconverted_ever' & message=='') %>% mutate(by=p.adjust(p.value,method='BY')) %>% filter(by<=.05) %>% select(CAG) %>% mutate_if(is.factor,as.character) 
presero_teddy_significant_ever_validated <- presero_teddy %>% filter(p.value<=.05 & term == 'seroconverted_ever' & message=='') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% presero_diabimmune_significant_ever$CAG)


#######checking in reverse

full_teddy_significant_ever <- full_teddy %>% filter(by<=.05 & term == 'seroconverted_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
full_diabimmune_significant_ever_validated <- full_diabimmune %>% filter(p.value<=.05 & term == 'seroconverted_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% full_teddy_significant_ever$CAG)

full_teddy_significant_at_sampling <- full_teddy %>% filter(by<=.05 & term == 'seroconverted_at_sampling') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
full_diabimmune_significant_at_sampling_validated <- full_diabimmune %>% filter(p.value<=.05 & term == 'seroconverted_at_sampling') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% full_teddy_significant_at_sampling$CAG)

presero_teddy_significant_ever <- presero_teddy %>% filter(by<=.05 & term == 'seroconverted_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) 
presero_diabimmune_significant_ever_validated <- presero_diabimmune %>% filter(p.value<=.05 & term == 'seroconverted_ever') %>% select(CAG) %>% mutate_if(is.factor,as.character) %>% filter(CAG %in% presero_teddy_significant_ever$CAG)




