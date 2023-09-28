#befa hunting

library(tidyverse)
library(broom)
library(ggpubr)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/tao_genes/')

data=as_tibble(read.csv('abundance_matrix_normalized.tsv',sep='\t')) %>% drop_na()
colnames(data)[1]='Run'

metadata=as_tibble(read.csv('teddy_metadata_20190821.csv'))
metadata=metadata %>% select(dbgap_maskid,Run,age_t1d,age_at_collection) 
metadata$age_t1d=replace_na(data=metadata$age_t1d, -1) 
metadata$t1d = 0
metadata$t1d[metadata$age_at_collection>=metadata$age_t1d & metadata$age_t1d>0]=1

merged=inner_join(data,metadata)
merged=merged %>% select(-c('Run','age_t1d','age_at_collection')) %>% group_by(dbgap_maskid,t1d) %>% summarise_all(funs(mean)) 
temp=merged %>% select(-t1d,dbgap_maskid)
out=temp %>% map(~tidy(t.test(. ~ merged$t1d)))
out=out[2:length(out)]
names=names(out)
out=do.call('rbind',out)
out$sequence_name=names
write.csv(out,'tao_genes_analysis_output.csv')

###comparing just pre diabetics

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/tao_genes/')

data=as_tibble(read.csv('abundance_matrix_normalized.tsv',sep='\t')) %>% drop_na()
colnames(data)[1]='Run'

metadata=as_tibble(read.csv('teddy_metadata_20190821.csv'))
metadata=metadata[metadata$age_at_collection<min(metadata$age_t1d[!is.na(metadata$age_t1d)]),]
metadata=metadata %>% select(dbgap_maskid,t1d,Run) 

merged=inner_join(data,metadata)
merged=merged %>% select(-Run)
merged=merged %>% group_by(dbgap_maskid,t1d) %>% summarise_all(funs(mean)) 
temp=merged %>% select(-t1d)
out=temp %>% map(~tidy(t.test(. ~ merged$t1d)))
out=out[2:length(out)]
names=names(out)
out=do.call('rbind',out)
out$sequence_name=names
write.csv(out,'tao_genes_analysis_output_young_patients.csv')


###useful plots

ggplot(merged,aes(x=as.factor(merged$t1d),group=merged$t1d,y=log(merged$Aeromonas_veronii_complete+0.00001))) +geom_jitter(alpha=.5)+geom_boxplot(outlier.shape = NA)  + ylab('Aeromonas veronii')
ggsave('Aeromonas_veronii_complete.pdf',height=6,width=6,units='in')
ggplot(merged,aes(x=as.factor(merged$t1d),group=merged$t1d,y=log(merged$Enterobacter_aerogenes_complete+0.00001)))+ geom_jitter(alpha=.5)+geom_boxplot(outlier.shape = NA)  + ylab('Enterococcus gallinarum')
ggsave('Enterobacter_aerogenes_complete.pdf',height=6,width=6,units='in')
ggplot(merged,aes(x=as.factor(merged$t1d),group=merged$t1d,y=log(merged$Enterococcus_gallinarum_complete+0.00001)))+geom_jitter(alpha=.5) +geom_boxplot(outlier.shape = NA)  + ylab('Enterobacter aerogenes')

ggplot(merged,aes(x=as.factor(merged$t1d),group=merged$t1d,y=log(merged$Aeromonas_veronii_sylf+0.00001))) +geom_jitter(alpha=.5)+geom_boxplot(outlier.shape = NA) + ylab('Aeromonas veronii, sylf domain')
ggsave('Aeromonas_veronii_sylf.pdf',height=6,width=6,units='in')
ggplot(merged,aes(x=as.factor(merged$t1d),group=merged$t1d,y=log(merged$Enterococcus_gallinarum_sylf+0.00001))) +geom_jitter(alpha=.5)+geom_boxplot(outlier.shape = NA)  + ylab('Enterococcus gallinarum, sylf domain')
ggsave('Enterococcus_gallinarum_sylf.pdf',height=6,width=6,units='in')
ggplot(merged,aes(x=as.factor(merged$t1d),group=merged$t1d,y=log(merged$Enterobacter_aerogenes_sylf+0.00001))) +geom_jitter(alpha=.5)+geom_boxplot(outlier.shape = NA)  + ylab('Enterobacter aerogenes, sylf domain')
ggsave('Enterobacter_aerogenes_sylf.pdf',height=6,width=6,units='in')






