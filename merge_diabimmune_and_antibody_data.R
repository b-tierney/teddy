#merge in old autoantibody data

library(tidyverse)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/diabimmune_data')

mdat=read.csv('diabimmune_metadata.tsv',sep='\t',stringsAsFactors = FALSE)
t1d=read.csv('t1d_metadata_from_old_analysis.csv',stringsAsFactors = FALSE)
t1d=t1d %>% select(-subjectID)
colnames(t1d)[1]='sampleID'

out=left_join(mdat,t1d,by='sampleID') 
out[out=='TRUE']=1
out[out=='FALSE']=1
out[is.na(out)]=0

write.csv(out,'diabimmune_metadata_with_aa_data.csv')