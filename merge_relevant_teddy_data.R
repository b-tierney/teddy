library(dplyr)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/TEDDY_V17/TEDDY_DATA/Analysis_Datasets')

m109 = read.sas7bdat('m_109_tvatanen_niddk_31may2012.sas7bdat')
m138 = read.sas7bdat('m_138_cstewart_31may2012_wgs_1.sas7bdat')
mapping = read.sas7bdat('../../../22417/mp109_138_dbgap_mapping.sas7bdat')
sra = read.csv('../../../22417/SraRunTable.txt',sep='\t')
sra = sra %>% filter(sra$Assay_Type=='WGS')

mapping_m109 = merge(mapping,m109,by.x='mp109_maskid',by.y='mask_id',how='inner')
mapping_m138 = merge(mapping,m138,by.x='mp138_maskid',by.y='maskid',how='left')

#check which individuals have the same number of samples
a=as.data.frame(table(mapping_m109$dbgap_maskid),stringsAsFactors = F)
a=a %>% arrange(a$Var1)
b=as.data.frame(table(sra$submitted_subject_id),stringsAsFactors = F)
b = b%>% filter(b$Var1 %in% a$Var1)
b=b %>% arrange(b$Var1)

which(a$Freq==b$Freq)