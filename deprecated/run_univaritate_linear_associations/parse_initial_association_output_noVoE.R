args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(qvalue)

input_folder=args[[1]]

prefix_name=basename(input_folder)

files=list.files(input_folder,full.names=TRUE)
files = files[grepl('association_output_full_healthy',files)]

data = list()
for(f in files){
data[[f]] = readRDS(f)
}

data_bound = bind_rows(data)
# calculate adjusted pvalues with all genes
data_bound = data_bound %>% mutate(BY_allgenes = p.adjust(p.value,method='BY'),BH_allgenes = p.adjust(p.value,method='BH'),bonferroni_allgenes = p.adjust(p.value,method='bonferroni'),qvalue=qvalue(p=p.value)$qvalues)
# remove the old adjusted pvalues
data_bound = data_bound %>% select(-c(bonferroni,BH,BY))
# rename columns
data_bound = data_bound %>% dplyr::rename(BY=BY_allgenes,BH=BH_allgenes,bonferroni=bonferroni_allgenes)

saveRDS(data_bound,paste(input_folder,'/',prefix_name,'_full_association_output_adjusted.rds',sep=''))

# calculate number of significant genes
sig_0.1_BY = sum(data_bound$BY < 0.1)
sig_0.1_bonferroni = sum(data_bound$bonferroni < 0.1)
sig_0.1_qvalue = sum(data_bound$qvalue < 0.1)

print(prefix_name)
print(paste("There are",sig_0.1_BY,"genes with BY less than 0.1"))
print(paste("There are",sig_0.1_bonferroni,"genes with bonferroni less than 0.1"))
print(paste("There are",sig_0.1_qvalue,"genes with qvalue less than 0.1"))
