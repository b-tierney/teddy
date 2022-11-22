#prep for voe pipeline

library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

abundance_file = args[[1]]
metadata_file = args[[2]]
mapping_file = args[[3]]
metadata_suffix = args[[4]]
system(paste('mkdir -p',gsub('.csv','',metadata_file)))
outputname = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),'_',gsub('.csv','.rds',abundance_file),sep='')
outputname_metadata = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),"_",metadata_suffix,'_metadata_filtered.rds',sep='')

print('Output loc:')
print(outputname)

print('Loading data...')

d_small = fread(abundance_file,sep=",",header=TRUE,data.table=FALSE,nrow=1) %>% select(-V1)
## aa.csv will have an extra row for the gene names so lets remove that
has_extra_row = d_small[1,1]=="genename"
if(has_extra_row == TRUE) {
  d = fread(abundance_file,sep=",",header=TRUE,data.table=FALSE,skip=1)
  d = d[,-1]
} else {
  d = fread(abundance_file,sep=",",header=TRUE,data.table=FALSE) %>% select(-V1)
}


metadata = read.csv(metadata_file,header=TRUE)
mapping = read.csv(mapping_file,header=T,sep=',') %>% select(V1, V2) %>% rename(SubjectID = V1,sampleID = V2)

sampleNames_with_mapping_info = intersect(mapping$sampleID,colnames(d))
columns_to_keep = match(sampleNames_with_mapping_info,colnames(d))
# keep columns with metadata and the gene name column
d = d[,c(1,columns_to_keep)]
# change column names to sample IDs
subjects_to_keep = mapping[match(colnames(d)[-1], mapping$sampleID),"SubjectID"]
new_colnames = c(colnames(d)[1],subjects_to_keep)
colnames(d) = new_colnames

d_dt = as.data.table(d)

d_melt = data.table::melt(d_dt,id.vars=c("genename"))
d_melt_avg = d_melt[,mean(value),by=.(genename,variable)]
d_avg_long = dcast(d_melt_avg,variable ~ genename, value.var="V1")
colnames(d_avg_long)[1] = "SubjectID"
d_avg_long[,SubjectID:=as.character(SubjectID)]

d_avg_long = as.data.frame(d_avg_long)
d_avg_long = d_avg_long %>% filter(SubjectID %in% metadata$SubjectID)
d_avg_long = as_tibble(d_avg_long)
metadata = metadata %>% filter(SubjectID %in% d_avg_long$SubjectID)
saveRDS(object = d_avg_long,outputname)
saveRDS(object = metadata %>% select(-X),outputname_metadata)
print('done')
