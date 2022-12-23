#prep for voe pipeline

library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

abundance_file = args[[1]]
metadata_file = args[[2]]
mapping_file = args[[3]]
training_subject_file = args[[4]]
testing_subject_file = args[[5]]
metadata_suffix = args[[6]]
system(paste('mkdir -p',gsub('.csv','',metadata_file)))
outputname_train = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),'_','train_',gsub('.csv','.rds',basename(abundance_file)),sep='')
outputname_test = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),'_','test_',gsub('.csv','.rds',basename(abundance_file)),sep='')
outputname_metadata_train = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),"_train_",metadata_suffix,'_metadata_filtered.rds',sep='')
outputname_metadata_test = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),"_test_",metadata_suffix,'_metadata_filtered.rds',sep='')

print('Output loc:')

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
training = read.table(training_subject_file,header=FALSE)
training = training[,1]
testing = read.table(testing_subject_file,header=FALSE)
testing = testing[,1]

sampleNames_with_mapping_info = intersect(mapping$sampleID,colnames(d))
columns_to_keep = match(sampleNames_with_mapping_info,colnames(d))
# keep columns with metadata and the gene name column
d = d[,c(1,columns_to_keep)]
# change column names to sample IDs
subjects_to_keep = mapping[match(colnames(d)[-1], mapping$sampleID),"SubjectID"]
new_colnames = c(colnames(d)[1],subjects_to_keep)
colnames(d) = new_colnames
# training subjects to keep
training = training[training%in%subjects_to_keep]
# testing subjects to keep
testing = testing[testing%in%subjects_to_keep]


d_dt = as.data.table(d)

d_melt = data.table::melt(d_dt,id.vars=c("genename"))
d_melt_avg = d_melt[,mean(value),by=.(genename,variable)]
d_avg_long = dcast(d_melt_avg,variable ~ genename, value.var="V1")
colnames(d_avg_long)[1] = "SubjectID"
d_avg_long[,SubjectID:=as.character(SubjectID)]

d_avg_long = as.data.frame(d_avg_long)
rownames(d_avg_long) = d_avg_long$SubjectID
d_avg_long_train = d_avg_long[training,]
d_avg_long_test = d_avg_long[testing,]

train_subjects_keep = intersect(d_avg_long_train$SubjectID,metadata$SubjectID)
test_subjects_keep = intersect(d_avg_long_test$SubjectID,metadata$SubjectID)

metadata_train = metadata[match(train_subjects_keep,metadata$SubjectID),]
metadata_test = metadata[match(test_subjects_keep,metadata$SubjectID),]

d_avg_long_train = d_avg_long_train[train_subjects_keep,]
d_avg_long_test = d_avg_long_test[test_subjects_keep,]

d_avg_long_train = as_tibble(d_avg_long_train)
d_avg_long_test = as_tibble(d_avg_long_test)

saveRDS(object = d_avg_long_train,outputname_train)
saveRDS(object = d_avg_long_test,outputname_test)

saveRDS(object = metadata_train %>% select(-X),outputname_metadata_train)
saveRDS(object = metadata_test %>% select(-X),outputname_metadata_test)

print('done')
