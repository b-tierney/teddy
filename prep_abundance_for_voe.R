#prep for voe pipeline

library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

abundance_file = args[[1]]
metadata_file = args[[2]]
mapping_file = args[[3]]

system(paste('mkdir',gsub('.csv','',metadata_file)))
outputname = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),'_',gsub('.csv','.rds',abundance_file),sep='')
outputname_metadata = paste(gsub('.csv','',metadata_file),'/',gsub('.csv','',metadata_file),'_metadata_filtered.rds',sep='')

print('Output loc:')
print(outputname)

average_feature = function(feature_name,data){
	data = data %>% select(SubjectID,feature_name) %>% group_by(SubjectID) %>% summarise(.groups='keep',gene = mean(get(feature_name))) %>% distinct
	colnames(data)[2]=feature_name
	data = data %>% ungroup %>% arrange(SubjectID)
	return(data)
}

print('Loading data...')

d = read.csv(abundance_file) %>% select(-X)
metadata = read.csv(metadata_file,header=TRUE)
mapping = read.csv(mapping_file,header=F,sep='\t') %>% select(V1, V2) %>% rename(SubjectID = V2,sampleID = V1)
d = d %>% as_tibble %>% column_to_rownames('genename')
d = t(d) %>% as.data.frame %>% rownames_to_column('sampleID')

d_merged = inner_join(mapping,d) %>% select(-sampleID)

print('Averaging gene data...')

output = list()
genes=colnames(d_merged %>% select(-SubjectID))
count=0
for(g in genes){
	count= count +1
	out = average_feature(g,d_merged) 
	subs = out %>% select(SubjectID) %>% arrange(SubjectID)
	output[[g]] = out %>% select(-SubjectID)
}

print('Writing to file...')

averaged_data = bind_cols(subs,output)

averaged_data = averaged_data %>% filter(SubjectID %in% metadata$SubjectID)
metadata = metadata %>% filter(SubjectID %in% averaged_data$SubjectID)

saveRDS(object = averaged_data,outputname)
saveRDS(object = metadata %>% select(-X),outputname_metadata)

print('done')