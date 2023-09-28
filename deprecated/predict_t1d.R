#attempt to predict diabetes onset

library(caret)
library(tidyverse)


average_feature = function(feature_name,data){
	data = data %>% select(SubjectID,feature_name) %>% group_by(SubjectID) %>% summarise(.groups='keep',gene = mean(get(feature_name))) %>% distinct
	colnames(data)[2]=feature_name
	data = data %>% ungroup %>% arrange(SubjectID)
	return(data)
}

#read in abundance data
abundances = read.csv('/n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/sig_gene_analysis/no_hla_no_post_TEDDY_diabimmune_sig_genes_norm_abundances.csv')
abundances = abundances %>% column_to_rownames('genename') 
abundances = t(abundances)
abundances = abundances %>% as.data.frame %>% rownames_to_column('sampleID')

#read in teddy metadata
teddymdat = readRDS('./teddy/healthy_any-T1D/healthy_any-T1D_metadata_filtered.rds')
teddymdat = teddymdat %>% select(SubjectID,condition)

#read in diabimmune metadata
dbmdat = readRDS('./diabimmune/healthy_any-T1D/healthy_any-T1D_metadata_filtered.rds')
dbmdat = dbmdat %>% select(SubjectID,condition)

#merge and filter metadata
mdat = bind_rows(teddymdat,dbmdat)

#read in sample subject mapping
mapping = readRDS('subject_sample_mapping.rds') 

d_merged = inner_join(abundances,mapping) %>% select(-sampleID)

#average data
output = list()
genes=colnames(d_merged %>% select(-SubjectID))
count=0
for(g in genes){
	count= count +1
	out = average_feature(g,d_merged) 
	subs = out %>% select(SubjectID) %>% arrange(SubjectID)
	output[[g]] = out %>% select(-SubjectID)
}

averaged_data = bind_cols(subs,output)

#read in teddy output
teddyout=readRDS('t1d_mas_output.rds')

genes_to_predict = teddyout %>% filter(comparison == 'healthy_pre-t1d') %>% select(feature) %>% unlist %>% unname

#abundances
averaged_data_subset = averaged_data %>% select(SubjectID,genes_to_predict)

#predict based on abundances alone for all t1d associated genes
db_abundances = averaged_data_subset %>% filter(SubjectID %in% dbmdat$SubjectID)
teddy_abundances = averaged_data_subset %>% filter(SubjectID %in% teddymdat$SubjectID)

db_abundances = inner_join(db_abundances,dbmdat) %>% select(-SubjectID)
teddy_abundances = inner_join(teddy_abundances,teddymdat) %>% select(-SubjectID)

rf_fit_teddy <- train(as.factor(condition) ~ ., data = teddy_abundances)

# predict the outcome on diabimmune
pred_db <- predict(rf_fit_teddy, db_abundances %>% select(-condition))
# compare predicted outcome and true outcome
confusionMatrix(pred_db, as.factor(db_abundances$condition))





