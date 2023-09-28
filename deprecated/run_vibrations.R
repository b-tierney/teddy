#compute vibrations

library(tidyverse)
library(quantvoe)

args = commandArgs(trailingOnly=TRUE)

dependent_variables = readRDS('healthy_any-T1D_kt.rds')
independent_variables = readRDS('healthy_any-T1D_metadata_filtered.rds')

colnames(dependent_variables)[1]='sampleID'
colnames(independent_variables)[1]='sampleID'

dependent_variables$sampleID = as.character(dependent_variables$sampleID)
independent_variables$sampleID = as.character(independent_variables$sampleID)

temp = dependent_variables %>% select(-sampleID)
temp_sum = temp + 0.000001
temp_logged = log(temp_sum)

dependent_variables = bind_cols(independent_variables %>% select(sampleID),temp_logged)

association_output_full = readRDS('full_association_output_adjusted.rds')

features_of_interest = association_output_full %>% dplyr::filter(BY<=as.numeric(0.1)) %>% dplyr::select(.data$feature) %>% unique %>% unlist %>% unname

features_of_interest = intersect(features_of_interest,colnames(dependent_variables))

if(length(features_of_interest)>0){
	dependent_variables = dependent_variables %>% select(sampleID,features_of_interest)

	bound_data = dplyr::tibble(dependent_variables=list(dependent_variables),independent_variables=list(independent_variables),dsid=1)

	vibration_output = compute_vibrations(bound_data,primary_variable = 'condition',features_of_interest = unname(unlist(features_of_interest)),max_vibration_num = 10000,max_vars_in_model = 20)

	saveRDS(vibration_output,paste('vibration_output_full_',args[[1]],sep=''))
}















