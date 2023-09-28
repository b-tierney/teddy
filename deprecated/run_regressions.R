library(dplyr)
library(quantvoe)

args = commandArgs(trailingOnly=TRUE)

dependent_variables = readRDS(args[[1]])
independent_variables = readRDS(args[[2]]) %>% select(-X)

colnames(dependent_variables)[1]='sampleID'
colnames(independent_variables)[1]='sampleID'

dependent_variables$sampleID = as.character(dependent_variables$sampleID)
independent_variables$sampleID = as.character(independent_variables$sampleID)

temp = dependent_variables %>% select(-sampleID)
temp_sum = temp + 0.000001
temp_logged = log(temp_sum)

dependent_variables = bind_cols(independent_variables %>% select(sampleID),temp_logged)

bound_data = dplyr::tibble(dependent_variables=list(dependent_variables),independent_variables=list(independent_variables),dsid=1)

association_output_full <- compute_initial_associations(bound_data, 'condition')

saveRDS(association_output_full,paste('association_output_full_',args[[1]],args[[2]],sep=''))