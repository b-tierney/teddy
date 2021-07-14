#combine diabimmune and TEDDY abundances

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# diabimmune parsed data locs = ../diabimmune_alignment_output/parsed_data/

# teddy parsed data locs = ../parsed_data/

# load in a file containing the names of the batched files
suffixes = read.csv('suffixes') %>% unlist %>% unname

# load in a file containing the prefixes of the analyses
p = args[[1]]

# for each prefix + batch combination, load both merge and save in a new location
system(paste('mkdir',p))
for(s in suffixes){
	db_file = paste('../diabimmune_alignment_output/parsed_data/',p,'/',p,'_',s,'.rds',sep='')
	td_file = paste('../parsed_data/',p,'/',p,'_',s,'.rds',sep='')
	db_file_data = readRDS(db_file)
	td_file_data = readRDS(td_file)
	merged_data = bind_rows(td_file_data,db_file_data)
	saveRDS(merged_data,paste(p,'/',p,'_',s,'.rds',sep=''))
}