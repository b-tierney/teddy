#parse diabimmune validation output

library(tidyverse)

data = unlist(unname(read.csv('locs',header=F)))

output_list = list()
for(d in data){
	a = readRDS(d)
	name = strsplit(d,'/') %>% map_chr(2)
	output_list[[name]] = a %>% mutate(comparison = name) %>% filter(p.value<.1)
}

output = bind_rows(output_list)

output = output %>% select(estimate,std.error,statistic,p.value,BY,feature,comparison) 

colnames(output)[1:5] <- paste("validation-", colnames(output)[1:5], sep = "")

teddy = readRDS('../association_vibration_output/t1d_mas_output.rds')

merged = left_join(teddy,output)

saveRDS(merged,'t1d_db_merged.rds')