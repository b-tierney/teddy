#parse association output

library(tidyverse)

files = list.files()
files = files[grepl('association_output_full',files)]

data = list()
for(f in files){
	data[[f]] = readRDS(f)$output
}

data_bound = bind_rows(data)

data_bound = data_bound %>% mutate(BY = p.adjust(p.value,method='BY'),BH = p.adjust(p.value,method='BH'),BY = p.adjust(p.value,method='bonferroni'))

saveRDS(data_bound,'full_association_output_adjusted.rds')