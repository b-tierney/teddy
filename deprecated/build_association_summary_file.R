#merge all association and vibration output into a single summary file

library(tidyverse)

files = list.files()
files = files[grep('healthy',files)]

output = list()
#for each folder
for(f in files){
	print(f)
	#load association data
	all_association_data = readRDS(paste(f,'/full_association_output_adjusted.rds',sep='')) %>% mutate(comparison=f) %>% filter(BY<.1)
	#create column indicating regression comparison
	if(nrow(all_association_data)!=0){
		#load vibration data
		vibs = readRDS(paste(f,'/full_vibration_output_summarized.rds',sep=''))
		#right join vibration data
		merged_data = right_join(all_association_data,vibs,by=c('feature'='dependent_feature'))
	output[[f]] = merged_data
	}
}

output_merged = bind_rows(output)

saveRDS(output_merged,'t1d_mas_output.rds')