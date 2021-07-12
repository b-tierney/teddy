#parse vibration output

library(tidyverse)
library(quantvoe)

files=list.files()
files = files[grep('vibration_output',files)]

data = list()
summarized_data = list()
for(f in files){
	data[[f]] = readRDS(f)[[1]]
	summarized_data[[f]] = quantvoe::analyze_voe_data(vibration_output = readRDS(f),constant_adjusters=NULL,confounder_analysis=FALSE)
}

data = bind_rows(data)
summarized_data = bind_rows(map(summarized_data, function(x) x[[1]]))

saveRDS(data,'full_vibration_output.rds')
saveRDS(summarized_data,'full_vibration_output_summarized.rds')