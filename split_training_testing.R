### create training and testing sets 

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

#load in merged metadata

metadata = readRDS(paste(args[[1]],sep=''))

# randomly select teddy training/testing

tdat = metadata %>% filter(dataset == 'TEDDY')
tdat_case = tdat %>% filter(condition == 1)
tdat_control = tdat %>% filter(condition == 0)
num_to_sample = ceiling(.2 * nrow(tdat_case))
td_case_subs = sample_n(tdat_case,num_to_sample) %>% select(SubjectID) %>% unlist %>% unname
td_control_subs = sample_n(tdat_control,num_to_sample) %>% select(SubjectID) %>% unlist %>% unname

# randomly select diabimmune training/testing

ddat = metadata %>% filter(dataset == 'diabimmune')
ddat_case = ddat %>% filter(condition == 1)
ddat_control = ddat %>% filter(condition == 0)
num_to_sample = ceiling(.2 * nrow(ddat_case))
db_case_subs = sample_n(ddat_case,num_to_sample) %>% select(SubjectID) %>% unlist %>% unname
db_control_subs = sample_n(ddat_control,num_to_sample) %>% select(SubjectID) %>% unlist %>% unname

# create training/testing ids
testing_subjects = c(td_case_subs,db_case_subs,td_control_subs,db_control_subs)

folderids = read.csv(args[[2]],header=F) %>% unlist %>% unname
for(f in folderids){
	batchfiles = list.files(f)
	for(ff in batchfiles){
		loc = paste(f,'/',ff,sep='')
		# load the file
		abundance_data = readRDS(loc)
		# get the training set
		training = abundance_data %>% filter(!(SubjectID %in% testing_subjects))
		# get the testing set
		testing = abundance_data %>% filter(SubjectID %in% testing_subjects)
		# save the training/testing sets in the correct output folder
		saveRDS(training,paste(f,'/training_',ff,sep=''))
		saveRDS(testing,paste(f,'/testing_',ff,sep=''))
	}
}