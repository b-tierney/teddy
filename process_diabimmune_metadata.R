#process diabimmune metadata

#SOMETHING IS WRONG WITH DIABETES AT SAMPLING VARIABLE, PROBABLY NEED TO REMAKE

library(tidyverse)

collapse_on_subject_id <- function(data){
  data = data %>% select(-sampleID)
  todrop=list()
  #count number of levels per value
  subjects = unique(data$SubjectID)
  #find those with more than 1
  for(s in subjects){
    data_sub_char = data %>% filter(SubjectID == s) %>% unique %>% select_if(is.character)
    data_sub_fac = data %>% filter(SubjectID == s) %>% unique %>% select_if(is.factor)
    todrop[[s]] = map(bind_cols(data_sub_char,data_sub_fac), function(x) length(unique(x))) %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% filter(.>1) %>% select(rowname) %>% unlist %>% unname
  }
  #remove columns and collapse data
  todrop = unique(unlist(unname(todrop)))
  if(length(todrop)>0){
    data = data %>% select(-all_of(todrop))
  }
  #also find columns that contain only one unique value
  to_drop_single_val = map(data, function(x) length(unique(x))) %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% filter(.==1) %>% select(rowname) %>% unlist %>% unname
  if(length(to_drop_single_val)>0){
    data = data %>% select(-all_of(to_drop_single_val))
  }
  return(data)
}

setwd('~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/diabimmune_data/')

metadata = read.csv('diabimmune_metadata_with_aa_data.csv')

#remove columns that are unnecessary (eg sampleIDs, etc)
toremove = c('Total_Reads','X')
metadata = metadata %>% select(-c(all_of(toremove)))  %>% rename(SubjectID = subjectID)

#set type of columns that could be numeric or factor
metadata$HLA_risk_class = as.factor(metadata$HLA_risk_class)

#all healthy vs post T1D
controls = metadata %>% filter(diabetes_at_sampling == 0,t1d_ever==0) %>% mutate(condition = 0)
cases = metadata %>% filter(diabetes_at_sampling == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_post-t1d.csv')

#all healthy vs post seroconversion
controls = metadata %>% filter(seroconverted_at_sampling == 0,seroconverted_ever==0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_at_sampling == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_at_sampling','seroconverted_ever','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_post-sero.csv')

#all healthy vs IAA
controls = metadata %>% filter(IAA == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(IAA_at_sampling == 1,IAA == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('IAA','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_miaa.csv')

#all healthy vs GAD
controls = metadata %>% filter(GADA == 0 ) %>% mutate(condition = 0)
cases = metadata %>% filter(GADA == 1, GADA_at_sampling == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('GADA','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_gad.csv')

#all healthy vs IA2A
controls = metadata %>% filter(IA2A==0) %>% mutate(condition = 0)
cases = metadata %>% filter(IA2A==1 , IA2A_at_sampling == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('IA2A','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_IA2A.csv')

#all healthy vs any T1D
controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_any-T1D.csv')


#all healthy vs any seroconversion
controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_any-sero.csv')


###DO THESE FOR OVERALL AND FOR 3 MONTH, 6 MONTH,1 YEAR, 1.5 YEAR

#OVERALL

#all healthy vs pre T1D
controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(diabetes_at_sampling == 0, t1d_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-t1d.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_ever == 1, seroconverted_at_sampling == 0) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-sero.csv')

#3 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
#cases = metadata %>% filter(diabetes_at_sampling == 0, t1d_ever == 1) %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
#data = data  %>% filter(age<=92)
#data = collapse_on_subject_id(data)
#write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-t1d-3month.csv')

#all healthy vs pre seroconversion
#controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
#cases = metadata %>% filter(seroconverted_ever == 1, seroconverted_at_sampling == 0) %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
#data = data  %>% filter(age<=92)
#data = collapse_on_subject_id(data)
#write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-sero-3month.csv')

#6 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(diabetes_at_sampling == 0, t1d_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age<=183)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-t1d-6month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_ever == 1, seroconverted_at_sampling == 0) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age<=183)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-sero-6month.csv')

#6 MONTH - 12 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(diabetes_at_sampling == 0, t1d_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age>=183, age<=365)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-t1d-6month-12month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_ever == 1, seroconverted_at_sampling == 0) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age>=183, age<=365)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-sero-6month-12month.csv')

#12 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(diabetes_at_sampling == 0, t1d_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age<=365)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-t1d-12month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_ever == 1, seroconverted_at_sampling == 0) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age<=365)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-sero-12month.csv')

#18 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(diabetes_at_sampling == 0, t1d_ever == 1) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_ever','diabetes_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age<=548)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-t1d-18month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(seroconverted_ever == 0) %>% mutate(condition = 0)
cases = metadata %>% filter(seroconverted_ever == 1, seroconverted_at_sampling == 0) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('seroconverted_ever','seroconverted_at_sampling','GADA_at_sampling','ZNT8A_at_sampling','IAA_at_sampling','ICA_at_sampling','IA2A_at_sampling'))
data = data  %>% filter(age<=548)
data = collapse_on_subject_id(data)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_diabimmune_metadata_for_regression/healthy_pre-sero-18month.csv')

