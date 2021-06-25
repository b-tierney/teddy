#process TEDDY metadata

###workflow
#load data
#remove useless columns or those with too many NAs or redundant ones
#force type for factor columns that look numeric
#for each comparison
  #select case and control samples
  #build a "condition" column marking cases vs controls
  #remove columns used to make condition if they encode redundant info
  #collapse metadata on subject ID
    #this function loops through each subject, finds factor or character columns that have multiple values for a given subject, and throws them out

library(tidyverse)

collapse_on_subject_id <- function(data){
  data = data %>% select(-Run)
  todrop=list()
  #count number of levels per value
  subjects = unique(data$SubjectID)
  #find those with more than 1
  for(s in subjects){
    data_sub_char = data %>% filter(SubjectID == s) %>% unique %>% select_if(is.character)
    data_sub_fac = data %>% filter(SubjectID == s) %>% unique %>% select_if(is.factor)
    data_sub_log = data %>% filter(SubjectID == s) %>% unique %>% select_if(is.logical)
    todrop[[s]] = map(bind_cols(data_sub_char,data_sub_fac,data_sub_log), function(x) length(unique(x))) %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% filter(.>1) %>% select(rowname) %>% unlist %>% unname
  }
  #remove columns and collapse data
  todrop = c('age_at_collection',unique(unlist(unname(todrop))))
  data = data %>% select(-all_of(todrop))
  data_num = data %>% {bind_cols(select_at(., "SubjectID"),select_if(., is.numeric))}  %>% group_by(SubjectID) %>% summarize_all(mean,na.rm=TRUE)
  data_nonnumeric = data %>% select(-c(data %>% select_if(is.numeric) %>% colnames)) %>% unique
  data = inner_join(data_num,data_nonnumeric)
  to_drop_single_val = map(data, function(x) length(unique(x))) %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% filter(.==1) %>% select(rowname) %>% unlist %>% unname
  if(length(to_drop_single_val)>0){
    data = data %>% select(-all_of(to_drop_single_val))
  }
  return(data)
}

setwd('~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/22417/for_validation/')

metadata = read.csv('teddy_metadata_20190821.csv')

#remove columns that are unnecessary (eg sampleIDs, etc)
toremove = c('X.y','MIAA_pos','dbgap_maskid','T1D_casecontrol_ind','gestational_age','IA_casecontrol_outcome','IA_casecontrol_ind','cc','mask_id','T1D_casecontrol_outcome','last_abx_type','m109_maskid','m138_maskid','num_different_formula_types','ongoing_other_formula','Antibioticsduringpregnancy','X.1','X','biospecimen_repository_sample_id','samplemaskid','Sex.y','HLA_Category.y','Breastmilk_Ever','AlphaDiv_OTUs','DMM_Cluster','sample_mask_id','X.x','Age_in_Months','fdr','delivery','birth_weight','ever_brstfed','brst_fed','maternal_diabetes','mom_antibiotic_use','mom_num_antibiotics','maternal_ab_exp','persistency_category','country','birth_month','time_to_excl_stop','probiotic_supplemental_ever','probiotic_start_week','probiotic_stop_week','breastfeeding','shannon_div','race_ethnicity','ongoing_fully_hydrolyzed_formula','ongoing_non_hydrolyzed_formula','ongoing_partially_hydrolyzed_for')
metadata = metadata %>% select(-c(all_of(toremove))) %>% rename(Sex=Sex.x,HLA_Category=HLA_Category.x,SubjectID=maskid)

#set type of columns that could be numeric or factor
metadata$hla_5grps = as.factor(metadata$hla_5grps)
metadata$hla_5grps_ref_DR44 = as.factor(metadata$hla_5grps_ref_DR44)

#all healthy vs post T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='After') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_post-t1d.csv')

#all healthy vs post seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection>=age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_post-sero.csv')

#all healthy vs MIAA
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection>=age_first_MIAA,!is.na(age_first_MIAA)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_miaa.csv')

#all healthy vs GAD
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection>=age_first_GAD,!is.na(age_first_GAD)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_gad.csv')

#all healthy vs IA2A
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection>=age_first_IA2A,!is.na(age_first_IA2A)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_IA2A.csv')

#all healthy vs any T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_any-T1D.csv')


#all healthy vs any seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_any-sero.csv')


###DO THESE FOR OVERALL AND FOR 3 MONTH, 6 MONTH,1 YEAR, 1.5 YEAR

#OVERALL

#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero.csv')

#3 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=92)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-3month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=92)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-3month.csv')

#6 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=183)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-6month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=183)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-6month.csv')

#6 MONTH - 12 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection>=183, age_at_collection<=365)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-6month-12month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection>=183, age_at_collection<=365)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-6month-12month.csv')

#12 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=365)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-12month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=365)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-12month.csv')

#12-18 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection>=365,age_at_collection<=548)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-12-18month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection>=365,age_at_collection<=548)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-12-18month.csv')

#18 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=548)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-18month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=548)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-18month.csv')

#18-24 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection>=548,age_at_collection<=730)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-18-24month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection>=548,age_at_collection<=730)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-18-24month.csv')

#24 MONTH
#all healthy vs pre T1D
controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=730)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-t1d-24month.csv')

#all healthy vs pre seroconversion
controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
data = bind_rows(cases,controls)
data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
data = data  %>% filter(age_at_collection<=730)
data = collapse_on_subject_id(data)
data$condition[data$condition>0]=1
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_pre-sero-24month.csv')

#all healthy vs HLA
data = collapse_on_subject_id(metadata)
data$condition = data$HLA_Category
data = data %>% select(-HLA_Category)
write.csv(data,'~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/TEDDY/processed_teddy_metadata_for_regression/healthy_HLA.csv')

