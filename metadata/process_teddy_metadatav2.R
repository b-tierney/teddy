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

setwd('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/')

metadata = read.csv('teddy_metadata_20190821.csv')

#remove columns that are unnecessary (eg sampleIDs, etc)
toremove = c('X.y','MIAA_pos','dbgap_maskid','T1D_casecontrol_ind','gestational_age','IA_casecontrol_outcome','IA_casecontrol_ind','cc','mask_id','T1D_casecontrol_outcome','last_abx_type','m109_maskid','m138_maskid','num_different_formula_types','ongoing_other_formula','Antibioticsduringpregnancy','X.1','X','biospecimen_repository_sample_id','samplemaskid','Sex.y','HLA_Category.y','Breastmilk_Ever','AlphaDiv_OTUs','DMM_Cluster','sample_mask_id','X.x','Age_in_Months','fdr','delivery','birth_weight','ever_brstfed','brst_fed','maternal_diabetes','mom_antibiotic_use','mom_num_antibiotics','maternal_ab_exp','persistency_category','country','birth_month','time_to_excl_stop','probiotic_supplemental_ever','probiotic_start_week','probiotic_stop_week','breastfeeding','shannon_div','race_ethnicity','ongoing_fully_hydrolyzed_formula','ongoing_non_hydrolyzed_formula','ongoing_partially_hydrolyzed_for')
metadata = metadata %>% select(-c(all_of(toremove))) %>% dplyr::rename(Sex=Sex.x,HLA_Category=HLA_Category.x,SubjectID=maskid)

#set type of columns that could be numeric or factor
metadata$hla_5grps = as.factor(metadata$hla_5grps)
metadata$hla_5grps_ref_DR44 = as.factor(metadata$hla_5grps_ref_DR44)


###DO THESE FOR OVERALL AND FOR 3 MONTH, 6 MONTH,1 YEAR, 1.5 YEAR

#OVERALL

# data is the metadata
# times is either '', a single time for 3, 6, 12... etc months. or a list of size 2 for 6-12 months for example
# HLA_status is either all_HLA, DR3_DR4_only, or not_DR3_DR4
# returns collapsed data frame and a corresponding mapping file
# filter_by_time is yes or no
get_metadata_healthy_pret1d = function(metadata,filter_by_time,times,HLA_status) {
  controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
  cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
  data = bind_rows(cases,controls)
  data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      data = data %>% filter(age_at_collection<=times[1])
    } else if (length(times) == 2) {
      data = data %>% filter(age_at_collection>=times[1], age_at_collection<=times[2])
    }
  }
  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "not_DR3_DR4") {
    # don't know why they are not eligible so just remove
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category!= "DR4*030X/0302*DR3*0501/0201")
  }
  # before collapsing create a mapping files of subjects to samples for microbiome mapping
  mapping_file = data %>% select(c('SubjectID','Run'))
  # just doing this to fit the column names of later script 
  colnames(mapping_file) <- c("V1","V2")
  data = collapse_on_subject_id(data)
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    write.csv(data,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression_march9_2022/healthy_pre-t1d-',time_out,"-",HLA_status,".csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression_march9_2022/healthy_pre-t1d-',time_out,"-",HLA_status,".mapping.csv",sep=""),row.names=FALSE)
  } else {
    write.csv(data,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression_march9_2022/healthy_pre-t1d-',HLA_status,".csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression_march9_2022/healthy_pre-t1d-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
  }
}

times_df = data.frame(c("",""),c(92,""),c(183,""),c(183,365),c(365,""),c(365,548),c(548,""),c(548,730),c(730,""))
times_df = t(times_df)
times_df = as.data.frame(times_df)

HLA_statuses = c("all_HLA","DR3_DR4_only","not_DR3_DR4")

for (y in 1:length(HLA_statuses)) {
  HLA_status_temp = HLA_statuses[y]
  for (x in 1:nrow(times_df)) {
    filter_time = "yes"
    my_times = times_df[x,]
    my_times = unlist(my_times)
    my_times = my_times[my_times!=""]
    my_times = as.numeric(my_times)
    if(length(my_times) == 0) {
      my_times = ""
      filter_time = "no"
    }
    get_metadata_healthy_pret1d(metadata=metadata,filter_by_time=filter_time,times=my_times,HLA_status=HLA_status_temp)
  }
}


# now get metadata for different antibodies

get_metadata_healthy_pre_antibody = function(metadata,filter_by_time,times,antibody_type) {
  if(antibody_type == "MIAA") {
    controls = metadata %>% filter(is.na(age_first_MIAA)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA)) %>% mutate(condition = 1)
  } else if (antibody_type == "GAD") {
    controls = metadata %>% filter(is.na(age_first_GAD)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD)) %>% mutate(condition = 1)
  } else if (antibody_type == "IA2A"){
    controls = metadata %>% filter(is.na(age_first_IA2A)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A)) %>% mutate(condition = 1)
  } else if (antibody_type == "seroconverters"){
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
  }
  data = bind_rows(cases,controls)
  data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      data = data %>% filter(age_at_collection<=times[1])
    } else if (length(times) == 2) {
      data = data %>% filter(age_at_collection>=times[1], age_at_collection<=times[2])
    }
  }
  # before collapsing create a mapping files of subjects to samples for microbiome mapping
  mapping_file = data %>% select(c('SubjectID','Run'))
  # just doing this to fit the column names of later script 
  colnames(mapping_file) <- c("V1","V2")
  data = collapse_on_subject_id(data)
  
  if(antibody_type == "MIAA") {
    main_cat_prefix = "healthy_pre-MIAA"
  } else if (antibody_type == "GAD") {
    main_cat_prefix = "healthy_pre-GAD"
  } else if (antibody_type == "IA2A") {
    main_cat_prefix = "healthy_pre-IA2A"
  } else if (antibody_type == "seroconverters") {
    main_cat_prefix = "healthy_pre-sero"
  }
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    write.csv(data,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_march9_2022/',main_cat_prefix,"-",time_out,"",".csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_march9_2022/',main_cat_prefix,"-",time_out,".mapping.csv",sep=""),row.names=FALSE)
  } else {
    write.csv(data,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_march9_2022/',main_cat_prefix,".csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_march9_2022/',main_cat_prefix,".mapping.csv",sep=""),row.names=FALSE)
  }
}

times_df = data.frame(c("",""),c(92,""),c(183,""),c(183,365),c(365,""),c(365,548),c(548,""),c(548,730),c(730,""))
times_df = t(times_df)
times_df = as.data.frame(times_df)

antibody_statuses = c("MIAA","GAD","IA2A","seroconverters")


for (y in 1:length(antibody_statuses)) {
  antibody_status_temp = antibody_statuses[y]
  for (x in 1:nrow(times_df)) {
    filter_time = "yes"
    my_times = times_df[x,]
    my_times = unlist(my_times)
    my_times = my_times[my_times!=""]
    my_times = as.numeric(my_times)
    if(length(my_times) == 0) {
      my_times = ""
      filter_time = "no"
    }
    get_metadata_healthy_pre_antibody(metadata=metadata,filter_by_time=filter_time,times=my_times,antibody_type=antibody_status_temp)
  }
}









#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d.csv')

#3 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection<=92)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-3month.csv')

#6 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection<=183)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-6month.csv')

#6 MONTH - 12 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection>=183, age_at_collection<=365)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-6month-12month.csv')


#12 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection<=365)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-12month.csv')

#12-18 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection>=365,age_at_collection<=548)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-12-18month.csv')

#18 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection<=548)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-18month.csv')

#18-24 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection>=548,age_at_collection<=730)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-18-24month.csv')

#24 MONTH
#all healthy vs pre T1D
#controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
#cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
#data = bind_rows(cases,controls)
#data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
#data = data  %>% filter(age_at_collection<=730)
#data = collapse_on_subject_id(data)
#data$condition[data$condition>0]=1
#write.csv(data,'~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/processed_teddy_metadata_for_regression/healthy_pre-t1d-24month.csv')
