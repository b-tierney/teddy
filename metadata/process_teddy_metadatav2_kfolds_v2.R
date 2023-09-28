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
  data = data %>% dplyr::select(-Run)
  #todrop=list()
  #count number of levels per value
  subjects = unique(data$SubjectID)
  todrop <- vector(mode = "list", length = length(subjects))
  #find those with more than 1
  todrop = lapply(subjects, function(s) {
    data_single_subj = unique(data[data$SubjectID%in%s,])
    col_classes = sapply(data_single_subj,class)
    data_sub_char = data_single_subj[,col_classes=="character"]
    data_sub_fac = data_single_subj[,col_classes=="factor"]
    data_sub_log = data_single_subj[,col_classes=="logical"]
    all_data_cols = cbind(data_sub_char,data_sub_fac,data_sub_log)
    unique_counts_per_col = apply(all_data_cols,2, function(mycol) length(unique(mycol)))
    return(names(unique_counts_per_col[unique_counts_per_col>1]))
  })
  #remove columns and collapse data
  todrop = c('age_at_collection',unique(unlist(unname(todrop))))
  data = data %>% dplyr::select(-all_of(todrop))
  data_num = data %>% {bind_cols(select_at(., "SubjectID"),select_if(., is.numeric))}  %>% group_by(SubjectID) %>% summarize_all(mean,na.rm=TRUE)
  data_nonnumeric = data %>% dplyr::select(-c(data %>% select_if(is.numeric) %>% colnames)) %>% unique
  data = inner_join(data_num,data_nonnumeric)
  to_drop_single_val = map(data, function(x) length(unique(x))) %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% filter(.==1) %>% dplyr::select(rowname) %>% unlist %>% unname
  if("getsCondition"%in%to_drop_single_val) {
    to_drop_single_val = to_drop_single_val[-match("getsCondition",to_drop_single_val)]
  }
  if("condition"%in%to_drop_single_val) {
    to_drop_single_val = to_drop_single_val[-match("condition",to_drop_single_val)]
  }
  
  if(length(to_drop_single_val)>0){
    data = data %>% dplyr::select(-all_of(to_drop_single_val))
  }
  return(data)
}

calculate_autoantibody_number = function(data) {
  data_sub = data %>% group_by(SubjectID) %>% slice_max(age_at_collection,n=1,with_ties = FALSE) %>% dplyr::select(c(age_at_collection,age_first_MIAA,age_first_GAD,age_first_IA2A)) %>% column_to_rownames(var = "SubjectID")
  number_autoantibodies_each_subject = apply(data_sub, 1, function(x) {
    oldest_age = x[1]
    had_MIAA <- oldest_age >= as.numeric(x["age_first_MIAA"])
    if(is.na(had_MIAA)) {
      had_MIAA = FALSE
    }
    had_GAD <- oldest_age >= as.numeric(x["age_first_GAD"])
    if(is.na(had_GAD)) {
      had_GAD = FALSE
    }
    had_IA2A <- oldest_age >= as.numeric(x["age_first_IA2A"])
    if(is.na(had_IA2A)) {
      had_IA2A = FALSE
    }
    num_autoantibodies_had = sum(as.numeric(c(had_MIAA,had_GAD,had_IA2A)))
    return(num_autoantibodies_had)
  })
  data$number_autoantibodies = number_autoantibodies_each_subject[match(data$SubjectID,names(number_autoantibodies_each_subject))]
  return(data)
}

setwd('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/')

metadata = read.csv('../analysis_old_assembly_alignment/teddy_old_analysis/teddy_metadata_20190821_with_GRS2.csv')
# don't include metadata we don't have GRS scores for
metadata = metadata[!is.na(metadata$GRS_score),]
# also only include samples we have abundance data for
abundance_samples = read.table("sample_we_have_abundance_data_for.txt",header=FALSE)
abundance_samples = abundance_samples[,1]
abundance_metadata_samples = intersect(metadata$Run,abundance_samples)
metadata = metadata[match(abundance_metadata_samples,metadata$Run),]

#remove columns that are unnecessary (eg sampleIDs, etc)
toremove = c('X.y','MIAA_pos','dbgap_maskid','T1D_casecontrol_ind','gestational_age','IA_casecontrol_outcome','IA_casecontrol_ind','cc','mask_id','T1D_casecontrol_outcome','last_abx_type','m109_maskid','m138_maskid','num_different_formula_types','ongoing_other_formula','Antibioticsduringpregnancy','X.1','X','biospecimen_repository_sample_id','samplemaskid','Sex.y','HLA_Category.y','AlphaDiv_OTUs','DMM_Cluster','sample_mask_id','X.x','Age_in_Months','delivery','birth_weight','ever_brstfed','maternal_diabetes','mom_antibiotic_use','mom_num_antibiotics','maternal_ab_exp','persistency_category','country','birth_month','time_to_excl_stop','probiotic_supplemental_ever','probiotic_start_week','probiotic_stop_week','breastfeeding','shannon_div','race_ethnicity','ongoing_fully_hydrolyzed_formula','ongoing_non_hydrolyzed_formula','ongoing_partially_hydrolyzed_for')
metadata = metadata %>% dplyr::select(-c(all_of(toremove))) %>% dplyr::rename(Sex=Sex.x,HLA_Category=HLA_Category.x,SubjectID=maskid,grs2=GRS_score)

#set type of columns that could be numeric or factor
metadata$hla_5grps = as.factor(metadata$hla_5grps)
metadata$hla_5grps_ref_DR44 = as.factor(metadata$hla_5grps_ref_DR44)
metadata$fdr = as.numeric(metadata$fdr)

###DO THESE FOR OVERALL AND FOR 3 MONTH, 6 MONTH,1 YEAR, 1.5 YEAR

#OVERALL
library(caret)
library(data.table)
# data is the metadata
# times is either '', a single time for 3, 6, 12... etc months. or a list of size 2 for 6-12 months for example
# HLA_status is either all_HLA, DR3_DR4_only, or not_DR3_DR4
# returns collapsed data frame and a corresponding mapping file
# filter_by_time is yes or no
get_metadata_healthy_pret1d = function(metadata,filter_by_time,times,HLA_status,proportion_train) {
  controls = metadata %>% filter(t1d == FALSE) %>% mutate(getsCondition = 0)
  censored_time = as.data.frame(as.data.table(controls)[,max(age_at_collection),by=SubjectID])
  controls$time = censored_time[match(controls$SubjectID,censored_time$SubjectID),"V1"]
  cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(getsCondition = 1)
  cases$time = cases$age_t1d
  data = bind_rows(cases,controls)

  #data = data %>% select(-c('T1D_Outcome','t1d','age_t1d','t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      # keep individuals that are younger than baseline
      data = data %>% filter(age_at_collection<=times[1])
      # only keep individuals who get disease or are censored after baseline though.
      data = data %>% filter(time>=times[1])
    } else if (length(times) == 2) {
      data = data %>% filter(age_at_collection>=times[1], age_at_collection<=times[2])
      data = data %>% filter(time>=times[2])
    }
  }
  # for each subject get number of autoantibodies at closest time to times[1]
  data = calculate_autoantibody_number(data)
  
  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR4_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR4*030X/0302")
  } else if (HLA_status == "DR4_DR8_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR8*0401/0402")
  } else if (HLA_status == "DR3_DR3_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR3*0501/0201*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR1_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR1*0101/0501")
  } else if (HLA_status == "DR4_DR13") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR13*0102/0604")
  }
  #data_subj_condition_orig = unique(data[,c("SubjectID","getsCondition")])
  #ctrl_num = sum(data_subj_condition_orig$getsCondition == 0)
  #case_num = sum(data_subj_condition_orig$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  #if(ctrl_num <= 9 | case_num <= 3) {
  #  return(c("healthy-preT1D",HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  #}
  
  # before collapsing create a mapping files of subjects to samples for microbiome mapping
  mapping_file = data %>% select(c('SubjectID','Run'))
  # just doing this to fit the column names of later script 
  colnames(mapping_file) <- c("V1","V2")
  data = collapse_on_subject_id(data)
  data$getsCondition = as.factor(data$getsCondition)
  
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 2) {
    return(c("healthy-preT1D",HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  
  set.seed(3456)
  trainIndex <- createDataPartition(data$getsCondition, p = proportion_train, 
                                    list = FALSE, 
                                    times = 1)
  train_df = data[trainIndex,]
  train_subjects <- data[ trainIndex,]$SubjectID
  test_df  <- data[-trainIndex,]
  folds <- createFolds(train_df$getsCondition, k=3, 
                       list = TRUE)
  
  train_subjects_1 = train_df[folds[[1]],c("SubjectID","getsCondition")]
  train_subjects_2 = train_df[folds[[2]],c("SubjectID","getsCondition")]
  train_subjects_3 = train_df[folds[[3]],c("SubjectID","getsCondition")]
  test_subjects = test_df[,c("SubjectID","getsCondition")]
  
  # also scale grs2 score
  grs2_df = train_df[,"grs2",drop=FALSE]
  preProcValues_grs <- preProcess(grs2_df, method = c("center","scale"))
  trainTransformed <- predict(preProcValues_grs, train_df[,"grs2",drop=FALSE])
  testTransformed <- predict(preProcValues_grs, test_df[,"grs2",drop=FALSE])
  train_df$grs2 = trainTransformed$grs2
  test_df$grs2 = testTransformed$grs2
  
  ctrl_num_train = sum(train_df$getsCondition == 0)
  case_num_train = sum(train_df$getsCondition == 1)
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num_train < 3 | case_num_train < 1) {
    return(c("healthy-preT1D",HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    } else {
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/healthy_pre-t1d-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    }
  return(c("healthy-preT1D",HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"willRun"))
}

times_df = data.frame(c("",""),c(92,""),c(183,""),c(183,365),c(365,""),c(365,548),c(548,""),c(548,730),c(730,""))
times_df = t(times_df)
times_df = as.data.frame(times_df)

HLA_statuses = c("all_HLA","DR3_DR4_only","DR4_DR4_only","DR4_DR8_only","DR3_DR3_only","DR4_DR1_only","DR4_DR13")

proportion_train_list = c(0.66,0.50)

train_test_nums_pre_T1D = matrix("",ncol=7,nrow=length(HLA_statuses)*nrow(times_df)*length(proportion_train_list))
counter = 1
for (y in 1:length(HLA_statuses)) {
  print(y)
  HLA_status_temp = HLA_statuses[y]
  for (x in 1:nrow(times_df)) {
    print(x)
    filter_time = "yes"
    my_times = times_df[x,]
    my_times = unlist(my_times)
    my_times = my_times[my_times!=""]
    my_times = as.numeric(my_times)
    if(length(my_times) == 0) {
      my_times = ""
      filter_time = "no"
    }
    for(j in 1:length(proportion_train_list)) {
      proportion_train = proportion_train_list[j]
      train_test_nums = get_metadata_healthy_pret1d(metadata=metadata,filter_by_time=filter_time,times=my_times,HLA_status=HLA_status_temp,proportion_train=proportion_train)
      train_test_nums_pre_T1D[counter,] = train_test_nums
      counter = counter + 1
    }
  }
}
colnames(train_test_nums_pre_T1D) = c("condition","HLA","time","proportion_samples_train","case_sample_size","ctrl_sample_size","will_or_wont_run_model")
train_test_nums_pre_T1D = as.data.frame(train_test_nums_pre_T1D)
# now get metadata for different antibodies

get_metadata_healthy_pre_antibody = function(metadata,filter_by_time,times,antibody_type,HLA_status,proportion_train) {
  if(antibody_type == "MIAA") {
    controls = metadata %>% filter(is.na(age_first_MIAA)) %>% mutate(getsCondition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA)) %>% mutate(getsCondition = 1)
    cases$time = cases$age_first_MIAA
  } else if (antibody_type == "GAD") {
    controls = metadata %>% filter(is.na(age_first_GAD)) %>% mutate(getsCondition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD)) %>% mutate(getsCondition = 1)
    cases$time = cases$age_first_GAD
  } else if (antibody_type == "IA2A"){
    controls = metadata %>% filter(is.na(age_first_IA2A)) %>% mutate(getsCondition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A)) %>% mutate(getsCondition = 1)
    cases$time = cases$age_first_IA2A
  } else if (antibody_type == "seroconverters"){
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(getsCondition = 0)
    cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(getsCondition = 1)
    cases$time = cases$age_mult_persist
  } else if (antibody_type == "triple_converters_vs_T1D") {
    controls = metadata %>% filter(t1d_sero_control == 'seroconverted') %>% filter(three_persist_conf == TRUE) %>% mutate(getsCondition = 0)
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(getsCondition = 1)
    cases$time = cases$age_t1d
  } else if (antibody_type == "serconverters_or_T1D") {
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(getsCondition = 0)
    cases_T1D = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(getsCondition = 1)
    cases_sero = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(getsCondition = 1)
    cases_run = unique(c(cases_T1D$Run,cases_sero$Run))
    cases = metadata[match(cases_run,metadata$Run),]  %>% mutate(getsCondition = 1)
    
    case_time = data.frame(cases$age_t1d,cases$age_mult_persist)
    cases$time = apply(case_time, 1, function(x) {
      time_sero = x["cases.age_mult_persist"]
      time_t1d = x["cases.age_t1d"]
      if(!is.na(time_sero) & !is.na(time_t1d)) {
        time = min(time_sero,time_t1d)
      } else if(!is.na(time_sero) & is.na(time_t1d)) {
        time = time_sero
      } else if(is.na(time_sero) & !is.na(time_t1d)) {
        time = time_t1d
      }
      return(time)
    })
  }
  censored_time = as.data.frame(as.data.table(controls)[,max(age_at_collection),by=SubjectID])
  controls$time = censored_time[match(controls$SubjectID,censored_time$SubjectID),"V1"]
  
  
  data = bind_rows(cases,controls)
  #data = data %>% select(-c('t1d_sero_control','age_first_MIAA','age_first_GAD','age_first_IA2A','age_first_pos','age_mult_persist'))
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      data = data %>% filter(age_at_collection<=times[1])
      data = data %>% filter(time>=times[1])
    } else if (length(times) == 2) {
      data = data %>% filter(age_at_collection>=times[1], age_at_collection<=times[2])
      data = data %>% filter(time>=times[2])
    }
  }
  # for each subject get number of autoantibodies at closest time to times[1]. 
  data = calculate_autoantibody_number(data)
  
  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR4_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR4*030X/0302")
  } else if (HLA_status == "DR4_DR8_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR8*0401/0402")
  } else if (HLA_status == "DR3_DR3_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR3*0501/0201*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR1_only") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR1*0101/0501")
  } else if (HLA_status == "DR4_DR13") {
    data = data %>% filter(HLA_Category!= "Not*Eligible")
    data = data %>% filter(HLA_Category== "DR4*030X/0302*DR13*0102/0604")
  }
  
  #data_subj_condition_orig = unique(data[,c("SubjectID","condition")])
  #ctrl_num = sum(data_subj_condition_orig$condition == 0)
  #case_num = sum(data_subj_condition_orig$condition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  #if(ctrl_num <= 9 | case_num <= 3) {
  #  return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  #}
  
  # before collapsing create a mapping files of subjects to samples for microbiome mapping
  mapping_file = data %>% select(c('SubjectID','Run'))
  # just doing this to fit the column names of later script 
  colnames(mapping_file) <- c("V1","V2")
  data = collapse_on_subject_id(data)
  data$getsCondition = as.factor(data$getsCondition)
  
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 2) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  

  set.seed(3456)
  trainIndex <- createDataPartition(data$getsCondition, p = proportion_train, 
                                    list = FALSE, 
                                    times = 1)
  train_df = data[trainIndex,]
  train_subjects <- data[ trainIndex,]$SubjectID
  test_df  <- data[-trainIndex,]
  folds <- createFolds(train_df$getsCondition, k=3, 
                       list = TRUE)
  
  train_subjects_1 = train_df[folds[[1]],c("SubjectID","getsCondition")]
  train_subjects_2 = train_df[folds[[2]],c("SubjectID","getsCondition")]
  train_subjects_3 = train_df[folds[[3]],c("SubjectID","getsCondition")]
  test_subjects = test_df[,c("SubjectID","getsCondition")]
  
  # also scale grs2 score
  grs2_df = train_df[,"grs2",drop=FALSE]
  preProcValues_grs <- preProcess(grs2_df, method = c("center","scale"))
  trainTransformed <- predict(preProcValues_grs, train_df[,"grs2",drop=FALSE])
  testTransformed <- predict(preProcValues_grs, test_df[,"grs2",drop=FALSE])
  train_df$grs2 = trainTransformed$grs2
  test_df$grs2 = testTransformed$grs2
  
  ctrl_num_train = sum(train_df$getsCondition == 0)
  case_num_train = sum(train_df$getsCondition == 1)
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)

  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num_train < 3 | case_num_train < 1) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  
  if(antibody_type == "MIAA") {
    main_cat_prefix = "healthy_pre-MIAA"
  } else if (antibody_type == "GAD") {
    main_cat_prefix = "healthy_pre-GAD"
  } else if (antibody_type == "IA2A") {
    main_cat_prefix = "healthy_pre-IA2A"
  } else if (antibody_type == "seroconverters") {
    main_cat_prefix = "healthy_pre-sero"
  } else if (antibody_type == "triple_converters_vs_T1D") {
    main_cat_prefix = "triple_converters_vs_T1D"
  } else if (antibody_type == "serconverters_or_T1D") {
    main_cat_prefix = "serconverters_or_T1D"
  }
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  } else {
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"willRun"))
}

times_df = data.frame(c("",""),c(92,""),c(183,""),c(183,365),c(365,""),c(365,548),c(548,""),c(548,730),c(730,""))
times_df = t(times_df)
times_df = as.data.frame(times_df)

antibody_statuses = c("MIAA","GAD","IA2A","seroconverters","triple_converters_vs_T1D","serconverters_or_T1D")
#antibody_statuses = c("serconverters_or_T1D")
proportion_train_list = c(0.66,0.50)

train_test_nums_pre_diff_conditions = matrix("",ncol=7,nrow=length(antibody_statuses)*length(HLA_statuses)*nrow(times_df)*length(proportion_train_list))
counter = 1
for (y in 1:length(antibody_statuses)) {
  antibody_status_temp = antibody_statuses[y]
  for(k in 1:length(HLA_statuses)) {
    HLA_status = HLA_statuses[k]
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
      for(j in 1:length(proportion_train_list)) {
        proportion_train = proportion_train_list[j]
        print(y)
        print(k)
        print(x)
        print(j)
        case_ctrl_counts = get_metadata_healthy_pre_antibody(metadata=metadata,filter_by_time=filter_time,times=my_times,antibody_type=antibody_status_temp,HLA_status=HLA_status,proportion_train = proportion_train)
        train_test_nums_pre_diff_conditions[counter,] = case_ctrl_counts
        counter = counter + 1
      }
    }
  }
}
colnames(train_test_nums_pre_diff_conditions) = c("condition","HLA","time","proportion_samples_train","case_sample_size","ctrl_sample_size","will_or_wont_run_model")
train_test_nums_pre_diff_conditions = as.data.frame(train_test_nums_pre_diff_conditions)
# essentially I am getting microbiome abundances of test cases here X months before T1D onset. By comparison I am taking all samples from a given control subject
get_metadata_healthy_days_before_condition = function(metadata,filter_by_time,times,antibody_type,HLA_status,proportion_train) {
  if(antibody_type == "MIAA") {
    controls = metadata %>% filter(is.na(age_first_MIAA)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_first_MIAA
    cases$age_condition = cases$age_first_MIAA
  } else if (antibody_type == "GAD") {
    controls = metadata %>% filter(is.na(age_first_GAD)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_first_GAD
    cases$age_condition = cases$age_first_GAD
  } else if (antibody_type == "IA2A"){
    controls = metadata %>% filter(is.na(age_first_IA2A)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_first_IA2A
    cases$age_condition = cases$age_first_IA2A
  } else if (antibody_type == "seroconverters"){
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_mult_persist
    cases$age_condition = cases$age_mult_persist
  } else if (antibody_type == "triple_converters_vs_T1D") {
    controls = metadata %>% filter(t1d_sero_control == 'seroconverted') %>% filter(three_persist_conf == TRUE) %>% mutate(condition = 0)
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    controls$age_condition = max(controls$age_first_MIAA,controls$age_first_IA2A,cases$age_first_GAD)
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "T1D") {
    controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    controls$age_condition = controls$age_t1d
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "serconverters_or_T1D") {
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases_T1D = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    cases_sero = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
    cases_run = unique(c(cases_T1D$Run,cases_sero$Run))
    cases = metadata[match(cases_run,metadata$Run),]  %>% mutate(condition = 1)
    controls$age_condition = NA
    cases$age_condition = apply(cases, 1, function(myrow) min(as.numeric(myrow["age_t1d"]),as.numeric(myrow["age_mult_persist"]),na.rm = TRUE))
  }
  
  cases$diff_in_num = cases$age_condition - cases$age_at_collection
  controls$diff_in_num = controls$age_condition - controls$age_at_collection
  cases =cases[cases$diff_in_num>0,]
  
  
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      cases = cases[cases$diff_in_num <= times[1],]
    } else if (length(times) == 2) {
      cases = cases[cases$diff_in_num <= times[2] & cases$diff_in_num >= times[1],]
    }
  }

  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR4_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR4*030X/0302")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR8_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR8*0401/0402")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR3_DR3_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR3*0501/0201*DR3*0501/0201")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR1_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR1*0101/0501")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR13") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR13*0102/0604")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  }
  
  
  controls = calculate_autoantibody_number(controls)
  cases = calculate_autoantibody_number(cases)
  # get average age at collection for each subject
  time_controls = as.data.frame(as.data.table(controls)[,mean(age_at_collection),by=SubjectID])
  controls$time = time_controls[match(controls$SubjectID,time_controls$SubjectID),"V1"]
  time_cases = as.data.frame(as.data.table(cases)[,mean(age_at_collection),by=SubjectID])
  cases$time = time_cases[match(cases$SubjectID,time_cases$SubjectID),"V1"]
  
  mapping_file = rbind(cases,controls) %>% select(c('SubjectID','Run'))
  
  # combine controls and cases temporarily so I can get to subject level
  data = rbind(controls,cases)
  data = collapse_on_subject_id(data)
  
  ctrl_num = sum(data$condition == 0)
  case_num = sum(data$condition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 3) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  
  controls = data[data$condition == 0,]
  cases = data[data$condition == 1,]
  
  data_subj_condition_case_ctrl = rbind(controls,cases)
  data_subj_condition_case_ctrl$condition = as.factor(data_subj_condition_case_ctrl$condition)
  # now that I have selected my cases I will select my controls by trying to choose samples with the same distribution of ages

  set.seed(3456)
  trainIndex <- createDataPartition(data_subj_condition_case_ctrl$condition, p = proportion_train, 
                                    list = FALSE, 
                                    times = 1)
  
  train_data = data_subj_condition_case_ctrl[trainIndex,]
  test_data = data_subj_condition_case_ctrl[-trainIndex,]
  
  case_train = train_data[train_data$condition == 1,]
  case_test = test_data[test_data$condition == 1,]
  ctrl_train = train_data[train_data$condition == 0,]
  ctrl_test = test_data[test_data$condition == 0,]

  case_mean = mean(case_train$time)
  case_sd = sd(case_train$time)
  
  ctrl_train_prob_vector = pnorm(ctrl_train$time,mean=case_mean,sd=case_sd)
  ctrl_test_prob_vector = pnorm(ctrl_test$time,mean=case_mean,sd=case_sd)
  ctrl_train$prob_vector = ctrl_train_prob_vector
  ctrl_test$prob_vector = ctrl_test_prob_vector

  #controls_temp = controls
  if(length(unique(case_train$SubjectID)) <= length(unique(ctrl_train$SubjectID))) {
    ctrl_train_samps = sample(ctrl_train$SubjectID,size=length(unique(case_train$SubjectID)),prob=ctrl_train$prob_vector,replace = FALSE)
    ctrl_test_samps = sample(ctrl_test$SubjectID,size=length(unique(case_test$SubjectID)),prob=ctrl_test$prob_vector,replace = FALSE)
  } else {
    ctrl_train_samps = sample(ctrl_train$SubjectID,size=length(unique(ctrl_train$SubjectID)),prob=ctrl_train$prob_vector,replace = FALSE)
    ctrl_test_samps = sample(ctrl_test$SubjectID,size=length(unique(ctrl_test$SubjectID)),prob=ctrl_test$prob_vector,replace = FALSE)
  }
  ctrl_train = ctrl_train[match(ctrl_train_samps,ctrl_train$SubjectID),]
  ctrl_test = ctrl_test[match(ctrl_test_samps,ctrl_test$SubjectID),]
  
  ctrl_train = ctrl_train[,-match(c("diff_in_num","prob_vector"),colnames(ctrl_train))]
  ctrl_test = ctrl_test[,-match(c("diff_in_num","prob_vector"),colnames(ctrl_test))]
  case_train = case_train[,-match(c("diff_in_num"),colnames(case_train))]
  case_test = case_test[,-match(c("diff_in_num"),colnames(case_test))]
  
  train_df2 = rbind(ctrl_train,case_train)
  test_df2 = rbind(ctrl_test,case_test)
  
  all_data2 = rbind(train_df2,test_df2)

  ctrl_num = sum(train_df2$condition == 0)
  case_num = sum(train_df2$condition == 1)
  ctrl_num_all = sum(all_data2$condition == 0)
  case_num_all = sum(all_data2$condition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 1) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num_all,ctrl_num_all,"wont_run"))
  }
  

  folds <- createFolds(train_df2$condition, k=3, 
                       list = TRUE)

  train_subjects_1 = train_df2[folds[[1]],c("SubjectID","condition")]
  train_subjects_2 = train_df2[folds[[2]],c("SubjectID","condition")]
  train_subjects_3 = train_df2[folds[[3]],c("SubjectID","condition")]
  test_subjects = test_df2[,c("SubjectID","condition")]
  train_subjects = train_df2[,c("SubjectID","condition")]
  
  # also scale grs2 score
  grs2_df = train_df2[,"grs2",drop=FALSE]
  preProcValues_grs <- preProcess(grs2_df, method = c("center","scale"))
  trainTransformed <- predict(preProcValues_grs, train_df2[,"grs2",drop=FALSE])
  testTransformed <- predict(preProcValues_grs, test_df2[,"grs2",drop=FALSE])
  train_df2$grs2 = train_df2$grs2
  test_df2$grs2 = test_df2$grs2

  if(antibody_type == "MIAA") {
    main_cat_prefix = "healthy_pre-MIAA"
  } else if (antibody_type == "GAD") {
    main_cat_prefix = "healthy_pre-GAD"
  } else if (antibody_type == "IA2A") {
    main_cat_prefix = "healthy_pre-IA2A"
  } else if (antibody_type == "seroconverters") {
    main_cat_prefix = "healthy_pre-sero"
  } else if (antibody_type == "triple_converters_vs_T1D") {
    main_cat_prefix = "triple_converters_vs_T1D"
  } else if (antibody_type == "T1D") {
    main_cat_prefix = "healthy_pre-t1d"
  } else if (antibody_type == "serconverters_or_T1D") {
    main_cat_prefix = "serconverters_or_T1D"
  }
  
  colnames(test_df2)[match("condition",colnames(test_df2))] = "getsCondition"
  colnames(train_df2)[match("condition",colnames(train_df2))] = "getsCondition"
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    time_out = paste(time_out,"_before_condition",sep="")
    write.csv(test_df2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(train_df2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  } else {
    write.csv(test_df2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(train_df2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"willRun"))
}

HLA_statuses = c("all_HLA","DR3_DR4_only","DR4_DR4_only","DR4_DR8_only","DR3_DR3_only","DR4_DR1_only","DR4_DR13")

times_df = data.frame(c(30,""),c(60,""),c(90,""),c(30,60),c(60,90))
times_df = t(times_df)
times_df = as.data.frame(times_df)

antibody_statuses = c("MIAA","GAD","IA2A","seroconverters","triple_converters_vs_T1D","T1D","serconverters_or_T1D")

proportion_train_list = c(0.66,0.50)

train_test_nums_pre_diff_before_conditions = matrix("",ncol=7,nrow=length(antibody_statuses)*length(HLA_statuses)*nrow(times_df)*length(proportion_train_list))
counter = 1
for (y in 1:length(antibody_statuses)) {
  antibody_status_temp = antibody_statuses[y]
  for(k in 1:length(HLA_statuses)) {
    HLA_status = HLA_statuses[k]
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
      for(j in 1:length(proportion_train_list)) {
        proportion_train = proportion_train_list[j]
        print(y)
        print(k)
        print(x)
        print(j)
        case_ctrl_samp_counts = get_metadata_healthy_days_before_condition(metadata=metadata,filter_by_time=filter_time,times=my_times,antibody_type=antibody_status_temp,HLA_status=HLA_status,proportion_train = proportion_train)
        train_test_nums_pre_diff_before_conditions[counter,] = case_ctrl_samp_counts
        counter = counter + 1
      }
    }
  }
}
colnames(train_test_nums_pre_diff_before_conditions) = c("condition","HLA","time","proportion_samples_train","case_sample_size","ctrl_sample_size","will_or_wont_run_model")
train_test_nums_pre_diff_before_conditions = as.data.frame(train_test_nums_pre_diff_before_conditions)

train_test_nums_pre_diff_before_conditions$cat = "before_condition"
train_test_nums_pre_diff_conditions$cat = "age_subject"
train_test_nums_pre_T1D$cat = "age_subject"

train_test_sample_nums_all = rbind(train_test_nums_pre_diff_before_conditions,train_test_nums_pre_diff_conditions,train_test_nums_pre_T1D)

dim(train_test_sample_nums_all[train_test_sample_nums_all$will_or_wont_run_model == "willRun",]) # 1065

train_test_sample_nums_all_willRun = train_test_sample_nums_all[train_test_sample_nums_all$will_or_wont_run_model == "willRun",]

train_test_sample_nums_all_willRun_without_cut = train_test_sample_nums_all_willRun[,-match("proportion_samples_train",colnames(train_test_sample_nums_all_willRun))]
train_test_sample_nums_all_willRun_without_cut = unique(train_test_sample_nums_all_willRun_without_cut)

write.csv(train_test_sample_nums_all,"train_test_sample_nums_all.csv")




# do another function where I compare very close to condition vs far from getting condition

get_metadata_healthy_days_before_condition_and_earlier = function(metadata,filter_by_time,times,antibody_type,HLA_status,proportion_train) {
  if(antibody_type == "MIAA") {
    cases = metadata %>% filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA))
    cases$age_condition = cases$age_first_MIAA
  } else if (antibody_type == "GAD") {
    cases = metadata %>% filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD))
    cases$age_condition = cases$age_first_GAD
  } else if (antibody_type == "IA2A"){
    cases = metadata %>% filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A))
    cases$age_condition = cases$age_first_IA2A
  } else if (antibody_type == "seroconverters"){
    cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist))
    cases$age_condition = cases$age_mult_persist
  } else if (antibody_type == "triple_converters_vs_T1D") {
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before')
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "T1D") {
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before')
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "serconverters_or_T1D") {
    cases_T1D = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') 
    cases_sero = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist))
    cases_run = unique(c(cases_T1D$Run,cases_sero$Run))
    cases = metadata[match(cases_run,metadata$Run),]
    cases$age_condition = apply(cases, 1, function(myrow) min(as.numeric(myrow["age_t1d"]),as.numeric(myrow["age_mult_persist"]),na.rm = TRUE))
  }

  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR4_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR4*030X/0302")
  } else if (HLA_status == "DR4_DR8_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR8*0401/0402")
  } else if (HLA_status == "DR3_DR3_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR3*0501/0201*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR1_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR1*0101/0501")
  } else if (HLA_status == "DR4_DR13") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR13*0102/0604")
  }
  
  cases$diff_in_num = cases$age_condition - cases$age_at_collection
  cases =cases[cases$diff_in_num>0,]
  
  
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      cases_close_to_condition = cases[cases$diff_in_num <= times[1],]
      cases_far_from_condition = cases[cases$diff_in_num > times[1],]
      num_days_analyzing = times[1]
    } else if (length(times) == 2) {
      cases_close_to_condition = cases[cases$diff_in_num <= times[2] & cases$diff_in_num >= times[1],]
      cases_far_from_condition = cases[cases$diff_in_num > times[2],]
      num_days_analyzing = times[2] - times[1]
    }
  }
  if(nrow(cases_close_to_condition)>0) {
    cases_close_to_condition$getsCondition = 1
  }
  if(nrow(cases_far_from_condition)>0) {
    cases_far_from_condition$getsCondition = 0
  }


  
  cases_far_from_condition = calculate_autoantibody_number(cases_far_from_condition)
  cases_close_to_condition = calculate_autoantibody_number(cases_close_to_condition)

  # combine controls and cases temporarily so I can get to subject level
  data = rbind(cases_far_from_condition,cases_close_to_condition)

  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 3) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  
  total_subject_number = length(unique(data$SubjectID))
  
  num_train_subjects = ceiling(total_subject_number*proportion_train)
  set.seed(3456)
  
  train_subjects = sample(unique(data$SubjectID),size=num_train_subjects,replace=FALSE)
  test_subjects = setdiff(unique(data$SubjectID),train_subjects)
  if(length(train_subjects) < 3 | length(test_subjects) < 3) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }

  train_metadata = data[data$SubjectID%in%train_subjects,]
  
  #train_subjects_split <- split(train_subjects,   ceiling(seq_along(train_subjects) / (num_train_subjects/3)))
  #train_subjects_1 = train_subjects_split[[1]]
  #train_subjects_2 = train_subjects_split[[2]]
  #train_subjects_3 = train_subjects_split[[3]]
  
  test_metadata = data[data$SubjectID%in%test_subjects,]
  
  folds <- createFolds(as.factor(train_metadata$getsCondition), k=3, 
                       list = TRUE)
  
  train_subjects_1 = train_metadata[folds[[1]],c("Run","getsCondition")]
  train_subjects_2 = train_metadata[folds[[2]],c("Run","getsCondition")]
  train_subjects_3 = train_metadata[folds[[3]],c("Run","getsCondition")]
  test_subjects = test_metadata[,c("Run","getsCondition")]
  train_subjects = train_metadata[,c("Run","getsCondition")]
  
  
  # also scale grs2 score
  grs2_df = train_metadata[,"grs2",drop=FALSE]
  preProcValues_grs <- preProcess(grs2_df, method = c("center","scale"))
  trainTransformed <- predict(preProcValues_grs, train_metadata[,"grs2",drop=FALSE])
  testTransformed <- predict(preProcValues_grs, test_metadata[,"grs2",drop=FALSE])
  train_metadata$grs2 = trainTransformed$grs2
  test_metadata$grs2 = testTransformed$grs2
  
  if(antibody_type == "MIAA") {
    main_cat_prefix = "healthy_pre-MIAA"
  } else if (antibody_type == "GAD") {
    main_cat_prefix = "healthy_pre-GAD"
  } else if (antibody_type == "IA2A") {
    main_cat_prefix = "healthy_pre-IA2A"
  } else if (antibody_type == "seroconverters") {
    main_cat_prefix = "healthy_pre-sero"
  } else if (antibody_type == "triple_converters_vs_T1D") {
    main_cat_prefix = "triple_converters_vs_T1D"
  } else if (antibody_type == "T1D") {
    main_cat_prefix = "healthy_pre-t1d"
  } else if (antibody_type == "serconverters_or_T1D") {
    main_cat_prefix = "serconverters_or_T1D"
  }
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    time_out = paste(time_out,"_before_vs_far_from_condition",sep="")
    write.csv(test_metadata,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(train_metadata,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  } else {
    write.csv(test_metadata,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(train_metadata,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"willRun"))
}

times_df = data.frame(c(30,""),c(60,""),c(90,""),c(30,60),c(60,90))
times_df = t(times_df)
times_df = as.data.frame(times_df)

antibody_statuses = c("MIAA","GAD","IA2A","seroconverters","T1D","serconverters_or_T1D")

proportion_train_list = c(0.66,0.50)
HLA_statuses = c("all_HLA","DR3_DR4_only","DR4_DR4_only","DR4_DR8_only","DR3_DR3_only","DR4_DR1_only","DR4_DR13")


train_test_nums_pre_diff_before_after_conditions = matrix("",ncol=7,nrow=length(antibody_statuses)*length(HLA_statuses)*nrow(times_df)*length(proportion_train_list))
counter = 1
for (y in 1:length(antibody_statuses)) {
  antibody_status_temp = antibody_statuses[y]
  for(k in 1:length(HLA_statuses)) {
    HLA_status = HLA_statuses[k]
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
      for(j in 1:length(proportion_train_list)) {
        proportion_train = proportion_train_list[j]
        print(y)
        print(k)
        print(x)
        print(j)
        case_ctrl_samp_counts = get_metadata_healthy_days_before_condition_and_earlier(metadata=metadata,filter_by_time=filter_time,times=my_times,antibody_type=antibody_status_temp,HLA_status=HLA_status,proportion_train = proportion_train)
        train_test_nums_pre_diff_before_after_conditions[counter,] = case_ctrl_samp_counts
        counter = counter + 1
      }
    }
  }
}













get_metadata_healthy_days_before_condition_and_earlier_subect_level = function(metadata,filter_by_time,times,antibody_type,HLA_status,proportion_train) {
  if(antibody_type == "MIAA") {
    cases = metadata %>% filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA))
    cases$age_condition = cases$age_first_MIAA
  } else if (antibody_type == "GAD") {
    cases = metadata %>% filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD))
    cases$age_condition = cases$age_first_GAD
  } else if (antibody_type == "IA2A"){
    cases = metadata %>% filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A))
    cases$age_condition = cases$age_first_IA2A
  } else if (antibody_type == "seroconverters"){
    cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist))
    cases$age_condition = cases$age_mult_persist
  } else if (antibody_type == "triple_converters_vs_T1D") {
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before')
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "T1D") {
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before')
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "serconverters_or_T1D") {
    cases_T1D = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') 
    cases_sero = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist))
    cases_run = unique(c(cases_T1D$Run,cases_sero$Run))
    cases = metadata[match(cases_run,metadata$Run),]
    cases$age_condition = apply(cases, 1, function(myrow) min(as.numeric(myrow["age_t1d"]),as.numeric(myrow["age_mult_persist"]),na.rm = TRUE))
  }
  
  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR4_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR4*030X/0302")
  } else if (HLA_status == "DR4_DR8_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR8*0401/0402")
  } else if (HLA_status == "DR3_DR3_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR3*0501/0201*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR1_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR1*0101/0501")
  } else if (HLA_status == "DR4_DR13") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR13*0102/0604")
  }
  
  cases$diff_in_num = cases$age_condition - cases$age_at_collection
  cases =cases[cases$diff_in_num>0,]
  
  
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      cases_close_to_condition = cases[cases$diff_in_num <= times[1],]
      cases_far_from_condition = cases[cases$diff_in_num > times[1],]
      num_days_analyzing = times[1]
    } else if (length(times) == 2) {
      cases_close_to_condition = cases[cases$diff_in_num <= times[2] & cases$diff_in_num >= times[1],]
      cases_far_from_condition = cases[cases$diff_in_num > times[2],]
      num_days_analyzing = times[2] - times[1]
    }
  }
  if(nrow(cases_close_to_condition)>0) {
    cases_close_to_condition$getsCondition = 1
  }
  if(nrow(cases_far_from_condition)>0) {
    cases_far_from_condition$getsCondition = 0
  }
  
  
  
  cases_far_from_condition = calculate_autoantibody_number(cases_far_from_condition)
  cases_close_to_condition = calculate_autoantibody_number(cases_close_to_condition)
  
  data = rbind(cases_far_from_condition,cases_close_to_condition)
  data$SubjectID = paste(data$SubjectID,data$getsCondition,sep="_")
  mapping_file = data %>% dplyr::select(c('SubjectID','Run'))
  # just doing this to fit the column names of later script 
  colnames(mapping_file) <- c("V1","V2")
  
  # for each subject ID calculate the average age_at_collection
  avg_times = as.data.table(data)[,mean(age_at_collection),by="SubjectID"]
  data$time = avg_times[match(data$SubjectID,avg_times$SubjectID),]$V1
  
  data = collapse_on_subject_id(data)
  data$getsCondition = as.factor(data$getsCondition)
  
  # lets make sure we have paired data. for every case we have a control
  
  data$SubjectID_old = sapply(strsplit(data$SubjectID,split="_"), function(x) x[[1]])
  num_subjects_with_both_conditions = table(data$SubjectID_old)
  num_subjects_with_both_conditions = names(num_subjects_with_both_conditions[num_subjects_with_both_conditions==2])
  
  data = data[data$SubjectID_old%in%num_subjects_with_both_conditions,]
  
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 2) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  # we want to partition by subject
  
  total_subjects = unique(data$SubjectID_old)
  train_subjects = sample(total_subjects,ceiling(length(total_subjects)*proportion_train))
  trainIndex = which(data$SubjectID_old%in%train_subjects)
  train_df = data[trainIndex,]
  train_subjects <- data[ trainIndex,]$SubjectID
  test_df  <- data[-trainIndex,]
  
  # now lets split by subject again
  
  total_train_subjects = unique(train_df$SubjectID_old)
  if(length(total_train_subjects) < 3) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  subject_training_folds = split(total_train_subjects, ceiling(seq_along(total_train_subjects) / (length(total_train_subjects)/3)))
  folds = lapply(subject_training_folds, function(fold) {
    which(train_df$SubjectID_old%in%fold)
  })
  
  
  train_subjects_1 = train_df[folds[[1]],c("SubjectID","getsCondition")]
  train_subjects_2 = train_df[folds[[2]],c("SubjectID","getsCondition")]
  train_subjects_3 = train_df[folds[[3]],c("SubjectID","getsCondition")]
  test_subjects = test_df[,c("SubjectID","getsCondition")]
  
  # also scale grs2 score
  grs2_df = train_df[,"grs2",drop=FALSE]
  preProcValues_grs <- preProcess(grs2_df, method = c("center","scale"))
  trainTransformed <- predict(preProcValues_grs, train_df[,"grs2",drop=FALSE])
  testTransformed <- predict(preProcValues_grs, test_df[,"grs2",drop=FALSE])
  train_df$grs2 = trainTransformed$grs2
  test_df$grs2 = testTransformed$grs2
  
  ctrl_num_train = sum(train_df$getsCondition == 0)
  case_num_train = sum(train_df$getsCondition == 1)
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num_train < 3 | case_num_train < 1) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  if(antibody_type == "MIAA") {
    main_cat_prefix = "healthy_pre-MIAA"
  } else if (antibody_type == "GAD") {
    main_cat_prefix = "healthy_pre-GAD"
  } else if (antibody_type == "IA2A") {
    main_cat_prefix = "healthy_pre-IA2A"
  } else if (antibody_type == "seroconverters") {
    main_cat_prefix = "healthy_pre-sero"
  } else if (antibody_type == "triple_converters_vs_T1D") {
    main_cat_prefix = "triple_converters_vs_T1D"
  } else if (antibody_type == "T1D") {
    main_cat_prefix = "healthy_pre-t1d"
  } else if (antibody_type == "serconverters_or_T1D") {
    main_cat_prefix = "serconverters_or_T1D"
  }
  
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    time_out = paste(time_out,"_before_vs_far_from_condition_subject_level",sep="")
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  } else {
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"willRun"))
}

times_df = data.frame(c(30,""),c(60,""),c(90,""),c(30,60),c(60,90))
times_df = t(times_df)
times_df = as.data.frame(times_df)

antibody_statuses = c("MIAA","GAD","IA2A","seroconverters","T1D","serconverters_or_T1D")

proportion_train_list = c(0.66,0.50)
HLA_statuses = c("all_HLA","DR3_DR4_only","DR4_DR4_only","DR4_DR8_only","DR3_DR3_only","DR4_DR1_only","DR4_DR13")


train_test_nums_pre_diff_before_after_conditions_subject_level = matrix("",ncol=7,nrow=length(antibody_statuses)*length(HLA_statuses)*nrow(times_df)*length(proportion_train_list))
counter = 1
for (y in 1:length(antibody_statuses)) {
  antibody_status_temp = antibody_statuses[y]
  for(k in 1:length(HLA_statuses)) {
    HLA_status = HLA_statuses[k]
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
      for(j in 1:length(proportion_train_list)) {
        proportion_train = proportion_train_list[j]
        print(y)
        print(k)
        print(x)
        print(j)
        case_ctrl_samp_counts = get_metadata_healthy_days_before_condition_and_earlier_subect_level(metadata=metadata,filter_by_time=filter_time,times=my_times,antibody_type=antibody_status_temp,HLA_status=HLA_status,proportion_train = proportion_train)
        train_test_nums_pre_diff_before_after_conditions_subject_level[counter,] = case_ctrl_samp_counts
        counter = counter + 1
      }
    }
  }
}


get_metadata_healthy_days_before_condition_version2 = function(metadata,filter_by_time,times,antibody_type,HLA_status,proportion_train) {
  if(antibody_type == "MIAA") {
    controls = metadata %>% filter(is.na(age_first_MIAA)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_first_MIAA
    cases$age_condition = cases$age_first_MIAA
  } else if (antibody_type == "GAD") {
    controls = metadata %>% filter(is.na(age_first_GAD)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_first_GAD
    cases$age_condition = cases$age_first_GAD
  } else if (antibody_type == "IA2A"){
    controls = metadata %>% filter(is.na(age_first_IA2A)) %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_first_IA2A
    cases$age_condition = cases$age_first_IA2A
  } else if (antibody_type == "seroconverters"){
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
    controls$age_condition = controls$age_mult_persist
    cases$age_condition = cases$age_mult_persist
  } else if (antibody_type == "triple_converters_vs_T1D") {
    controls = metadata %>% filter(t1d_sero_control == 'seroconverted') %>% filter(three_persist_conf == TRUE) %>% mutate(condition = 0)
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    controls$age_condition = max(controls$age_first_MIAA,controls$age_first_IA2A,cases$age_first_GAD)
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "T1D") {
    controls = metadata %>% filter(t1d == FALSE) %>% mutate(condition = 0)
    cases = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    controls$age_condition = controls$age_t1d
    cases$age_condition = cases$age_t1d
  } else if (antibody_type == "serconverters_or_T1D") {
    controls = metadata %>% filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases_T1D = metadata %>% filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    cases_sero = metadata %>% filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
    cases_run = unique(c(cases_T1D$Run,cases_sero$Run))
    cases = metadata[match(cases_run,metadata$Run),]  %>% mutate(condition = 1)
    controls$age_condition = NA
    cases$age_condition = apply(cases, 1, function(myrow) min(as.numeric(myrow["age_t1d"]),as.numeric(myrow["age_mult_persist"]),na.rm = TRUE))
  }
  
  cases$diff_in_num = cases$age_condition - cases$age_at_collection
  controls$diff_in_num = controls$age_condition - controls$age_at_collection
  cases =cases[cases$diff_in_num>0,]
  
  
  if(filter_by_time == "yes") {
    if(length(times) == 1) {
      cases = cases[cases$diff_in_num <= times[1],]
    } else if (length(times) == 2) {
      cases = cases[cases$diff_in_num <= times[2] & cases$diff_in_num >= times[1],]
    }
  }
  
  if(HLA_status == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR4_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR4*030X/0302")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR8_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR8*0401/0402")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR3_DR3_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR3*0501/0201*DR3*0501/0201")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR1_only") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR1*0101/0501")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA_status == "DR4_DR13") {
    cases = cases %>% filter(HLA_Category!= "Not*Eligible")
    cases = cases %>% filter(HLA_Category== "DR4*030X/0302*DR13*0102/0604")
    controls = controls %>% filter(HLA_Category!= "Not*Eligible")
    controls = controls %>% filter(HLA_Category== "DR4*030X/0302*DR3*0501/0201")
  }
  
  
  controls = calculate_autoantibody_number(controls)
  cases = calculate_autoantibody_number(cases)
  # get average age at collection for each subject
  #time_controls = as.data.frame(as.data.table(controls)[,mean(age_at_collection),by=SubjectID])
  #controls$time = time_controls[match(controls$SubjectID,time_controls$SubjectID),"V1"]
  #time_cases = as.data.frame(as.data.table(cases)[,mean(age_at_collection),by=SubjectID])
  #cases$time = time_cases[match(cases$SubjectID,time_cases$SubjectID),"V1"]
  
  # now lets get an aged matched control for each case
  
  data = rbind(controls,cases)
  library(MatchIt)
  
  m.out1 <- matchit(condition ~ age_at_collection, data = data,
                    method = "optimal", distance = "glm")
  
  pairs = m.out1$match.matrix
  
  controls = data[as.numeric(pairs[,1]),]

  # get average age at collection for each subject
  time_controls = as.data.frame(as.data.table(controls)[,mean(age_at_collection),by=SubjectID])
  controls$time = time_controls[match(controls$SubjectID,time_controls$SubjectID),"V1"]
  time_cases = as.data.frame(as.data.table(cases)[,mean(age_at_collection),by=SubjectID])
  cases$time = time_cases[match(cases$SubjectID,time_cases$SubjectID),"V1"]
    
  data = rbind(controls,cases)
  
  
  mapping_file = rbind(cases,controls) %>% select(c('SubjectID','Run'))
  
  
  colnames(mapping_file) <- c("V1","V2")
  data = collapse_on_subject_id(data)

  colnames(data)[match("condition",colnames(data))] = "getsCondition"
  data$getsCondition = as.factor(data$getsCondition)
  
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num < 3 | case_num < 2) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  
  set.seed(3456)
  trainIndex <- createDataPartition(data$getsCondition, p = proportion_train, 
                                    list = FALSE, 
                                    times = 1)
  train_df = data[trainIndex,]
  train_subjects <- data[ trainIndex,]$SubjectID
  test_df  <- data[-trainIndex,]
  folds <- createFolds(train_df$getsCondition, k=3, 
                       list = TRUE)
  
  train_subjects_1 = train_df[folds[[1]],c("SubjectID","getsCondition")]
  train_subjects_2 = train_df[folds[[2]],c("SubjectID","getsCondition")]
  train_subjects_3 = train_df[folds[[3]],c("SubjectID","getsCondition")]
  test_subjects = test_df[,c("SubjectID","getsCondition")]
  
  # also scale grs2 score
  grs2_df = train_df[,"grs2",drop=FALSE]
  preProcValues_grs <- preProcess(grs2_df, method = c("center","scale"))
  trainTransformed <- predict(preProcValues_grs, train_df[,"grs2",drop=FALSE])
  testTransformed <- predict(preProcValues_grs, test_df[,"grs2",drop=FALSE])
  train_df$grs2 = trainTransformed$grs2
  test_df$grs2 = testTransformed$grs2
  
  ctrl_num_train = sum(train_df$getsCondition == 0)
  case_num_train = sum(train_df$getsCondition == 1)
  ctrl_num = sum(data$getsCondition == 0)
  case_num = sum(data$getsCondition == 1)
  
  # if we have too few control or case samples return number of train and test samples and don't continue
  if(ctrl_num_train < 3 | case_num_train < 1) {
    return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"wont_run"))
  }
  
  
  if(antibody_type == "MIAA") {
    main_cat_prefix = "healthy_pre-MIAA"
  } else if (antibody_type == "GAD") {
    main_cat_prefix = "healthy_pre-GAD"
  } else if (antibody_type == "IA2A") {
    main_cat_prefix = "healthy_pre-IA2A"
  } else if (antibody_type == "seroconverters") {
    main_cat_prefix = "healthy_pre-sero"
  } else if (antibody_type == "triple_converters_vs_T1D") {
    main_cat_prefix = "triple_converters_vs_T1D"
  } else if (antibody_type == "serconverters_or_T1D") {
    main_cat_prefix = "serconverters_or_T1D"
  } else if (antibody_type == "T1D") {
    main_cat_prefix = "healthy_pre-t1d"
  }
  
  if(filter_by_time == "yes") {
    times_in_months = round(times/30)
    time_out = paste(sapply(times_in_months,function(x) paste(x,"month",sep="")),collapse="-")
    time_out = paste(time_out,"_before_condition_case_vs_control",sep="")
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',time_out,"-",HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,"-",time_out,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  } else {
    write.csv(train_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_train_metadata.csv",sep=""))
    write.csv(test_df,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,"_test_metadata.csv",sep=""))
    write.csv(mapping_file,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,".mapping.csv",sep=""),row.names=FALSE)
    write.table(train_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_1,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(train_subjects_2,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(train_subjects_3,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".train_subjects_3.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(test_subjects,paste('~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/processed_teddy_metadata_for_regression_kfolds_v2/',main_cat_prefix,'-',HLA_status,"_proportion_train_",proportion_train,".test_subjects.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  return(c(antibody_type,HLA_status,paste(times,collapse="-"),proportion_train,case_num,ctrl_num,"willRun"))
}

times_df = data.frame(c(30,""),c(60,""),c(90,""),c(30,60),c(60,90))
times_df = t(times_df)
times_df = as.data.frame(times_df)

antibody_statuses = c("MIAA","GAD","IA2A","seroconverters","T1D")

proportion_train_list = c(0.66,0.50)
#HLA_statuses = c("all_HLA","DR3_DR4_only","DR4_DR4_only","DR4_DR8_only","DR3_DR3_only","DR4_DR1_only","DR4_DR13")
HLA_statuses = c("all_HLA")

#train_test_nums_pre_diff_before_after_conditions_subject_level_case_vs_control = matrix("",ncol=7,nrow=length(antibody_statuses)*length(HLA_statuses)*nrow(times_df)*length(proportion_train_list))
#counter = 1
for (y in 1:length(antibody_statuses)) {
  antibody_status_temp = antibody_statuses[y]
  for(k in 1:length(HLA_statuses)) {
    HLA_status = HLA_statuses[k]
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
      for(j in 1:length(proportion_train_list)) {
        proportion_train = proportion_train_list[j]
        print(y)
        print(k)
        print(x)
        print(j)
        case_ctrl_samp_counts = get_metadata_healthy_days_before_condition_version2(metadata=metadata,filter_by_time=filter_time,times=my_times,antibody_type=antibody_status_temp,HLA_status=HLA_status,proportion_train = proportion_train)
        #train_test_nums_pre_diff_before_after_conditions_subject_level[counter,] = case_ctrl_samp_counts
        #counter = counter + 1
      }
    }
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
