# script to go through diabimmune and teddy data and comnbine them

setwd("~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/analysis/")

library(tidyverse)
library(plyr)
#all healthy vs post T1D

all_healthy_vs_postT1d_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_post-t1d_metadata_filtered.rds")
all_healthy_vs_postT1d_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_post-t1d_metadata_filtered.rds")

# remove maternal age from diabimmune. remove. not found in teddy. remove HLA risk. risk scores in teddy and 
# diabimmune do not seem equivalent. remove location cause there is no urban or rural category in teddy. there is 
# living on a farm with animals but thats not the same. also we don't have metadata for ZNT8A and ICA in teddy
# also remove gestational diabetes. There is info on this in TEDDY but they also include other types of diabetes
toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","gestational_diabetes")
all_healthy_vs_postT1d_diabimmune = all_healthy_vs_postT1d_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,IA2A_pos=IA2A,GAD_pos=GADA,MIAA_pos=IAA)
# create seroconverted_ever category in teddy data
all_healthy_vs_postT1d_teddy$seroconverted_ever = as.numeric(!is.na(all_healthy_vs_postT1d_teddy$age_mult_persist))
# create MIAA variable 
all_healthy_vs_postT1d_teddy$MIAA_pos = as.numeric(!is.na(all_healthy_vs_postT1d_teddy$age_first_MIAA))
# convert boolean to numeric
all_healthy_vs_postT1d_teddy$GAD_pos = as.numeric(all_healthy_vs_postT1d_teddy$GAD_pos)
all_healthy_vs_postT1d_teddy$IA2A_pos = as.numeric(all_healthy_vs_postT1d_teddy$IA2A_pos)
all_healthy_vs_postT1d_teddy$Geographical_Location[all_healthy_vs_postT1d_teddy$Geographical_Location=="Finland"] = "FIN"

# add column saying which df is which
all_healthy_vs_postT1d_teddy = all_healthy_vs_postT1d_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_postT1d_diabimmune = all_healthy_vs_postT1d_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_postT1d = rbind.fill(all_healthy_vs_postT1d_teddy,all_healthy_vs_postT1d_diabimmune)

saveRDS(all_healthy_vs_postT1d,"processed_diab_teddy_metadata_for_regression/healthy_post-t1d.rds")


#all healthy vs post seroconversion
all_healthy_vs_postSero_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_post-sero_metadata_filtered.rds")
all_healthy_vs_postSero_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_post-sero_metadata_filtered.rds")
# remove diabetes_at_sampling cause that is a sample measurement, not a subject one. 
toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes","diabetes_at_sampling")

all_healthy_vs_postSero_diabimmune = all_healthy_vs_postSero_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,IA2A_pos=IA2A,GAD_pos=GADA)

all_healthy_vs_postSero_teddy$t1d_ever = as.numeric(!is.na(all_healthy_vs_postSero_teddy$age_t1d))
all_healthy_vs_postSero_teddy$GAD_pos = as.numeric(all_healthy_vs_postSero_teddy$GAD_pos)
all_healthy_vs_postSero_teddy$IA2A_pos = as.numeric(all_healthy_vs_postSero_teddy$IA2A_pos)
all_healthy_vs_postSero_teddy$Geographical_Location[all_healthy_vs_postSero_teddy$Geographical_Location=="Finland"] = "FIN"

all_healthy_vs_postSero_teddy = all_healthy_vs_postSero_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_postSero_diabimmune = all_healthy_vs_postSero_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_postSero = rbind.fill(all_healthy_vs_postSero_teddy,all_healthy_vs_postSero_diabimmune)
saveRDS(all_healthy_vs_postSero,"processed_diab_teddy_metadata_for_regression/healthy_post-sero.rds")

#all healthy vs IAA
all_healthy_vs_IAA_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_miaa_metadata_filtered.rds")
all_healthy_vs_IAA_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_miaa_metadata_filtered.rds")

toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","gestational_diabetes")
all_healthy_vs_IAA_diabimmune = all_healthy_vs_IAA_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,IA2A_pos=IA2A,GAD_pos=GADA)

all_healthy_vs_IAA_teddy$t1d_ever = as.numeric(!is.na(all_healthy_vs_IAA_teddy$age_t1d))
all_healthy_vs_IAA_teddy$Geographical_Location = as.character(all_healthy_vs_IAA_teddy$Geographical_Location)
all_healthy_vs_IAA_teddy$Geographical_Location[all_healthy_vs_IAA_teddy$Geographical_Location=="Finland"] = "FIN"
all_healthy_vs_IAA_teddy$seroconverted_ever = as.numeric(all_healthy_vs_IAA_teddy$two_or_more_persistent)
all_healthy_vs_IAA_teddy$GAD_pos = as.numeric(all_healthy_vs_IAA_teddy$GAD_pos)
all_healthy_vs_IAA_teddy$IA2A_pos = as.numeric(all_healthy_vs_IAA_teddy$IA2A_pos)

all_healthy_vs_IAA_teddy = all_healthy_vs_IAA_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_IAA_diabimmune = all_healthy_vs_IAA_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_IAA = rbind.fill(all_healthy_vs_IAA_teddy,all_healthy_vs_IAA_diabimmune)
saveRDS(all_healthy_vs_IAA,"processed_diab_teddy_metadata_for_regression/healthy_miaa.rds")

#all healthy vs GAD
all_healthy_vs_GAD_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_gad_metadata_filtered.rds")
all_healthy_vs_GAD_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_gad_metadata_filtered.rds")
toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes")
all_healthy_vs_GAD_diabimmune = all_healthy_vs_GAD_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,IA2A_pos=IA2A)

all_healthy_vs_GAD_teddy$t1d_ever = as.numeric(!is.na(all_healthy_vs_GAD_teddy$age_t1d))
all_healthy_vs_GAD_teddy$Geographical_Location[all_healthy_vs_GAD_teddy$Geographical_Location=="Finland"] = "FIN"
all_healthy_vs_GAD_teddy$seroconverted_ever = as.numeric(all_healthy_vs_GAD_teddy$two_or_more_persistent)
all_healthy_vs_GAD_teddy$IA2A_pos = as.numeric(all_healthy_vs_GAD_teddy$IA2A_pos)

all_healthy_vs_GAD_teddy = all_healthy_vs_GAD_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_GAD_diabimmune = all_healthy_vs_GAD_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_GAD = rbind.fill(all_healthy_vs_GAD_teddy,all_healthy_vs_GAD_diabimmune)
saveRDS(all_healthy_vs_GAD,"processed_diab_teddy_metadata_for_regression/healthy_gad.rds")


#all healthy vs IA2A
all_healthy_vs_IA2A_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_IA2A_metadata_filtered.rds")
all_healthy_vs_IA2A_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_IA2A_metadata_filtered.rds")
toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes")
all_healthy_vs_IA2A_diabimmune = all_healthy_vs_IA2A_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA)

all_healthy_vs_IA2A_teddy$t1d_ever = as.numeric(!is.na(all_healthy_vs_IA2A_teddy$age_t1d))
all_healthy_vs_IA2A_teddy$Geographical_Location[all_healthy_vs_IA2A_teddy$Geographical_Location=="Finland"] = "FIN"
all_healthy_vs_IA2A_teddy$seroconverted_ever = as.numeric(all_healthy_vs_IA2A_teddy$two_or_more_persistent)
all_healthy_vs_IA2A_teddy$GAD_pos = as.numeric(all_healthy_vs_IA2A_teddy$GAD_pos)

all_healthy_vs_IA2A_teddy = all_healthy_vs_IA2A_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_IA2A_diabimmune = all_healthy_vs_IA2A_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_IA2A = rbind.fill(all_healthy_vs_IA2A_teddy,all_healthy_vs_IA2A_diabimmune)
saveRDS(all_healthy_vs_IA2A,"processed_diab_teddy_metadata_for_regression/healthy_IA2A.rds")

#all healthy vs any T1D
all_healthy_vs_anyT1D_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_any-T1D_metadata_filtered.rds")
all_healthy_vs_anyT1D_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_any-T1D_metadata_filtered.rds")
toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes")
all_healthy_vs_anyT1D_diabimmune = all_healthy_vs_anyT1D_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA,IA2A_pos=IA2A)

all_healthy_vs_anyT1D_teddy$Geographical_Location[all_healthy_vs_anyT1D_teddy$Geographical_Location=="Finland"] = "FIN"
all_healthy_vs_anyT1D_teddy$seroconverted_ever = as.numeric(all_healthy_vs_anyT1D_teddy$two_or_more_persistent)
all_healthy_vs_anyT1D_teddy$IA2A_pos = as.numeric(all_healthy_vs_anyT1D_teddy$IA2A_pos)
all_healthy_vs_anyT1D_teddy$GAD_pos = as.numeric(all_healthy_vs_anyT1D_teddy$GAD_pos)

all_healthy_vs_anyT1D_teddy = all_healthy_vs_anyT1D_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_anyT1D_diabimmune = all_healthy_vs_anyT1D_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_anyT1D = rbind.fill(all_healthy_vs_anyT1D_teddy,all_healthy_vs_anyT1D_diabimmune)
saveRDS(all_healthy_vs_anyT1D,"processed_diab_teddy_metadata_for_regression/healthy_any-T1D.rds")

#all healthy vs any seroconversion
all_healthy_vs_anySero_diabimmune = readRDS("processed_diabimmune_metadata_for_regression/healthy_any-sero_metadata_filtered.rds")
all_healthy_vs_anySero_teddy = readRDS("processed_teddy_metadata_for_regression/healthy_any-sero_metadata_filtered.rds")

toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes")
all_healthy_vs_anySero_diabimmune = all_healthy_vs_anySero_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA,IA2A_pos=IA2A)

all_healthy_vs_anySero_teddy$Geographical_Location[all_healthy_vs_anySero_teddy$Geographical_Location=="Finland"] = "FIN"
all_healthy_vs_anySero_teddy$t1d_ever = as.numeric(!is.na(all_healthy_vs_anySero_teddy$age_t1d))
all_healthy_vs_anySero_teddy$IA2A_pos = as.numeric(all_healthy_vs_anySero_teddy$IA2A_pos)
all_healthy_vs_anySero_teddy$GAD_pos = as.numeric(all_healthy_vs_anySero_teddy$GAD_pos)

all_healthy_vs_anySero_teddy = all_healthy_vs_anySero_teddy %>% mutate(dataset = "TEDDY")
all_healthy_vs_anySero_diabimmune = all_healthy_vs_anySero_diabimmune %>% mutate(dataset = "diabimmune")

all_healthy_vs_anySero = rbind.fill(all_healthy_vs_anySero_teddy,all_healthy_vs_anySero_diabimmune)
saveRDS(all_healthy_vs_anySero,"processed_diab_teddy_metadata_for_regression/healthy_any-sero.rds")

# now do pre T1D overall and at different times

health_preT1D_file_list = c("healthy_pre-t1d_metadata_filtered.rds","healthy_pre-t1d-6month_metadata_filtered.rds","healthy_pre-t1d-6month-12month_metadata_filtered.rds","healthy_pre-t1d-12month_metadata_filtered.rds","healthy_pre-t1d-12-18month_metadata_filtered.rds","healthy_pre-t1d-18month_metadata_filtered.rds","healthy_pre-t1d-18-24month_metadata_filtered.rds","healthy_pre-t1d-24month_metadata_filtered.rds")

for(x in health_preT1D_file_list) {
  print(x)
  all_healthy_vs_preT1D_diabimmune = readRDS(paste("processed_diabimmune_metadata_for_regression/",x,sep=""))
  all_healthy_vs_preT1D_teddy = readRDS(paste("processed_teddy_metadata_for_regression/",x,sep=""))
  # seroconverted at sampling is not at subject level so I think we should remove this, right?
  if("seroconverted_at_sampling"%in%colnames(all_healthy_vs_preT1D_diabimmune)) {
    toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes","seroconverted_at_sampling")
  } else {
    toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","ZNT8A","ICA","IAA","gestational_diabetes")
  }
  # I'm guessing this is necessary cause IA2A is the same value for every subject in diabimmune data
  if("IA2A"%in%colnames(all_healthy_vs_preT1D_diabimmune)) {
    all_healthy_vs_preT1D_diabimmune = all_healthy_vs_preT1D_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA,IA2A_pos=IA2A)
  } else {
    all_healthy_vs_preT1D_diabimmune = all_healthy_vs_preT1D_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA)
  }
  all_healthy_vs_preT1D_teddy$Geographical_Location = as.character(all_healthy_vs_preT1D_teddy$Geographical_Location)
  all_healthy_vs_preT1D_teddy$Geographical_Location[all_healthy_vs_preT1D_teddy$Geographical_Location=="Finland"] = "FIN"
  all_healthy_vs_preT1D_teddy$seroconverted_ever = as.numeric(all_healthy_vs_preT1D_teddy$two_or_more_persistent)
  if("IA2A"%in%colnames(all_healthy_vs_preT1D_diabimmune)) {
    all_healthy_vs_preT1D_teddy$IA2A_pos = as.numeric(all_healthy_vs_preT1D_teddy$IA2A_pos)
  }
  all_healthy_vs_preT1D_teddy$GAD_pos = as.numeric(all_healthy_vs_preT1D_teddy$GAD_pos)
  
  all_healthy_vs_preT1D_teddy = all_healthy_vs_preT1D_teddy %>% mutate(dataset = "TEDDY")
  all_healthy_vs_preT1D_diabimmune = all_healthy_vs_preT1D_diabimmune %>% mutate(dataset = "diabimmune")
  
  all_healthy_vs_preT1D = rbind.fill(all_healthy_vs_preT1D_teddy,all_healthy_vs_preT1D_diabimmune)
  saveRDS(all_healthy_vs_preT1D,paste("processed_diab_teddy_metadata_for_regression/",x,sep=""))
}

# now do pre seroconversion overall and at different times

health_preSero_file_list = c("healthy_pre-sero_metadata_filtered.rds","healthy_pre-sero-6month_metadata_filtered.rds","healthy_pre-sero-6month-12month_metadata_filtered.rds","healthy_pre-sero-12month_metadata_filtered.rds","healthy_pre-sero-12-18month_metadata_filtered.rds","healthy_pre-sero-18month_metadata_filtered.rds","healthy_pre-sero-18-24month_metadata_filtered.rds","healthy_pre-sero-24month_metadata_filtered.rds")

for(x in health_preSero_file_list) {
  print(x)
  all_healthy_vs_preSero_diabimmune = readRDS(paste("processed_diabimmune_metadata_for_regression/",x,sep=""))
  all_healthy_vs_preSero_teddy = readRDS(paste("processed_teddy_metadata_for_regression/",x,sep=""))
  toremove_diabimmune = c("mom_age_at_birth","HLA_risk_class","location","IAA","gestational_diabetes")
  if("ZNT8A"%in%colnames(all_healthy_vs_preSero_diabimmune)) {
    toremove_diabimmune = c(toremove_diabimmune,"ZNT8A")
  }
  if("ICA"%in%colnames(all_healthy_vs_preSero_diabimmune)) {
    toremove_diabimmune = c(toremove_diabimmune,"ICA")
  }
  # remove cause this is sample level data?
  if("diabetes_at_sampling"%in%colnames(all_healthy_vs_preSero_diabimmune)) {
    toremove_diabimmune = c(toremove_diabimmune,"diabetes_at_sampling")
  }
  if("IA2A"%in%colnames(all_healthy_vs_preSero_diabimmune)) {
    all_healthy_vs_preSero_diabimmune = all_healthy_vs_preSero_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA,IA2A_pos=IA2A)
  } else {
    all_healthy_vs_preSero_diabimmune = all_healthy_vs_preSero_diabimmune %>% select(-c(all_of(toremove_diabimmune))) %>% dplyr::rename(Geographical_Location = country,Sex=gender,GAD_pos=GADA)
  }
  all_healthy_vs_preSero_teddy$Geographical_Location = as.character(all_healthy_vs_preSero_teddy$Geographical_Location)
  all_healthy_vs_preSero_teddy$Geographical_Location[all_healthy_vs_preSero_teddy$Geographical_Location=="Finland"] = "FIN"
  all_healthy_vs_preSero_teddy$t1d_ever = as.numeric(!is.na(all_healthy_vs_preSero_teddy$age_t1d))
  if("IA2A"%in%colnames(all_healthy_vs_preSero_diabimmune)) {
    all_healthy_vs_preSero_teddy$IA2A_pos = as.numeric(all_healthy_vs_preSero_teddy$IA2A_pos)
  }
  all_healthy_vs_preSero_teddy$GAD_pos = as.numeric(all_healthy_vs_preSero_teddy$GAD_pos)

  all_healthy_vs_preSero_teddy = all_healthy_vs_preSero_teddy %>% mutate(dataset = "TEDDY")
  all_healthy_vs_preSero_diabimmune = all_healthy_vs_preSero_diabimmune %>% mutate(dataset = "diabimmune")
  
  all_healthy_vs_preSero = rbind.fill(all_healthy_vs_preSero_teddy,all_healthy_vs_preSero_diabimmune)
  saveRDS(all_healthy_vs_preSero,paste("processed_diab_teddy_metadata_for_regression/",x,sep=""))
}

