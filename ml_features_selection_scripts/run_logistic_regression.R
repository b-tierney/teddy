library(caret)
library(data.table)
library(reticulate)
library(pROC)
library(timeROC)
library(glmnet)
library(doMC)
library(survival)
library(purrr)
library(MLeval)
library(PRROC)

#use_condaenv("/home/sez10/miniconda3_2/envs/python_env")
#source_python('/n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/read_df.py')

args = commandArgs(trailingOnly=TRUE)
set.seed(123)

metadata_train = args[1]
metadata_test = args[2]
feature_list = args[3]

metadata_train = noquote(metadata_train)
metadata_test = noquote(metadata_test)

print(metadata_train)
print(metadata_test)

run_logistic_regression = function(metadata_train,metadata_test,feature_list) {

  
  # get the begin time and end time
  if(grepl("24month",basename(metadata_train))) {
    begin_time = 730
  } else if(grepl("18month",basename(metadata_train))) {
    begin_time = 548
  } else if(grepl("12month",basename(metadata_train))) {
    begin_time = 365
  } else if(grepl("6month",basename(metadata_train))) {
    begin_time = 183
  } else if(grepl("3month",basename(metadata_train))) {
    begin_time = 92
  } else {
    begin_time = 7300
  }
  day_endv <- c(365.25,365.25*3,5*365.25,8*365.25)
  if(begin_time == 7300) {
    day_end = 7300
  }
  day_end = begin_time + day_endv

  metadata_train_df = read.csv(metadata_train)
  metadata_test_df = read.csv(metadata_test)
  rownames(metadata_train_df) = metadata_train_df$SubjectID
  rownames(metadata_test_df) = metadata_test_df$SubjectID
  
  train_counts = table(metadata_train_df$getsCondition)
  test_counts = table(metadata_test_df$getsCondition)
  
  train_ctrl_counts = train_counts["0"]
  if(is.na(train_ctrl_counts)) {
    train_ctrl_counts = 0
  }
  train_case_counts = train_counts["1"]
  if(is.na(train_case_counts)) {
    train_case_counts = 0
  }
  test_ctrl_counts = test_counts["0"]
  if(is.na(test_ctrl_counts)) {
    test_ctrl_counts = 0
  }
  test_case_counts = test_counts["1"]
  if(is.na(test_case_counts)) {
    test_case_counts = 0
  }
  sample_counts = c(train_ctrl=train_ctrl_counts,train_case=train_case_counts,total_train=sum(c(train_ctrl_counts,train_case_counts),na.rm=TRUE),test_ctrl=test_ctrl_counts,test_case=test_case_counts,total_test=sum(c(test_ctrl_counts,test_case_counts),na.rm=TRUE))
  names(sample_counts) = c("train_ctrl","train_case","total_train","test_ctrl","test_case","total_test")

  num_unique_antibodies = length(unique(metadata_train_df$number_autoantibodies))
  if(num_unique_antibodies>1) {
    metadata_train_df$number_autoantibodies = as.factor(metadata_train_df$number_autoantibodies)
    antibody_dummies = model.matrix(~number_autoantibodies, metadata_train_df)
    antibody_dummies = antibody_dummies[,-1,drop=FALSE]
    #if(ncol(antibody_dummies)==1) {
    #  colnames(antibody_dummies) = paste(colnames(antibody_dummies),seq(1,ncol(antibody_dummies)),sep="")
    #}
    metadata_train_df = metadata_train_df[,-match("number_autoantibodies",colnames(metadata_train_df))]
    metadata_train_df = cbind(metadata_train_df,antibody_dummies)
  } else {
    # remove number_autoantibodies cause it has a variance of 0. not a useful predictor
    if("number_autoantibodies"%in%colnames(metadata_train_df)) {
      metadata_train_df = metadata_train_df[,-match("number_autoantibodies",colnames(metadata_train_df))]
    }
  }

  num_unique_antibodies_test = length(unique(metadata_test_df$number_autoantibodies))
  if(num_unique_antibodies_test>1) {
    metadata_test_df$number_autoantibodies = as.factor(metadata_test_df$number_autoantibodies)
    antibody_dummies = model.matrix(~number_autoantibodies, metadata_test_df)
    antibody_dummies = antibody_dummies[,-1,drop=FALSE]
    #if(ncol(antibody_dummies)==1) {
    #  colnames(antibody_dummies) = paste(colnames(antibody_dummies),seq(1,ncol(antibody_dummies)),sep="")
    #}
    metadata_test_df = metadata_test_df[,-match("number_autoantibodies",colnames(metadata_test_df))]
    metadata_test_df = cbind(metadata_test_df,antibody_dummies)
  } else {
    if("number_autoantibodies"%in%colnames(metadata_test_df)) {
      metadata_test_df = metadata_test_df[,-match("number_autoantibodies",colnames(metadata_test_df))]
    }
  }


  # its possible that the number of antibodies in train and test are different so we have to prepare for this
  missing_vars = setdiff(colnames(metadata_train_df),colnames(metadata_test_df))
  if(length(missing_vars) >0) {
    for(x in 1:length(missing_vars)) {
      newVar = rep(0,nrow(metadata_test_df))
      metadata_test_df = cbind(newVar,metadata_test_df)
      colnames(metadata_test_df)[1] = missing_vars[x]
    }
  }
  # now make sure test metadata has same columns as train metadata
  metadata_test_df = metadata_test_df[,match(colnames(metadata_train_df),colnames(metadata_test_df))]
  
  metadata_train_df$getsCondition = as.factor(metadata_train_df$getsCondition)
  metadata_test_df$getsCondition = as.factor(metadata_test_df$getsCondition)
  
  feature_vec = strsplit(feature_list,split=",")[[1]]
  if("number_autoantibodies"%in%feature_vec & num_unique_antibodies>1) {
    metadata_cols_to_take = colnames(metadata_train_df)[grep("number_autoantibodies",colnames(metadata_train_df))]
    feature_vec = feature_vec[-match("number_autoantibodies",feature_vec)]
    feature_vec = c(feature_vec,metadata_cols_to_take)
  } else if("number_autoantibodies"%in%feature_vec & num_unique_antibodies==1) {
    feature_vec = feature_vec[-match("number_autoantibodies",feature_vec)]
  }
  feature_vec = intersect(feature_vec,colnames(metadata_train_df))


  metadata_train_df = metadata_train_df[,c("time","getsCondition",feature_vec)]
  metadata_test_df = metadata_test_df[,c("time","getsCondition",feature_vec)]

  train_logistic_regression = function(train_ctrl_counts,train_case_counts,metadata_train_df,metadata_test_df,weighting) {
    if(weighting == TRUE) {
      weight_0 = 1-sum(metadata_train_df[,"getsCondition"]==0)/nrow(metadata_train_df)
      weight_1 = 1-sum(metadata_train_df[,"getsCondition"]==1)/nrow(metadata_train_df)
      weights = rep(0,nrow(metadata_train_df))
      weights[metadata_train_df[,"getsCondition"]==0] = weight_0
      weights[metadata_train_df[,"getsCondition"]==1] = weight_1
    } else {
      weights = rep(1,nrow(metadata_train_df))
    }
   
    if(train_ctrl_counts>1 & train_case_counts>1) {
      metadata_train_df$getsCondition = as.numeric(as.character(metadata_train_df$getsCondition))
      metadata_test_df$getsCondition = as.numeric(as.character(metadata_test_df$getsCondition))
      my_formula <- paste0("~",paste(feature_vec, collapse = " + "))
      my_formula <- as.formula(paste("getsCondition",my_formula,sep=""))
      model_train_out = glm(my_formula,family=binomial(link='logit'),data=metadata_train_df,weights=weights)
      predict_prob_train = predict(model_train_out,newdata=metadata_train_df,type = "response")
      predict_prob_train = as.data.frame(predict_prob_train)
      predict_prob_test = predict(model_train_out,newdata=metadata_test_df,type = "response")
      predict_prob_test = as.data.frame(predict_prob_test)
      predictions_train = predict_prob_train
      predictions_test = predict_prob_test
      predictions_train[predictions_train[,1]>0.5,1] = 1
      predictions_train[predictions_train[,1]<=0.5,1] = 0
      predictions_test[predictions_test[,1]>0.5,1] = 1
      predictions_test[predictions_test[,1]<=0.5,1] = 0
      train_coefs = model_train_out$coefficients
      predict_prob_train_with_labels = predict_prob_train
      predict_prob_test_with_labels = predict_prob_test
      colnames(predict_prob_train_with_labels) = c("case")
      colnames(predict_prob_test_with_labels) = c("case")
      predict_prob_train_with_labels = as.data.frame(predict_prob_train_with_labels)
      predict_prob_test_with_labels = as.data.frame(predict_prob_test_with_labels)
      predict_prob_train_with_labels$ctrl = 1-predict_prob_train_with_labels$case
      predict_prob_test_with_labels$ctrl = 1-predict_prob_test_with_labels$case
      predict_prob_test_with_labels$obs = as.character(metadata_test_df$getsCondition)
      predict_prob_train_with_labels$obs = as.character(metadata_train_df$getsCondition)
      predict_prob_train_with_labels$obs[predict_prob_train_with_labels$obs==1] = "case"
      predict_prob_train_with_labels$obs[predict_prob_train_with_labels$obs==0] =	"ctrl"
      predict_prob_test_with_labels$obs[predict_prob_test_with_labels$obs==1] = "case"
      predict_prob_test_with_labels$obs[predict_prob_test_with_labels$obs==0] = "ctrl"

      # get AUCs
      tryCatch({
        train_roc_out = evalm(predict_prob_train_with_labels,showplots=FALSE)
        train_roc_out_stats = train_roc_out$optres$Group1
        ROC.T_train = rep("",nrow(train_roc_out_stats)*2)
        names(ROC.T_train) = unlist(lapply(rownames(train_roc_out_stats), function(x) paste(x,colnames(train_roc_out_stats),sep="_")))
        counter = 1
        for(x in 1:nrow(train_roc_out_stats)) {
          myrow = train_roc_out_stats[x,]
          ROC.T_train[counter] = myrow[1,1]
          counter = counter + 1
          ROC.T_train[counter] = myrow[1,2]
          counter = counter + 1
        }
        # now do bootstrapping
        train_pr_all = pr.curve(scores.class0 = predict_prob_train[,1],weights.class0=as.numeric(as.character(metadata_train_df$getsCondition)))
        train_pr_all = train_pr_all$auc.integral
        train_pr_vec = rep(0,100)
        for(x in 1:100) {
          index_to_sample = sample(x = seq(1,nrow(predict_prob_train)), size  = nrow(predict_prob_train),replace=TRUE)	  
	  sampled_probs = predict_prob_train[index_to_sample,1]
	  sampled_labels = as.numeric(as.character(metadata_train_df$getsCondition))[index_to_sample]
	  pr_train <- pr.curve(scores.class0 = sampled_probs, weights.class0 = sampled_labels)
	  pr_train = pr_train$auc.integral
	  train_pr_vec[x] = pr_train
        }
        train_pr_CI = paste(quantile(train_pr_vec,c(0.025,0.975),na.rm=TRUE),collapse="-")
        ROC.T_train["AUC-PR_Score"] = train_pr_all
        ROC.T_train["AUC-PR_CI"] = train_pr_CI

        test_roc_out = evalm(predict_prob_test_with_labels,showplots=FALSE)
        test_roc_out_stats = test_roc_out$optres$Group1
        ROC.T_test = rep("",nrow(test_roc_out_stats)*2)
        names(ROC.T_test) = unlist(lapply(rownames(test_roc_out_stats), function(x) paste(x,colnames(test_roc_out_stats),sep="_")))
        counter = 1
        for(x in 1:nrow(test_roc_out_stats)) {
          myrow = test_roc_out_stats[x,]
          ROC.T_test[counter] = myrow[1,1]
          counter = counter + 1
          ROC.T_test[counter] = myrow[1,2]
          counter = counter + 1
        }
        # now do bootstrapping
        test_pr_all = pr.curve(scores.class0 = predict_prob_test[,1],weights.class0=as.numeric(as.character(metadata_test_df$getsCondition)))
        test_pr_all = test_pr_all$auc.integral
        test_pr_vec = rep(0,100)
        for(x in 1:100) {
          index_to_sample = sample(x = seq(1,nrow(predict_prob_test)), size  = nrow(predict_prob_test),replace=TRUE)
          sampled_probs = predict_prob_test[index_to_sample,1]
          sampled_labels = as.numeric(as.character(metadata_test_df$getsCondition))[index_to_sample]
          pr_test <- pr.curve(scores.class0 = sampled_probs, weights.class0 = sampled_labels)
          pr_test = pr_test$auc.integral
          test_pr_vec[x] = pr_test
        }
        test_pr_CI = paste(quantile(test_pr_vec,c(0.025,0.975),na.rm=TRUE),collapse="-")
        ROC.T_test["AUC-PR_Score"] = test_pr_all
        ROC.T_test["AUC-PR_CI"] = test_pr_CI
        final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec)
        return(final_results)	
      },error = function(e) {
       ROC.T_train = rep(NA,30)
       ROC.T_test = rep(NA,nrow(test_roc_out_stats)*2)
       names(ROC.T_test) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
       names(ROC.T_train) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
       final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec)
       return(final_results)
       })
    } else {
       ROC.T_train = rep(NA,30)
       ROC.T_test = rep(NA,30)
       names(ROC.T_train) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
       names(ROC.T_test) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
       train_coefs = NA
       final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec)
       return(final_results)
      }
    }
    poss_train_logistic_regression = possibly(.f = train_logistic_regression, otherwise = NA)
    weighted_results = poss_train_logistic_regression(train_ctrl_counts,train_case_counts,metadata_train_df,metadata_test_df,weighting=TRUE)  
    unweighted_results = poss_train_logistic_regression(train_ctrl_counts,train_case_counts,metadata_train_df,metadata_test_df,weighting=FALSE)

    return(list(weighted_results,unweighted_results))
  }

model_results = run_logistic_regression(metadata_train=metadata_train,metadata_test=metadata_test,feature_list=feature_list)

#out_path = dirname(train_abundance_data)
#out_prefix = basename(out_path)
#if(abundance_type == "species") {
#  out_prefix = paste(out_prefix,"species",sep="_")
#} else if(abundance_type == "pathway") {
#  out_prefix = paste(out_prefix,"pathway",sep="_")
#}

out_path = dirname(metadata_train)
out_prefix = basename(metadata_train)
out_prefix = gsub("_train_metadata.csv","",out_prefix)
out_path = paste(out_path,out_prefix,sep="/")


output_weighted = paste(out_path,'/',out_prefix,"output_plain_logistic_regression_feature_list_",feature_list,"_weighted.rds",sep="")
output_unweighted = paste(out_path,'/',out_prefix,"output_plain_logistic_regression_feature_list_",feature_list,".rds",sep="")


saveRDS(model_results[[1]],output_weighted)
saveRDS(model_results[[2]],output_unweighted)


