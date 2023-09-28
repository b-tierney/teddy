library(caret)
library(data.table)
library(reticulate)
library(pROC)
library(timeROC)
library(glmnet)
library(doMC)
library(survival)
library(caTools)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
set.seed(123)
metadata_train = args[1]
metadata_test = args[2]
feature_list = args[3]

metadata_train = noquote(metadata_train)
metadata_test = noquote(metadata_test)

run_cox_regression = function(metadata_train,metadata_test,feature_list) {
  
  # get the begin time and end time
  if(grepl("24month",basename(metadata_test))) {
    begin_time = 730
  } else if(grepl("18month",basename(metadata_test))) {
    begin_time = 548
  } else if(grepl("12month",basename(metadata_test))) {
    begin_time = 365
  } else if(grepl("6month",basename(metadata_test))) {
    begin_time = 183
  } else if(grepl("3month",basename(metadata_test))) {
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
  train_case_counts = train_counts["1"]
  test_ctrl_counts = test_counts["0"]
  test_case_counts = test_counts["1"]

  sample_counts = c(train_ctrl=train_ctrl_counts,train_case=train_case_counts,total_train=train_ctrl_counts+train_case_counts,test_ctrl=test_ctrl_counts,test_case=test_case_counts,total_test=test_ctrl_counts+test_case_counts)
  names(sample_counts) = c("train_ctrl","train_case","total_train","test_ctrl","test_case","total_test")

  if("number_autoantibodies"%in%colnames(metadata_train_df)) {
    metadata_train_df$number_autoantibodies = as.factor(metadata_train_df$number_autoantibodies)
  }
  if("number_autoantibodies"%in%colnames(metadata_test_df)) {
    metadata_test_df$number_autoantibodies = as.factor(metadata_test_df$number_autoantibodies)
  }

  feature_vec = strsplit(feature_list,split=",")[[1]]
  feature_vec = intersect(feature_vec,colnames(metadata_train_df))

  # its possible there are features in the train set that are not in the test set
  
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
  

  metadata_train_df = metadata_train_df[,c("time","getsCondition",feature_vec)]
  metadata_test_df = metadata_test_df[,c("time","getsCondition",feature_vec)]

  my_formula <- paste0("~",paste(feature_vec, collapse = " + "))
  my_formula <- as.formula(paste("Surv(time,getsCondition)",my_formula,sep=""))

  train_cox_regression = function(my_formula,metadata_train_df,metadata_test_df,day_end,weighting) {
    if(weighting == TRUE) {
      weight_0 = 1-sum(metadata_train_df[,"getsCondition"]==0)/nrow(metadata_train_df)
      weight_1 = 1-sum(metadata_train_df[,"getsCondition"]==1)/nrow(metadata_train_df)
      weights = rep(0,nrow(metadata_train_df))
      weights[metadata_train_df[,"getsCondition"]==0] = weight_0
      weights[metadata_train_df[,"getsCondition"]==1] = weight_1
    } else {
      weights = rep(1,nrow(metadata_train_df))
    }
    model_train_min <- coxph(my_formula,data=metadata_train_df,control = coxph.control(iter.max = 50),weights=weights)
    predict_train = predict(model_train_min,newdata=metadata_train_df,type = "risk")
    predict_test = predict(model_train_min,newdata=metadata_test_df,type = "risk")
    train_coefs = model_train_min$coefficients

    # get AUCs
    ROC.T_train <- timeROC(T = metadata_train_df[,"time"],
        delta = metadata_train_df[,"getsCondition"],marker = predict_train,
        cause = 1, weighting = "marginal",
        times = day_end,
        iid = TRUE) # set to false for now
    ROC.T_test <- timeROC(T = metadata_test_df[,"time"],
          delta = metadata_test_df[,"getsCondition"],marker = predict_test,
          cause = 1, weighting = "marginal",
          times = day_end,
          iid = TRUE) # set to false for now

    #now get AUCs for precision recall
    calc_precision_recall = function(risks,all_data) {
      if(min(risks) >0) {
        risks2 = c(0,risks)
      } else {
        risks2 = risks
      }
      times_temp = all_data[,"time"]
      delta_temp = all_data[,"getsCondition"]
      train_precision_recall_vals = lapply(risks2,function(relative_risk) {
        train_stats_temp = SeSpPPVNPV(cutpoint=relative_risk,T=times_temp,delta = delta_temp,marker = risks,cause = 1, weighting = "marginal",times = day_end,iid = FALSE)
        precision_vals = train_stats_temp$PPV
        recall_vals = train_stats_temp$TP
        return(list(precision_vals,recall_vals))
      })
      precision_vals = do.call("rbind",lapply(train_precision_recall_vals,function(x) x[[1]]))
      recall_vals = do.call("rbind",lapply(train_precision_recall_vals,function(x) x[[2]]))
      precision_vals = as.data.frame(precision_vals)
      recall_vals = as.data.frame(recall_vals)
      precision_vals = precision_vals[order(risks2,decreasing=TRUE),]
      recall_vals = recall_vals[order(risks2,decreasing=TRUE),]
      # calculate trapezoid rule
      precision_recall_aucs = sapply(seq(1,ncol(recall_vals)), function(mycount) {
         recalls_temp = recall_vals[,mycount]
         precision_temp = precision_vals[,mycount]
         NA_vals_precision = which(is.na(precision_temp))
         NA_vals_recall = which(is.na(recalls_temp))
         NA_vals = c(NA_vals_precision,NA_vals_recall)
         if(length(NA_vals) > 0) {
            recalls_temp = recalls_temp[-NA_vals]
            precision_temp = precision_temp[-NA_vals]
         }
         if(length(precision_temp)>0) {
           precision_recall_auc = trapz(recalls_temp,precision_temp)
           return(precision_recall_auc)
         } else {
           return(NA)
        }
      })
      return(precision_recall_aucs)
    }

    pr_all_train = calc_precision_recall(predict_train,metadata_train_df)
    pr_all_test = calc_precision_recall(predict_test,metadata_test_df)

    pr_mat_train = matrix(NA,ncol=4,nrow=100)
    for (x in 1:100) {
      index_to_sample = sample(x = seq(1,length(predict_train)), size  = length(predict_train),replace=TRUE)
      sampled_risk = predict_train[index_to_sample]
      sampled_data = metadata_train_df[index_to_sample,,drop=FALSE]
      pr_temp = calc_precision_recall(sampled_risk,sampled_data)
      pr_mat_train[x,] = pr_temp
    }

    pr_mat_test = matrix(NA,ncol=4,nrow=100)
    for (x in 1:100) {
      index_to_sample = sample(x = seq(1,length(predict_test)), size  = length(predict_test),replace=TRUE)
      sampled_risk = predict_test[index_to_sample]
      sampled_data = metadata_test_df[index_to_sample,,drop=FALSE]
      pr_temp = calc_precision_recall(sampled_risk,sampled_data)
      pr_mat_test[x,] = pr_temp
    }
    # now calculate 2.5% and 97.5% quantiles
    train_pr_CI = apply(pr_mat_train, 2, function(mycol) {
      myquantiles = quantile(mycol,c(0.025,0.975),na.rm=TRUE)
      lowerCI = myquantiles[1]
      upperCI = myquantiles[2]
      CI = paste(lowerCI,upperCI,sep="-")
      return(CI)
    })

    test_pr_CI = apply(pr_mat_test, 2, function(mycol) {
      myquantiles = quantile(mycol,c(0.025,0.975),na.rm=TRUE)
      lowerCI = myquantiles[1]
      upperCI = myquantiles[2]
      CI = paste(lowerCI,upperCI,sep="-")
      return(CI)
    })

    train_pr_cis = cbind(pr_all_train,train_pr_CI)
    test_pr_cis = cbind(pr_all_test,test_pr_CI)
    rownames(train_pr_cis) = day_end
    rownames(test_pr_cis) = day_end

    all_pr_cis = cbind(train_pr_cis,test_pr_cis)

    final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec,all_pr_cis)
    return(final_results)
 }
    poss_train_cox_regression = possibly(.f = train_cox_regression, otherwise = NA)
    weighted_cox = poss_train_cox_regression(my_formula,metadata_train_df,metadata_test_df,day_end,weighting=TRUE)
    unweighted_cox = poss_train_cox_regression(my_formula,metadata_train_df,metadata_test_df,day_end,weighting=FALSE)
    return(list(weighted_cox,unweighted_cox))
}

model_results = run_cox_regression(metadata_train=metadata_train,metadata_test=metadata_test,feature_list=feature_list)

out_path = dirname(metadata_train)
out_prefix = basename(metadata_train)
out_prefix = gsub("_train_metadata.csv","",out_prefix)
out_path = paste(out_path,out_prefix,sep="/")

output_weighted = paste(out_path,'/',out_prefix,"output_cox_regression_time_to_event","_feature_list_",feature_list,"_weighted.rds",sep="")
output_unweighted = paste(out_path,'/',out_prefix,"output_cox_regression_time_to_event","_feature_list_",feature_list,".rds",sep="")
saveRDS(model_results[[1]],output_weighted)
saveRDS(model_results[[2]],output_unweighted)


