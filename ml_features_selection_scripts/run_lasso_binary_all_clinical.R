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

test_abundance_data = args[1]
train_abundance_data = args[2]
metadata_train = args[3]
metadata_test = args[4]
train_sujbects1 = args[5]
train_sujbects2 = args[6]
train_sujbects3 = args[7]
test_subjects = args[8]
feature_list = args[9]
features_to_keep = args[10]
microbiome_features = args[11]
number_threads = as.numeric(args[12])
loss_function = args[13]

test_abundance_data = noquote(test_abundance_data)
train_abundance_data = noquote(train_abundance_data)
metadata_train = noquote(metadata_train)
metadata_test = noquote(metadata_test)
train_sujbects1 = noquote(train_sujbects1)
train_sujbects2 = noquote(train_sujbects2)
train_sujbects3 = noquote(train_sujbects3)
test_subjects = noquote(test_subjects)

print(test_abundance_data)
print(train_abundance_data)
print(metadata_train)
print(metadata_test)
print(train_sujbects1)
print(train_sujbects2)
print(train_sujbects3)
print(test_subjects)
print(feature_list)

feature_vec_1 = strsplit(feature_list,split=",")[[1]]

if("microbiome"%in%feature_vec_1) {
  if(grepl("species",basename(test_abundance_data))) {
    abundance_type = "species"
  } else if(grepl("pathabundance",basename(test_abundance_data))) {
    abundance_type = "pathway"
  } else {
    abundance_type = "gene"
  }
} else {
  abundance_type = "none"
}

run_lasso = function(test_abundance_data,train_abundance_data,metadata_train,metadata_test,train_sujbects1,train_sujbects2,train_sujbects3,test_subjects,feature_list,features_to_keep,microbiome_features,abundance_type) {

  
  # get the begin time and end time
  if(grepl("24month",basename(test_abundance_data))) {
    begin_time = 730
  } else if(grepl("18month",basename(test_abundance_data))) {
    begin_time = 548
  } else if(grepl("12month",basename(test_abundance_data))) {
    begin_time = 365
  } else if(grepl("6month",basename(test_abundance_data))) {
    begin_time = 183
  } else if(grepl("3month",basename(test_abundance_data))) {
    begin_time = 92
  } else {
    begin_time = 7300
  }
  day_endv <- c(365.25,365.25*3,5*365.25,8*365.25)
  if(begin_time == 7300) {
    day_end = 7300
  }
  day_end = begin_time + day_endv

  feature_vec = strsplit(feature_list,split=",")[[1]]
  if("microbiome"%in%feature_vec) {
    
    if(abundance_type == "gene") {
      test_abundance_data_df = fread(test_abundance_data,data.table=FALSE)
      test_gene_names = read.table(paste(dirname(test_abundance_data),"/",basename(dirname(test_abundance_data)),"filtered_abundance_test_names.txt",sep=""),header=FALSE)[,1]
      rownames(test_abundance_data_df) = test_gene_names

      train_abundance_data_df = fread(train_abundance_data,data.table=FALSE)
      train_gene_names = read.table(paste(dirname(train_abundance_data),"/",basename(dirname(train_abundance_data)),"filtered_abundance_train_names.txt",sep=""),header=FALSE)[,1]
      rownames(train_abundance_data_df) = train_gene_names
  
    } else {
      test_abundance_data_df = fread(test_abundance_data,data.table=FALSE)
      train_abundance_data_df = fread(train_abundance_data,data.table=FALSE)
      test_names = read.table(gsub("_test.csv","_test_names.csv",test_abundance_data),header=FALSE)[,1]
      train_names = read.table(gsub("_train.csv","_train_names.csv",train_abundance_data),header=FALSE)[,1]
      rownames(test_abundance_data_df) = test_names
      rownames(train_abundance_data_df) = train_names
    }

    if(microbiome_features != "all" & microbiome_features == "ttest_sig") {
      if(abundance_type == "gene") {
        sig_genes_file = paste(dirname(test_abundance_data),"/",basename(dirname(test_abundance_data)),"output_ttest_results.rds",sep="")
      } else if(abundance_type == "species") {
        sig_genes_file = paste(dirname(test_abundance_data),"/",basename(dirname(test_abundance_data)),"_speciesoutput_ttest_results.rds",sep="")
      } else if(abundance_type == "pathway") {
        sig_genes_file = paste(dirname(test_abundance_data),"/",basename(dirname(test_abundance_data)),"_pathwayoutput_ttest_results.rds",sep="")
      }
      sig_genes_df = readRDS(sig_genes_file)
      sig_genes_df = names(sig_genes_df[[2]])
      train_abundance_data_df = train_abundance_data_df[,sig_genes_df,drop=FALSE]
      test_abundance_data_df = test_abundance_data_df[,sig_genes_df,drop=FALSE]
    }

    if(abundance_type != "gene") {
      colnames(train_abundance_data_df) = gsub("[|]","_",colnames(train_abundance_data_df))
      colnames(train_abundance_data_df) = gsub(" ","_",colnames(train_abundance_data_df))

      colnames(test_abundance_data_df) = gsub("[|]","_",colnames(test_abundance_data_df))
      colnames(test_abundance_data_df) = gsub(" ","_",colnames(test_abundance_data_df))
    }
  }
  metadata_train_df = read.csv(metadata_train)
  metadata_test_df = read.csv(metadata_test)
  rownames(metadata_train_df) = metadata_train_df$SubjectID
  rownames(metadata_test_df) = metadata_test_df$SubjectID
  train_sujbects1_df = read.table(train_sujbects1)
  train_sujbects2_df = read.table(train_sujbects2)
  train_sujbects3_df = read.table(train_sujbects3)
  test_subjects_df = read.table(test_subjects)
  
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


  train_sujbects1_to_keep = intersect(rownames(metadata_train_df),train_sujbects1_df[,1])
  train_sujbects2_to_keep = intersect(rownames(metadata_train_df),train_sujbects2_df[,1])
  train_sujbects3_to_keep = intersect(rownames(metadata_train_df),train_sujbects3_df[,1])
  test_subjects_to_keep = intersect(rownames(metadata_test_df),test_subjects_df[,1])
  
  train_sujbects1_df = train_sujbects1_df[match(train_sujbects1_to_keep,train_sujbects1_df[,1]),,drop=FALSE]
  train_sujbects2_df = train_sujbects2_df[match(train_sujbects2_to_keep,train_sujbects2_df[,1]),,drop=FALSE]
  train_sujbects3_df = train_sujbects3_df[match(train_sujbects3_to_keep,train_sujbects3_df[,1]),,drop=FALSE]
  test_subjects_df = test_subjects_df[match(test_subjects_to_keep,test_subjects_df[,1]),,drop=FALSE]
  
  # combine metadata and abundance data
  # first order metadata correctly
  if("microbiome"%in%feature_vec) {
    metadata_train_df = metadata_train_df[rownames(train_abundance_data_df),]
    metadata_test_df = metadata_test_df[rownames(test_abundance_data_df),]
  }

  categorical_features = c("fdr","num_persistent","number_autoantibodies","Maternal_PreEclampsia_Toxemia","Maternal_Weight_Gain_Aagaard","Sex","Maternal_Medication","Preterm","Celiac_Disease","Breastmilk_Ever","Maternal_Probiotic",
  "Maternal_Antibiotics","Geographical_Location","Birth_Mode","HLA_Category","Maternal_Diabetes","Maternal_BMI_Category","Race_ethnicity","Maternal_Diabetes_Medication","Maternal_Diabetes_Medication_Nam","delivery_simple",
  "brst_fed","hla_DR4","hla_DR3","hla_DR34","hla_5grps","hla_5grps_ref_DR44","t1d","single_persist_conf","two_persist_conf","three_persist_conf","two_or_more_persistent","GAD_pos","IA2A_pos","t1d_sero_control","Insulin",
  "Metformin","Glyburide","antihypertensives")

  # change don't know to NA 
  if("Maternal_PreEclampsia_Toxemia"%in%colnames(metadata_train_df)) {
    metadata_train_df$Maternal_PreEclampsia_Toxemia[which(metadata_train_df$Maternal_PreEclampsia_Toxemia=="Don't Know")] = NA
  }
  features_to_keep_vec = strsplit(features_to_keep,split=",")[[1]]
  feature_vec_temporary = feature_vec

  for (feature_temp in intersect(feature_vec,colnames(metadata_train_df))) {
    unique_stats = unique(metadata_train_df[,feature_temp])
    unique_stats = unique_stats[!is.na(unique_stats)]
    num_unique_stats = length(unique_stats)
    if(num_unique_stats>1) {
      if(feature_temp%in%categorical_features) {
        if(class(metadata_train_df[,feature_temp]) == "logical") {
           metadata_train_df[,feature_temp] = as.numeric(metadata_train_df[,feature_temp])
        }
        metadata_train_df[,feature_temp] = as.factor(metadata_train_df[,feature_temp])
        dummy_mat = model.matrix(as.formula(paste("~",feature_temp,sep="")), metadata_train_df)
        dummy_mat = dummy_mat[,-1,drop=FALSE]
        dummy_mat = dummy_mat[match(rownames(metadata_train_df),rownames(dummy_mat)),,drop=FALSE]
        metadata_train_df = metadata_train_df[,-match(feature_temp,colnames(metadata_train_df))]
        metadata_train_df = cbind(metadata_train_df,dummy_mat)
        feature_vec_temporary = feature_vec_temporary[-match(feature_temp,feature_vec_temporary)]
        feature_vec_temporary = c(feature_vec_temporary,colnames(dummy_mat))
        if(feature_temp%in%features_to_keep_vec) {
          features_to_keep_vec = features_to_keep_vec[-match(feature_temp,features_to_keep_vec)]
          features_to_keep_vec = c(features_to_keep_vec,colnames(dummy_mat))
        }
      }
    } else {
      feature_vec_temporary = feature_vec_temporary[-match(feature_temp,feature_vec_temporary)] 
    }
  }

  for (feature_temp in intersect(feature_vec,colnames(metadata_test_df))) {
    unique_stats = unique(metadata_test_df[,feature_temp])
    unique_stats = unique_stats[!is.na(unique_stats)]
    num_unique_stats = length(unique_stats)
    if(num_unique_stats>1) {
      if(feature_temp%in%categorical_features) {
        if(class(metadata_test_df[,feature_temp]) == "logical") {
          metadata_test_df[,feature_temp] = as.numeric(metadata_test_df[,feature_temp])
        }
        metadata_test_df[,feature_temp] = as.factor(metadata_test_df[,feature_temp])
	dummy_mat = model.matrix(as.formula(paste("~",feature_temp,sep="")), metadata_test_df)
	dummy_mat = dummy_mat[,-1,drop=FALSE]
	dummy_mat = dummy_mat[match(rownames(metadata_test_df),rownames(dummy_mat)),,drop=FALSE]
	metadata_test_df = metadata_test_df[,-match(feature_temp,colnames(metadata_test_df))]
	metadata_test_df = cbind(metadata_test_df,dummy_mat)
      }	
    }
  }

  feature_vec =	feature_vec_temporary

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

  if("microbiome"%in%feature_vec & length(feature_vec) > 1) {
    metadata_to_get = c("getsCondition",intersect(feature_vec,colnames(metadata_train_df)))
    all_train_data = cbind(train_abundance_data_df,metadata_train_df[,metadata_to_get])
    all_test_data = cbind(test_abundance_data_df,metadata_test_df[,metadata_to_get])
  } else if("microbiome"%in%feature_vec & length(feature_vec) == 1) {
    metadata_to_get = c("getsCondition")
    all_train_data = cbind(train_abundance_data_df,getsCondition=metadata_train_df[,metadata_to_get])
    all_test_data = cbind(test_abundance_data_df,getsCondition=metadata_test_df[,metadata_to_get])
  } else if(!"microbiome"%in%feature_vec) {
    metadata_to_get = c("getsCondition",intersect(feature_vec,colnames(metadata_train_df)))
    all_train_data = metadata_train_df[,metadata_to_get]
    all_test_data = metadata_test_df[,metadata_to_get]
  }

  #train_sujbects1_index1 = match(train_sujbects1_df[,1],rownames(all_train_data))
  #train_sujbects2_index2 = match(train_sujbects2_df[,1],rownames(all_train_data))
  #train_sujbects3_index3 = match(train_sujbects3_df[,1],rownames(all_train_data))

  #indexes = list(train_sujbects1_index1,train_sujbects2_index2,train_sujbects3_index3)
  
  foldIDs = rep(0,nrow(all_train_data))
  foldIDs[rownames(all_train_data)%in%train_sujbects1_df[,1]] = 1
  foldIDs[rownames(all_train_data)%in%train_sujbects2_df[,1]] =	2
  foldIDs[rownames(all_train_data)%in%train_sujbects3_df[,1]] =	3

  train_y = all_train_data[,"getsCondition"]
  test_y = all_test_data[,"getsCondition"]
  train_y = as.numeric(as.character(train_y))
  test_y = as.numeric(as.character(test_y))
  train_x = all_train_data[,-match(c("getsCondition"),colnames(all_train_data)),drop=FALSE]
  test_x = all_test_data[,-match(c("getsCondition"),colnames(all_test_data)),drop=FALSE] 

  train_x = as.matrix(train_x)
  test_x = as.matrix(test_x)
  
  penalty_factor = rep(1,ncol(train_x))
  penalty_factor[colnames(train_x)%in%features_to_keep_vec] = 0

  train_lasso = function(train_ctrl_counts,train_case_counts,all_train_data,all_test_data,train_x,train_y,test_x,test_y,weighting) {
    if(weighting == TRUE) {
      weight_0 = 1-sum(train_y==0)/length(train_y)
      weight_1 = 1-sum(train_y==1)/length(train_y)
      weights = rep(0,length(train_y))
      weights[train_y==0] = weight_0
      weights[train_y==1] = weight_1
    } else {
      weights = rep(1,length(train_y))
    }
   
    if(train_ctrl_counts>1 & train_case_counts>1) {
      if(ncol(train_x)>1) {
        if(number_threads>1) {
          library(doMC)
          registerDoMC(cores = number_threads)
          cv_output <- cv.glmnet(train_x, train_y, nfolds = 3,family="binomial",penalty.factor=penalty_factor,parallel=TRUE,foldid=foldIDs,weights=weights)
        } else {
          cv_output <- cv.glmnet(train_x, train_y, nfolds = 3,family="binomial",penalty.factor=penalty_factor,parallel=FALSE,foldid=foldIDs,weights=weights)
        }
        print("End lasso cross validation")
        best_lam <- cv_output$lambda.min
        train_coefs = coef(cv_output,s=best_lam)
        train_coefs = train_coefs[abs(train_coefs[,1])>0,,drop=FALSE]
        predict_prob_train = predict(cv_output, newx = train_x, s = "lambda.min",type="response")
        predict_prob_test = predict(cv_output, newx = test_x, s = "lambda.min",type="response")
        predictions_train = predict(cv_output, newx = train_x, s = "lambda.min",type="class")
        predictions_test = predict(cv_output, newx = test_x, s = "lambda.min",type="class")
        class(predictions_test) = "numeric"
        class(predictions_train) = "numeric"
      } else {
        all_train_data$getsCondition = as.numeric(as.character(all_train_data$getsCondition))
        all_test_data$getsCondition = as.numeric(as.character(all_test_data$getsCondition))
        my_formula <- paste0("~",paste(colnames(train_x), collapse = " + "))
        my_formula <- as.formula(paste("getsCondition",my_formula,sep=""))
        model_train_out = glm(my_formula,family=binomial(link='logit'),data=all_train_data,weights=weights)
        predict_prob_train = predict(model_train_out,newdata=all_train_data,type = "response")
        predict_prob_train = as.data.frame(predict_prob_train)
        predict_prob_test = predict(model_train_out,newdata=all_test_data,type = "response")
        predict_prob_test = as.data.frame(predict_prob_test)
        predictions_train = predict_prob_train
        predictions_test = predict_prob_test
        predictions_train[predictions_train[,1]>0.5,1] = 1
        predictions_train[predictions_train[,1]<=0.5,1] = 0
        predictions_test[predictions_test[,1]>0.5,1] = 1
        predictions_test[predictions_test[,1]<=0.5,1] = 0
        train_coefs = model_train_out$coefficients
      }

      predict_prob_train_with_labels = predict_prob_train
      predict_prob_test_with_labels = predict_prob_test
      colnames(predict_prob_train_with_labels) = c("case")
      colnames(predict_prob_test_with_labels) = c("case")
      predict_prob_train_with_labels = as.data.frame(predict_prob_train_with_labels)
      predict_prob_test_with_labels = as.data.frame(predict_prob_test_with_labels)
      predict_prob_train_with_labels$ctrl = 1-predict_prob_train_with_labels$case
      predict_prob_test_with_labels$ctrl = 1-predict_prob_test_with_labels$case
      predict_prob_test_with_labels$obs = as.character(all_test_data$getsCondition)
      predict_prob_train_with_labels$obs = as.character(all_train_data$getsCondition)
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
	train_pr_all = pr.curve(scores.class0 = predict_prob_train[,1],weights.class0=as.numeric(as.character(all_train_data$getsCondition)))
	train_pr_all = train_pr_all$auc.integral
	train_pr_vec = rep(0,100)
 	for(x in 1:100) {
          index_to_sample = sample(x = seq(1,nrow(predict_prob_train)), size  = nrow(predict_prob_train),replace=TRUE)	  
	  sampled_probs = predict_prob_train[index_to_sample,1]
	  sampled_labels = as.numeric(as.character(all_train_data$getsCondition))[index_to_sample]
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
        test_pr_all = pr.curve(scores.class0 = predict_prob_test[,1],weights.class0=as.numeric(as.character(all_test_data$getsCondition)))
        test_pr_all = test_pr_all$auc.integral
        test_pr_vec = rep(0,100)
        for(x in 1:100) {
          index_to_sample = sample(x = seq(1,nrow(predict_prob_test)), size  = nrow(predict_prob_test),replace=TRUE)
          sampled_probs = predict_prob_test[index_to_sample,1]
          sampled_labels = as.numeric(as.character(all_test_data$getsCondition))[index_to_sample]
          pr_test <- pr.curve(scores.class0 = sampled_probs, weights.class0 = sampled_labels)
          pr_test = pr_test$auc.integral
          test_pr_vec[x] = pr_test
        }
        test_pr_CI = paste(quantile(test_pr_vec,c(0.025,0.975),na.rm=TRUE),collapse="-")
        ROC.T_test["AUC-PR_Score"] = test_pr_all
        ROC.T_test["AUC-PR_CI"] = test_pr_CI
	final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec,predict_prob_test_with_labels)
	return(final_results)	
      },error = function(e) {
       ROC.T_train = rep(NA,30)
       ROC.T_test = rep(NA,nrow(test_roc_out_stats)*2)
       names(ROC.T_test) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
       names(ROC.T_train) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
       final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec,predict_prob_test_with_labels)
       return(final_results)
      })
    } else {
      ROC.T_train = rep(NA,30)
      ROC.T_test = rep(NA,30)
      names(ROC.T_train) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
      names(ROC.T_test) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
      train_coefs = NA
      predict_prob_test_with_labels = NA
      final_results = list(ROC.T_train,ROC.T_test,train_coefs,sample_counts,feature_vec,predict_prob_test_with_labels)
      return(final_results)
    }
  }

  poss_train_lasso = possibly(.f = train_lasso, otherwise = NA)
  weighted_results = poss_train_lasso(train_ctrl_counts,train_case_counts,all_train_data,all_test_data,train_x,train_y,test_x,test_y,weighting=TRUE)  
  unweighted_results = poss_train_lasso(train_ctrl_counts,train_case_counts,all_train_data,all_test_data,train_x,train_y,test_x,test_y,weighting=FALSE)

  return(list(weighted_results,unweighted_results))
}

model_results = run_lasso(test_abundance_data=test_abundance_data,train_abundance_data=train_abundance_data,metadata_train=metadata_train,metadata_test=metadata_test,train_sujbects1=train_sujbects1,train_sujbects2=train_sujbects2,train_sujbects3=train_sujbects3,test_subjects=test_subjects,feature_list=feature_list,features_to_keep=features_to_keep,microbiome_features=microbiome_features,abundance_type=abundance_type)

out_path = dirname(train_abundance_data)
out_prefix = basename(out_path)
if(abundance_type == "species") {
  out_prefix = paste(out_prefix,"species",sep="_")
} else if(abundance_type == "pathway") {
  out_prefix = paste(out_prefix,"pathway",sep="_")
}

if(length(feature_vec_1)>5 & "microbiome"%in%feature_vec_1) {
  output_weighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_","all_clinical_plus_microbiome","_weighted.rds",sep="")
  output_unweighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_","all_clinical_plus_microbiome",".rds",sep="")
} else if(length(feature_vec_1)>5 & !"time"%in%feature_vec_1) {
  output_weighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_","all_clinical_no_time","_weighted.rds",sep="")
  output_unweighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_","all_clinical_no_time",".rds",sep="")
} else if(length(feature_vec_1)>5) {
  output_weighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_","all_clinical","_weighted.rds",sep="")
  output_unweighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_","all_clinical",".rds",sep="")
} else {
  output_weighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_",feature_list,"_weighted.rds",sep="")
  output_unweighted = paste(out_path,'/',out_prefix,"output_lasso_logistic_regression_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_",feature_list,".rds",sep="")
}

saveRDS(model_results[[1]],output_weighted)
saveRDS(model_results[[2]],output_unweighted)


