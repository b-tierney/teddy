library(caret)
library(data.table)
library(pROC)
library(timeROC)
library(glmnet)
library(doMC)
library(survival)
library(randomForestSRC)
library(pec)
library(pROC)
library(caTools)
library(parallel)

## modification to randomForestSRC function to correct for NA values

tune.rfsrc = function (formula, data, mtryStart = ncol(data)/2, nodesizeTry = c(1:9, 
    seq(10, 100, by = 5)), ntreeTry = 100, sampsize = function(x) {
    min(x * 0.632, max(150, x^(3/4)))
}, nsplit = 1, stepFactor = 1.25, improve = 0.001, strikeout = 3, 
    maxIter = 25, trace = FALSE, doBest = FALSE, ...) 
{
    if (improve < 0) {
        stop("improve must be non-negative.")
    }
    if (stepFactor <= 1) {
        stop("stepFactor must be great than 1.")
    }
    if (missing(formula)) {
        stop("a formula must be supplied (only supervised forests allowed).")
    }
    stump <- rfsrc(formula, data, nodedepth = 0, perf.type = "none", 
        save.memory = TRUE, ntree = 1, splitrule = "random")
    n <- stump$n
    yvar.names <- stump$yvar.names
    data <- data.frame(stump$yvar, stump$xvar)
    colnames(data)[1:length(yvar.names)] <- yvar.names
    rm(stump)
    if (is.function(sampsize)) {
        ssize <- sampsize(n)
    }
    else {
        ssize <- sampsize
    }
    if ((2 * ssize) < n) {
        tst <- sample(1:n, size = ssize, replace = FALSE)
        trn <- setdiff(1:n, tst)
        newdata <- data[tst, , drop = FALSE]
    }
    else {
        trn <- 1:n
        newdata <- NULL
    }
    res <- list()
    counter1 <- 0
    for (nsz in nodesizeTry) {
        counter1 <- counter1 + 1
        if (is.null(newdata)) {
            o <- rfsrc.fast(formula, data, ntree = ntreeTry, 
                mtry = mtryStart, nodesize = nsz, sampsize = ssize, 
                nsplit = nsplit, ...)
            errorOld <- mean(get.mv.error(o, TRUE), na.rm = TRUE)
        }
        else {
            o <- rfsrc.fast(formula, data[trn, , drop = FALSE], 
                ntree = ntreeTry, mtry = mtryStart, nodesize = nsz, 
                sampsize = ssize, nsplit = nsplit, forest = TRUE, 
                perf.type = "none", save.memory = FALSE, ...)
            errorOld <- mean(get.mv.error(predict(o, newdata, 
                perf.type = "default"), TRUE), na.rm = TRUE)
        }
        mtryStart <- o$mtry
        mtryMax <- length(o$xvar.names)
        if (trace) {
            cat("nodesize = ", nsz, " mtry =", mtryStart, "error =", 
                paste(100 * round(errorOld, 4), "%", sep = ""), 
                "\n")
        }
        if (mean(o$leaf.count) <= 2) {
            break
        }
        oobError <- list()
        oobError[[1]] <- errorOld
        names(oobError)[1] <- mtryStart
        for (direction in c("left", "right")) {
            if (trace) 
                cat("Searching", direction, "...\n")
            Improve <- 1.1 * improve
            mtryBest <- mtryStart
            mtryCur <- mtryStart
            counter2 <- 1
            strikes <- 0
            while (counter2 <= maxIter && (Improve >= improve || 
                (Improve < 0 && strikes < strikeout))) {
                counter2 <- counter2 + 1
                if (Improve < 0) {
                  strikes <- strikes + 1
                }
                mtryOld <- mtryCur
                if (direction == "left") {
                  mtryCur <- max(1, min(ceiling(mtryCur/stepFactor), 
                    mtryCur - 1))
                }
                else {
                  mtryCur <- min(mtryMax, max(floor(mtryCur * 
                    stepFactor), mtryCur + 1))
                }
                if (mtryCur == mtryOld) {
                  break
                }
                if (is.null(newdata)) {
                  errorCur <- mean(get.mv.error(rfsrc.fast(formula, 
                    data, ntree = ntreeTry, mtry = mtryCur, nodesize = nsz, 
                    sampsize = ssize, nsplit = nsplit, ...), 
                    TRUE), na.rm = TRUE)
                }
                else {
                  errorCur <- mean(get.mv.error(predict(rfsrc.fast(formula, 
                    data[trn, , drop = FALSE], ntree = ntreeTry, 
                    mtry = mtryCur, nodesize = nsz, sampsize = ssize, 
                    nsplit = nsplit, forest = TRUE, perf.type = "none", 
                    save.memory = FALSE, ...), newdata, perf.type = "default"), 
                    TRUE), na.rm = TRUE)
                }
                if (trace) {
                  cat("nodesize = ", nsz, " mtry =", mtryCur, 
                    " error =", paste(100 * round(errorCur, 4), 
                      "%", sep = ""), "\n")
                }
                oobError[[as.character(mtryCur)]] <- errorCur
                Improve <- 1 - errorCur/errorOld
                if(is.na(Improve)) {
                  Improve = 1
                }
                if (trace) {
                  cat(Improve, improve, "\n")
                }
                if (Improve > improve) {
                  errorOld <- errorCur
                  mtryBest <- mtryCur
                }
            }
        }
        mtry <- sort(as.numeric(names(oobError)))
        err <- unlist(oobError[as.character(mtry)])
        res[[counter1]] <- cbind(nodesize = nsz, mtry = mtry, 
            err = err)
    }
    if (is.null(res)) {
        stop("NULL results - something is wrong, check parameter (tuning) settings\n")
    }
    res <- do.call(rbind, res)
    res <- res[order(res[, 1], res[, 2]), ]
    rownames(res) <- 1:nrow(res)
    opt.idx <- which.min(res[, 3])
    rf <- NULL
    if (doBest) {
        rf <- rfsrc.fast(formula, data, mtry = res[opt.idx, 2], 
            nodesize = res[opt.idx, 1], sampsize = sampsize, 
            nsplit = nsplit)
    }
    list(results = res, optimal = res[opt.idx, -3], rf = rf)
}


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

#options(mc.cores = number_threads)

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

if(grepl("species",basename(test_abundance_data))) {
  abundance_type = "species"
} else if(grepl("pathabundance",basename(test_abundance_data))) {
  abundance_type = "pathway"
} else {
  abundance_type = "gene"
}

run_survival_RF = function(test_abundance_data,train_abundance_data,metadata_train,metadata_test,train_sujbects1,train_sujbects2,train_sujbects3,test_subjects,feature_list,features_to_keep,microbiome_features,abundance_type) {

  
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
    sig_genes_df = names(sort(sig_genes_df[[2]]))
    if(length(sig_genes_df)>100) {
      sig_genes_df = sig_genes_df[1:100]
    }
    train_abundance_data_df = train_abundance_data_df[,sig_genes_df,drop=FALSE]
    test_abundance_data_df = test_abundance_data_df[,sig_genes_df,drop=FALSE]

#    train_abundance_data_df_mat = as.matrix(train_abundance_data_df)
#    cor_mat = cor(train_abundance_data_df_mat)
#    correlations_found = findCorrelation(cor_mat,cutoff=0.95)
#    if(length(correlations_found>0)) {
#       train_abundance_data_df = train_abundance_data_df[,-correlations_found]
#       test_abundance_data_df = test_abundance_data_df[,-correlations_found]
#    }
   }

  if(abundance_type != "gene") {
    colnames(train_abundance_data_df) = gsub("[|]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub(" ","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub(":","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("-","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[(]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[)]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[']","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[,]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^1","ONE_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^2","TWO_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^3","THREE_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^4","FOUR_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^5","FIVE_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^6","SIX_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^7","SEVEN_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^8","EIGHT_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("^9","NINE_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("&","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub(";","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[[]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[]]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[/]","_",colnames(train_abundance_data_df))
    colnames(train_abundance_data_df) = gsub("[+]","_",colnames(train_abundance_data_df))

    colnames(test_abundance_data_df) = gsub("[|]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub(" ","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub(":","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("-","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[(]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[)]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[']","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[,]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^1","ONE_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^2","TWO_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^3","THREE_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^4","FOUR_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^5","FIVE_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^6","SIX_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^7","SEVEN_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^8","EIGHT_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("^9","NINE_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("&","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub(";","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[[]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[]]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[/]","_",colnames(test_abundance_data_df))
    colnames(test_abundance_data_df) = gsub("[+]","_",colnames(test_abundance_data_df))

  }
  


  metadata_train_df = read.csv(metadata_train)
  metadata_test_df = read.csv(metadata_test)
  rownames(metadata_train_df) = metadata_train_df$SubjectID
  rownames(metadata_test_df) = metadata_test_df$SubjectID
  #train_sujbects1_df = read.table(train_sujbects1)
  #train_sujbects2_df = read.table(train_sujbects2)
  #train_sujbects3_df = read.table(train_sujbects3)
  test_subjects_df = read.table(test_subjects)
  
  train_counts = table(metadata_train_df$getsCondition)
  test_counts = table(metadata_test_df$getsCondition)
  
  train_ctrl_counts = train_counts["0"]
  train_case_counts = train_counts["1"]
  test_ctrl_counts = test_counts["0"]
  test_case_counts = test_counts["1"]

  sample_counts = c(train_ctrl=train_ctrl_counts,train_case=train_case_counts,total_train=train_ctrl_counts+train_case_counts,test_ctrl=test_ctrl_counts,test_case=test_case_counts,total_test=test_ctrl_counts+test_case_counts)
  names(sample_counts) = c("train_ctrl","train_case","total_train","test_ctrl","test_case","total_test")


  #train_sujbects1_to_keep = intersect(rownames(metadata_train_df),train_sujbects1_df[,1])
  #train_sujbects2_to_keep = intersect(rownames(metadata_train_df),train_sujbects2_df[,1])
  #train_sujbects3_to_keep = intersect(rownames(metadata_train_df),train_sujbects3_df[,1])
  test_subjects_to_keep = intersect(rownames(metadata_test_df),test_subjects_df[,1])
  
  #train_sujbects1_df = train_sujbects1_df[match(train_sujbects1_to_keep,train_sujbects1_df[,1]),]
  #train_sujbects2_df = train_sujbects2_df[match(train_sujbects2_to_keep,train_sujbects2_df[,1]),]
  #train_sujbects3_df = train_sujbects3_df[match(train_sujbects3_to_keep,train_sujbects3_df[,1]),]
  test_subjects_df = test_subjects_df[match(test_subjects_to_keep,test_subjects_df[,1]),,drop=FALSE]
  
  # combine metadata and abundance data
  # first order metadata correctly
  metadata_train_df = metadata_train_df[rownames(train_abundance_data_df),]
  metadata_test_df = metadata_test_df[rownames(test_abundance_data_df),]

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

  if("microbiome"%in%feature_vec & length(feature_vec) > 1) {
    metadata_to_get = c("time","getsCondition",intersect(feature_vec,colnames(metadata_train_df)))
    all_train_data = cbind(train_abundance_data_df,metadata_train_df[,metadata_to_get])
    all_test_data = cbind(test_abundance_data_df,metadata_test_df[,metadata_to_get])
  } else if("microbiome"%in%feature_vec & length(feature_vec) == 1) {
    metadata_to_get = c("time","getsCondition")
    all_train_data = cbind(train_abundance_data_df,metadata_train_df[,metadata_to_get])
    all_test_data = cbind(test_abundance_data_df,metadata_test_df[,metadata_to_get])
  } else if(!"microbiome"%in%feature_vec) {
    metadata_to_get = c("time","getsCondition",intersect(feature_vec,colnames(metadata_train_df)))
    all_train_data = metadata_train_df[,metadata_to_get]
    all_test_data = metadata_test_df[,metadata_to_get]
  }

  features_to_keep_vec = strsplit(features_to_keep,split=",")[[1]]
  if("number_autoantibodies"%in%features_to_keep_vec & num_unique_antibodies>1) {
    metadata_cols_to_take = colnames(metadata_train_df)[grep("number_autoantibodies",colnames(metadata_train_df))]
    features_to_keep_vec = features_to_keep_vec[-match("number_autoantibodies",features_to_keep_vec)]
    features_to_keep_vec = c(features_to_keep_vec,metadata_cols_to_take)
  } else if("number_autoantibodies"%in%features_to_keep_vec & num_unique_antibodies==1) {
    features_to_keep_vec = features_to_keep_vec[-match("number_autoantibodies",features_to_keep_vec)]  
  }

  #train_sujbects1_index1 = match(train_sujbects1_df[,1],rownames(all_train_data))
  #train_sujbects2_index2 = match(train_sujbects2_df[,1],rownames(all_train_data))
  #train_sujbects3_index3 = match(train_sujbects3_df[,1],rownames(all_train_data))

  #indexes = list(train_sujbects1_index1,train_sujbects2_index2,train_sujbects3_index3)
  
  #foldIDs = rep(0,nrow(all_train_data))
  #foldIDs[rownames(all_train_data)%in%train_sujbects1_df[,1]] = 1
  #foldIDs[rownames(all_train_data)%in%train_sujbects2_df[,1]] = 2
  #foldIDs[rownames(all_train_data)%in%train_sujbects3_df[,1]] = 3

  colnames(all_train_data)[match("getsCondition",colnames(all_train_data))] = "status"
  colnames(all_test_data)[match("getsCondition",colnames(all_test_data))] = "status"
  
  all_train_data$time = as.numeric(all_train_data$time)
  all_test_data$time = as.numeric(all_test_data$time)
  all_train_data$status = as.numeric(as.character(all_train_data$status))
  all_test_data$status	= as.numeric(as.character(all_test_data$status))

  features_no_time_no_status = colnames(all_train_data)
  features_no_time_no_status= features_no_time_no_status[-match(c("time","status"),features_no_time_no_status)]

  my_formula <- paste0("~",paste(features_no_time_no_status, collapse = " + "))
  my_formula <- as.formula(paste("Surv(time,status)",my_formula,sep=""))

  # using OOB error tune random forest with training data
  ntree100 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 100)
  ntree500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 500)
  ntree1000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1000)
  ntree1500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1500)
  ntree2000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 2000)

  # get min error from each
  ntree100_res = as.data.frame(ntree100$results)
  ntree100_res$ntree = 100
  ntree500_res = as.data.frame(ntree500$results)
  ntree500_res$ntree = 500
  ntree1000_res = as.data.frame(ntree1000$results)
  ntree1000_res$ntree = 1000
  ntree1500_res = as.data.frame(ntree1500$results)
  ntree1500_res$ntree = 1500
  ntree2000_res = as.data.frame(ntree2000$results)
  ntree2000_res$ntree = 2000
  ntree_all_res = rbind(ntree100_res,ntree500_res,ntree1000_res,ntree1500_res,ntree2000_res)
  min_error_index = which.min(ntree_all_res$err)
  optimal_nodesize = ntree_all_res[min_error_index,1]
  optimal_mtry = ntree_all_res[min_error_index,2]
  optimal_ntree = ntree_all_res[min_error_index,4]

  rf_res = rfsrc(my_formula,data=all_train_data,ntree=optimal_ntree,mtry=optimal_mtry,nodesize=optimal_nodesize,importance=TRUE)
  importance_variables = rf_res$importance


  predicted_risk_train = rf_res$predicted
  # get predicted risk for test data
  predicted_risk_test_all = predict(rf_res,newdata=all_test_data)      
  predicted_risk_test = predicted_risk_test_all$predicted

  # get AUCs
  ROC.T_train <- tryCatch({
     ROC.T_train <- timeROC(T = all_train_data[,"time"],
        delta = all_train_data[,"status"],marker = predicted_risk_train,
        cause = 1, weighting = "marginal",
        times = day_end,
        iid = TRUE) # set to false for now
  },error = function(e) {
     ROC.T_train <- timeROC(T = all_train_data[,"time"],
        delta = all_train_data[,"status"],marker = predicted_risk_train,
        cause = 1, weighting = "marginal",
        times = day_end,
        iid = FALSE) # set to false for now
  })
  ROC.T_test <-  tryCatch({
     ROC.T_test <- timeROC(T = all_test_data[,"time"],
          delta = all_test_data[,"status"],marker = predicted_risk_test,
          cause = 1, weighting = "marginal",
          times = day_end,
          iid = TRUE) # set to false for now
  },error = function(e) {
     ROC.T_test <- timeROC(T = all_test_data[,"time"],
          delta = all_test_data[,"status"],marker = predicted_risk_test,
          cause = 1, weighting = "marginal",
          times = day_end,
          iid = FALSE) # set to false for now
  })

  #now get AUCs for precision recall
  calc_precision_recall = function(risks,all_data) {
    if(min(risks) >0) {
      risks2 = c(0,risks)
    } else {
      risks2 = risks
    }
    train_precision_recall_vals = lapply(risks2,function(relative_risk) {
      train_stats_temp = SeSpPPVNPV(cutpoint=relative_risk,T=all_data[,"time"],delta = all_data[,"status"],marker = risks,cause = 1, weighting = "marginal",times = day_end,iid = FALSE)
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

  pr_all_train = calc_precision_recall(predicted_risk_train,all_train_data)
  pr_all_test = calc_precision_recall(predicted_risk_test,all_test_data)

  # now bootstrap. with replacement!!!
  pr_mat_train = matrix(NA,ncol=4,nrow=100)
  for (x in 1:100) {
    index_to_sample = sample(x = seq(1,length(predicted_risk_train)), size  = length(predicted_risk_train),replace=TRUE)
    sampled_risk = predicted_risk_train[index_to_sample]
    sampled_data = all_train_data[index_to_sample,]
    pr_temp = calc_precision_recall(sampled_risk,sampled_data)
    pr_mat_train[x,] = pr_temp
  }

  pr_mat_test = matrix(NA,ncol=4,nrow=100)
  for (x in 1:100) {
    index_to_sample = sample(x = seq(1,length(predicted_risk_test)), size = length(predicted_risk_test),replace=TRUE)
    sampled_risk = predicted_risk_test[index_to_sample]
    sampled_data = all_test_data[index_to_sample,]
    pr_temp = calc_precision_recall(sampled_risk,sampled_data)
    pr_mat_train[x,] = pr_temp
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
    myquantiles	= quantile(mycol,c(0.025,0.975),na.rm=TRUE)
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

  # get C-index
  if(is.null(predicted_risk_test_all$yvar[,1])) {
    c_index_test = NA
  } else {
    c_index_test = get.cindex(predicted_risk_test_all$yvar[,1], predicted_risk_test_all$yvar[,2], predicted_risk_test_all$predicted)
  }
  if(is.null(rf_res$yvar[,1])) {
    c_index_train = NA
  } else {
    c_index_train = get.cindex(rf_res$yvar[,1], rf_res$yvar[,2], rf_res$predicted.oob)
  }
  c_indexes = c(train=c_index_train,test=c_index_test)

  final_results = list(ROC.T_train,ROC.T_test,all_pr_cis,c_indexes,importance_variables,sample_counts,feature_vec)

  return(final_results)
}

model_results = run_survival_RF(test_abundance_data=test_abundance_data,train_abundance_data=train_abundance_data,metadata_train=metadata_train,metadata_test=metadata_test,train_sujbects1=train_sujbects1,train_sujbects2=train_sujbects2,train_sujbects3=train_sujbects3,test_subjects=test_subjects,feature_list=feature_list,features_to_keep=features_to_keep,microbiome_features=microbiome_features,abundance_type=abundance_type)

out_path = dirname(train_abundance_data)
out_prefix = basename(out_path)

if(abundance_type == "species") {
  out_prefix = paste(out_prefix,"species",sep="_")
} else if(abundance_type == "pathway") {
  out_prefix = paste(out_prefix,"pathway",sep="_")
}


output = paste(out_path,'/',out_prefix,"output_random_forest_time_to_event_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_",feature_list,".rds",sep="")


saveRDS(model_results,output)