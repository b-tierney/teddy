library(caret)
library(data.table)
library(reticulate)
library(pROC)
library(timeROC)
library(glmnet)
library(doMC)
library(survival)
library(randomForestSRC)
library(pec)
library(pROC)
library(PRROC)
library(MLeval)
library(purrr)

tune.rfsrc = function (formula, data, mtryStart = ncol(data)/2, nodesizeTry = c(1:9, 
    seq(10, 100, by = 5)), ntreeTry = 100, sampsize = function(x) {
    min(x * 0.632, max(150, x^(3/4)))
}, nsplit = 1, stepFactor = 1.25, improve = 0.001, strikeout = 3, 
    maxIter = 25, trace = FALSE, doBest = FALSE, imbalanced=FALSE,...) 
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
    } else {
        ssize <- sampsize
    }

    if ((2 * ssize) < n & imbalanced == FALSE) { # only do if not imbalanced data
        tst <- sample(1:n, size = ssize, replace = FALSE)
        trn <- setdiff(1:n, tst)
        newdata <- data[tst, , drop = FALSE]
    } else {
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
        } else {
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
                } else {
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
                } else {
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
                if(is.na(Improve)) { # this occurs rarely but good to correct for
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

run_RF = function(test_abundance_data,train_abundance_data,metadata_train,metadata_test,train_sujbects1,train_sujbects2,train_sujbects3,test_subjects,feature_list,features_to_keep,microbiome_features,abundance_type) {

  
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
    
    #train_abundance_data_df_mat = as.matrix(train_abundance_data_df)
    #cor_mat = cor(train_abundance_data_df_mat)
    #correlations_found = findCorrelation(cor_mat,cutoff=0.95)
    #if(length(correlations_found>0)) {
    #  train_abundance_data_df = train_abundance_data_df[,-correlations_found]
    #  test_abundance_data_df = test_abundance_data_df[,-correlations_found]
    #}
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
  #all_train_data$status = as.numeric(as.character(all_train_data$status))
  #all_test_data$status	= as.numeric(as.character(all_test_data$status))

  features_no_time_no_status = colnames(all_train_data)
  features_no_time_no_status= features_no_time_no_status[-match(c("time","status"),features_no_time_no_status)]

  my_formula <- as.formula(paste0("status ~",paste(features_no_time_no_status, collapse = " + ")))
  #my_formula <- as.formula(paste("status",my_formula,sep=""))

  train_rf = function(my_formula,all_train_data,all_test_data, weightingMethod=NA) {
  
    # using OOB error tune random forest with training data
    if(is.na(weightingMethod)) {
      ntree100 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 100)
      ntree500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 500)
      ntree1000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1000)
      ntree1500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1500)
      ntree2000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 2000)
    } else if(weightingMethod == "BRF") {
      ntree100 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 100,case.wt = randomForestSRC:::make.wt(all_train_data$status),perf.type = "gmean",imbalanced=TRUE)
      ntree500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 500,case.wt = randomForestSRC:::make.wt(all_train_data$status),perf.type = "gmean",imbalanced=TRUE)
      ntree1000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1000,case.wt = randomForestSRC:::make.wt(all_train_data$status),perf.type = "gmean",imbalanced=TRUE)
      ntree1500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1500,case.wt = randomForestSRC:::make.wt(all_train_data$status),perf.type = "gmean",imbalanced=TRUE)
      ntree2000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 2000,case.wt = randomForestSRC:::make.wt(all_train_data$status),perf.type = "gmean",imbalanced=TRUE)
   } else if(weightingMethod == "RFQ") {
      ntree100 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 100,rfq=TRUE,perf.type="gmean",imbalanced=TRUE)
      ntree500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 500,rfq=TRUE,perf.type="gmean",imbalanced=TRUE)
      ntree1000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1000,rfq=TRUE,perf.type="gmean",imbalanced=TRUE)
      ntree1500 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 1500,rfq=TRUE,perf.type="gmean",imbalanced=TRUE)
      ntree2000 = tune.rfsrc(my_formula, data=all_train_data,ntreeTry = 2000,rfq=TRUE,perf.type="gmean",imbalanced=TRUE)
   }


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
    
    if(is.na(weightingMethod)) {
      rf_res = rfsrc(my_formula,data=all_train_data,ntree=optimal_ntree,mtry=optimal_mtry,nodesize=optimal_nodesize,importance=TRUE)
    } else if(weightingMethod=="BRF") {
      rf_res = imbalanced(my_formula,data=all_train_data,method="brf",splitrule = "auc",ntree=optimal_ntree,mtry=optimal_mtry,nodesize=optimal_nodesize,importance=TRUE)
    } else if(weightingMethod=="RFQ") {
      rf_res = imbalanced(my_formula,data=all_train_data,method="rfq",splitrule = "auc",ntree=optimal_ntree,mtry=optimal_mtry,nodesize=optimal_nodesize,importance=TRUE)
    }
    importance_variables = rf_res$importance

    print("End cross validation")

    # get predicted probabilities for having condition
    predict_prob_train = rf_res$predicted
    predict_prob_test_mod = predict(rf_res,newdata=all_test_data)
    predict_prob_test = predict_prob_test_mod$predicted
  
    predictions_train = colnames(rf_res$predicted)[apply(rf_res$predicted,1, which.max)]
    predictions_test = colnames(predict_prob_test_mod$predicted)[apply(predict_prob_test_mod$predicted,1, which.max)]

    
    predict_prob_train_with_labels = as.data.frame(predict_prob_train)
    predict_prob_test_with_labels = as.data.frame(predict_prob_test)
    colnames(predict_prob_train_with_labels) = c("ctrl","case")
    colnames(predict_prob_test_with_labels) = c("ctrl","case")
    predict_prob_train_with_labels$obs = as.character(all_train_data$status)
    predict_prob_train_with_labels$obs[predict_prob_train_with_labels$obs==0] = "ctrl"
    predict_prob_train_with_labels$obs[predict_prob_train_with_labels$obs==1] =	"case"
    predict_prob_test_with_labels$obs = as.character(all_test_data$status)
    predict_prob_test_with_labels$obs[predict_prob_test_with_labels$obs==0] =	"ctrl"
    predict_prob_test_with_labels$obs[predict_prob_test_with_labels$obs==1] = "case"

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
      train_pr_all = pr.curve(scores.class0 = predict_prob_train[,"1"],weights.class0=as.numeric(as.character(all_train_data$status)))
      train_pr_all = train_pr_all$auc.integral
      train_pr_vec = rep(0,100)
      for(x in 1:100) {
        index_to_sample = sample(x = seq(1,nrow(predict_prob_train)), size  = nrow(predict_prob_train),replace=TRUE)
	sampled_probs = predict_prob_train[index_to_sample,"1"]
	sampled_labels = as.numeric(as.character(all_train_data$status))[index_to_sample]
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
     test_pr_all = pr.curve(scores.class0 = predict_prob_test[,"1"],weights.class0=as.numeric(as.character(all_test_data$status)))
     test_pr_all = test_pr_all$auc.integral
     test_pr_vec = rep(0,100)
     for(x in 1:100) {
       index_to_sample = sample(x = seq(1,nrow(predict_prob_test)), size  = nrow(predict_prob_test),replace=TRUE)
       sampled_probs = predict_prob_test[index_to_sample,"1"]
       sampled_labels = as.numeric(as.character(all_test_data$status))[index_to_sample]
       pr_test <- pr.curve(scores.class0 = sampled_probs, weights.class0 = sampled_labels)
       pr_test = pr_test$auc.integral
       test_pr_vec[x] = pr_test
     }
     test_pr_CI = paste(quantile(test_pr_vec,c(0.025,0.975),na.rm=TRUE),collapse="-")
     ROC.T_test["AUC-PR_Score"] = test_pr_all
     ROC.T_test["AUC-PR_CI"] = test_pr_CI
     return(list(ROC.T_train,ROC.T_test,importance_variables,sample_counts,feature_vec))

   },error = function(e) {
     ROC.T_test = rep(NA,30)
     ROC.T_train = rep(NA,30)
     names(ROC.T_train) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
     names(ROC.T_test) = c("SENS_Score","SENS_CI","SPEC_Score","SPEC_CI","MCC_Score","MCC_CI","Informedness_Score","Informedness_CI","PREC_Score","PREC_CI","NPV_Score","NPV_CI","FPR_Score","FPR_CI","F1_Score","F1_CI","TP_Score","TP_CI","FP_Score","FP_CI","TN_Score","TN_CI","FN_Score","FN_CI","AUC-ROC_Score","AUC-ROC_CI","AUC-PR_Score","AUC-PR_CI","AUC-PRG_Score","AUC-PRG_CI")
     return(list(ROC.T_train,ROC.T_test,importance_variables,sample_counts,feature_vec))
   })
  }
  poss_train_rf = possibly(.f = train_rf, otherwise = NA)
  train_test_results_no_weight = poss_train_rf(my_formula,all_train_data,all_test_data,weightingMethod=NA)    
  train_test_results_BRF = poss_train_rf(my_formula,all_train_data,all_test_data,weightingMethod="BRF")
  train_test_results_RFQ = poss_train_rf(my_formula,all_train_data,all_test_data,weightingMethod="RFQ")

  return(list(train_test_results_no_weight,train_test_results_BRF,train_test_results_RFQ))
}



model_results = run_RF(test_abundance_data=test_abundance_data,train_abundance_data=train_abundance_data,metadata_train=metadata_train,metadata_test=metadata_test,train_sujbects1=train_sujbects1,train_sujbects2=train_sujbects2,train_sujbects3=train_sujbects3,test_subjects=test_subjects,feature_list=feature_list,features_to_keep=features_to_keep,microbiome_features=microbiome_features,abundance_type=abundance_type)

out_path = dirname(train_abundance_data)
out_prefix = basename(out_path)

if(abundance_type == "species") {
  out_prefix = paste(out_prefix,"species",sep="_")
} else if(abundance_type == "pathway") {
  out_prefix = paste(out_prefix,"pathway",sep="_")
}


output_no_weights = paste(out_path,'/',out_prefix,"output_random_forest_binary_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_",feature_list,".rds",sep="")

output_BRF = paste(out_path,'/',out_prefix,"output_random_forest_binary_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_",feature_list,"_BRF_weighted.rds",sep="")

output_RFQ = paste(out_path,'/',out_prefix,"output_random_forest_binary_",loss_function,"_loss","_microbiome_selection_method_",microbiome_features,"_feature_list_",feature_list,"_RFQ_weighted.rds",sep="")

saveRDS(model_results[[1]],output_no_weights)
saveRDS(model_results[[2]],output_BRF)
saveRDS(model_results[[3]],output_RFQ)

