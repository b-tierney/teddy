library(dplyr)
library(RNOmni)
library(glmnet)
library(survival)
library(timeROC)
library(caret)
library(parallel)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

input_folder = args[[1]]
proportion_cutoff = args[[2]]
thread_number = as.numeric(args[[3]])

metadata_train = read.csv(paste(input_folder,"_train_metadata.csv",sep=""))
metadata_test = read.csv(paste(input_folder,"_test_metadata.csv",sep=""))


proportion_cutoff = as.numeric(proportion_cutoff)
all_data = list.files(input_folder,full.names=TRUE)

baseline_files_locs = grep("_baseline",all_data)
if(length(baseline_files_locs)>0) {
all_data = all_data[-baseline_files_locs]
}

filtered_files = grep("filtered_",all_data)
if(length(filtered_files)>0) {
  all_data = all_data[-filtered_files]
}

abundance_train = all_data[grep("train_[a-z][a-z].rds",all_data)]
abundance_test = all_data[grep("test_[a-z][a-z].rds",all_data)]

abundance_data_train_test = data.frame(abundance_train,abundance_test)

#filter abundance data
filtered_abundance_data = lapply(seq(1,nrow(abundance_data_train_test)), function(x) {
  abundance_files = unlist(abundance_data_train_test[x,])
  train_abundance_file = abundance_files[1]
  test_abundance_file = abundance_files[2]
  print(train_abundance_file)
  dependent_variable_train = readRDS(train_abundance_file)
  dependent_variable_test = readRDS(test_abundance_file)
  colnames(dependent_variable_train)[1]='sampleID'
  colnames(dependent_variable_test)[1] = 'sampleID'
  dependent_variable_train$sampleID = as.character(dependent_variable_train$sampleID)
  dependent_variable_test$sampleID = as.character(dependent_variable_test$sampleID)
  metadata_abundance_samples_train = intersect(metadata_train$SubjectID,dependent_variable_train$sampleID)
  metadata_abundance_samples_test = intersect(metadata_test$SubjectID,dependent_variable_test$sampleID)
  dependent_variable_train = dependent_variable_train[match(metadata_abundance_samples_train,dependent_variable_train$sampleID),]  
  dependent_variable_test = dependent_variable_test[match(metadata_abundance_samples_test,dependent_variable_test$sampleID),]
  y_vars_train = dependent_variable_train$sampleID
  y_vars_test = dependent_variable_test$sampleID
  dependent_variable_train = dependent_variable_train[,-match("sampleID",colnames(dependent_variable_train))]
  dependent_variable_test = dependent_variable_test[,-match("sampleID",colnames(dependent_variable_test))]
  dependent_variable_train = as.matrix(dependent_variable_train)
  dependent_variable_test = as.matrix(dependent_variable_test)
  rownames(dependent_variable_train) = y_vars_train
  rownames(dependent_variable_test) = y_vars_test
  
  number_samples_gt_zero_abundance_each_gene = colSums(dependent_variable_train > 0)
  proportion_samples_with_gt_zero_abundance_each_gene = number_samples_gt_zero_abundance_each_gene/nrow(dependent_variable_train)
  genes_to_keep = which(proportion_samples_with_gt_zero_abundance_each_gene >= proportion_cutoff)
  dependent_variable_train = dependent_variable_train[,genes_to_keep,drop=FALSE]
  dependent_variable_test = dependent_variable_test[,colnames(dependent_variable_train),drop=FALSE]
  
  return(list(dependent_variable_train,dependent_variable_test))
})

train_data_filtered = lapply(filtered_abundance_data, function(abund) abund[[1]])
test_data_filtered = lapply(filtered_abundance_data,  function(abund) abund[[2]])


train_data_filtered = do.call("cbind",train_data_filtered)
test_data_filtered = do.call("cbind",test_data_filtered)

train_data_filtered = as.data.frame(train_data_filtered)
test_data_filtered = as.data.frame(test_data_filtered)

#trainout = paste(input_folder,"/",basename(input_folder),"filtered_abundance_train.csv",sep="")
#testout = paste(input_folder,"/",basename(input_folder),"filtered_abundance_test.csv",sep="")

trainout_transformed = paste(input_folder,"/",basename(input_folder),"filtered_transformed_abundance_train.csv",sep="")
testout_transformed = paste(input_folder,"/",basename(input_folder),"filtered_transformed_abundance_test.csv",sep="")

trainout_subjectNames = paste(input_folder,"/",basename(input_folder),"filtered_abundance_train_names.txt",sep="")
testout_subjectNames = paste(input_folder,"/",basename(input_folder),"filtered_abundance_test_names.txt",sep="")

train_data_filtered_transformed = do.call("cbind",apply(train_data_filtered,2, function(mycol) RankNorm(mycol),simplify=FALSE))

test_data_filtered_transformed = do.call("cbind",apply(test_data_filtered,2, function(mycol) RankNorm(mycol),simplify=FALSE))

train_data_filtered_transformed = as.data.frame(train_data_filtered_transformed)
test_data_filtered_transformed = as.data.frame(test_data_filtered_transformed)

#fwrite(train_data_filtered,trainout,row.names=FALSE)
#fwrite(test_data_filtered,testout,row.names=FALSE)

fwrite(train_data_filtered_transformed,trainout_transformed,row.names=FALSE)
fwrite(test_data_filtered_transformed,testout_transformed,row.names=FALSE)


write.table(rownames(train_data_filtered),file=trainout_subjectNames,col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(rownames(test_data_filtered),file=testout_subjectNames,col.names=FALSE,row.names=FALSE,quote=FALSE)
