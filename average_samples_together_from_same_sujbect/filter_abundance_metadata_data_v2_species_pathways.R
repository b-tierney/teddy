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

abundance_train_pathways = paste(input_folder,"/",basename(input_folder),"_train_TEDDY_pathabundance-cpm.tsv",sep="")
abundance_test_pathways = paste(input_folder,"/",basename(input_folder),"_test_TEDDY_pathabundance-cpm.tsv",sep="")
abundance_train_species  = paste(input_folder,"/",basename(input_folder),"_train_TEDDY_species_abundances.tsv",sep="")
abundance_test_species = paste(input_folder,"/",basename(input_folder),"_test_TEDDY_species_abundances.tsv",sep="")


abundance_data_train_test_species = c(abundance_train_species,abundance_test_species)
abundance_data_train_test_pathways = c(abundance_train_pathways,abundance_test_pathways)

filter_normalize_data = function(abundance_files) {
  train_abundance_file = abundance_files[1]
  test_abundance_file = abundance_files[2]
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
  train_data_filtered = dependent_variable_train[,genes_to_keep,drop=FALSE]
  test_data_filtered = dependent_variable_test[,colnames(train_data_filtered),drop=FALSE]
  train_data_filtered = as.data.frame(train_data_filtered)
  test_data_filtered = as.data.frame(test_data_filtered)
  trainout_transformed = paste(input_folder,"/",gsub(".tsv","",basename(train_abundance_file)),"_filtered_transformed_abundance_train.csv",sep="")
  testout_transformed = paste(input_folder,"/",gsub(".tsv","",basename(test_abundance_file)),"_filtered_transformed_abundance_test.csv",sep="")
  
  train_data_filtered_transformed = do.call("cbind",apply(train_data_filtered,2, function(mycol) RankNorm(mycol),simplify=FALSE))
  test_data_filtered_transformed = do.call("cbind",apply(test_data_filtered,2, function(mycol) RankNorm(mycol),simplify=FALSE))
  train_data_filtered_transformed = as.data.frame(train_data_filtered_transformed)
  test_data_filtered_transformed = as.data.frame(test_data_filtered_transformed)
  fwrite(train_data_filtered_transformed,trainout_transformed,row.names=FALSE)
  fwrite(test_data_filtered_transformed,testout_transformed,row.names=FALSE)

  trainout_subjectNames = paste(input_folder,"/",gsub(".tsv","",basename(train_abundance_file)),"_filtered_transformed_abundance_train_names.csv",sep="")
  testout_subjectNames = paste(input_folder,"/",gsub(".tsv","",basename(test_abundance_file)),"_filtered_transformed_abundance_test_names.csv",sep="")

  write.table(rownames(train_data_filtered),file=trainout_subjectNames,col.names=FALSE,row.names=FALSE,quote=FALSE)
  write.table(rownames(test_data_filtered),file=testout_subjectNames,col.names=FALSE,row.names=FALSE,quote=FALSE)
}

filter_normalize_data(abundance_data_train_test_species)
filter_normalize_data(abundance_data_train_test_pathways)