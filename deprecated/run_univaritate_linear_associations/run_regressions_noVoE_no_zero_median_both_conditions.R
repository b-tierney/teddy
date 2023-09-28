library(dplyr)
#library(quantvoe)

args = commandArgs(trailingOnly=TRUE)

if(!file.exists(paste(args[[3]],"/",'association_output_full_',args[[1]],args[[2]],sep=''))) {
  dependent_variables = readRDS(args[[1]])
  independent_variables = readRDS(args[[2]])

  #min_nonzero_val = readRDS("min_val")
  
  colnames(dependent_variables)[1]='sampleID'
  colnames(independent_variables)[1]='sampleID'

  dependent_variables$sampleID = as.character(dependent_variables$sampleID)
  independent_variables$sampleID = as.character(independent_variables$sampleID)

  independent_variables = independent_variables[match(dependent_variables$sampleID,independent_variables$sampleID),]

  temp = dependent_variables[,-match("sampleID",colnames(dependent_variables))]
  temp = as.matrix(temp)

  temp_t1d_only= temp[independent_variables$condition==1,]
  temp_healthy_only = temp[independent_variables$condition==0,]
  median_t1d = apply(temp_t1d_only,2,median)
  median_healthy = apply(temp_healthy_only,2,median)
  genes_to_remove = which(median_healthy == 0 &median_t1d == 0)

  if(length(genes_to_remove) > 0) {
    temp = temp[,-genes_to_remove]
  }
  # find min non zero value in each column
  #temp_sum = temp + min_nonzero_val
  #temp_sum = temp + 0.000001
  #temp_logged = log(temp_sum)

  #dependent_variables = bind_cols(independent_variables %>% select(sampleID),temp_logged)

  condition = independent_variables$condition  
  # get list of gene names
  y_s = colnames(temp)
  # run regressions
  my_output = lapply(y_s, function(x) {  
    dependent_var = temp[,x]
    # get smallest non zero value
    minimum_nonzero = min(dependent_var[dependent_var>0])
    dependent_var = dependent_var + minimum_nonzero
    dependent_var_logged = log(dependent_var)
    true_model = glm(dependent_var_logged ~ condition,family=gaussian())
    out_info = summary(true_model)$coefficients["condition",]
    out_info['feature'] = x
    return(out_info)
  })
  my_output = bind_rows(my_output)
  colnames(my_output) = c("estimate","std.error","statistic","p.value","feature")
  my_output = cbind(term=rep("condition",nrow(my_output)),my_output)
  my_output$estimate = as.numeric(my_output$estimate)
  my_output$std.error = as.numeric(my_output$std.error)
  my_output$statistic = as.numeric(my_output$statistic)
  my_output$p.value = as.numeric(my_output$p.value)
  my_output$bonferroni = p.adjust(my_output$p.value,method="bonferroni")
  my_output$BH = p.adjust(my_output$p.value,method="BH")
  my_output$BY = p.adjust(my_output$p.value,method="BY")
  my_output = as_tibble(my_output)

  saveRDS(my_output,paste(args[[3]],"/",'association_output_full_',args[[1]],args[[2]],sep=''))
}
