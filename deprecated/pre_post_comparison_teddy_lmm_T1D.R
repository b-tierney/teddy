#process teddy metadata

library(lme4)
library(tidyverse)
library(lmerTest)
library(broom.mixed)
options(warn=2)

metadata = read.csv('teddy_metadata_20190821.csv',stringsAsFactors = FALSE)

metadata = metadata %>% select(Run,maskid,age_at_collection,two_or_more_persistent,IA_Outcome,t1d,T1D_Outcome)

lookup <- c("Before" = 0, "After" = 1, "Never" = 0)

metadata$IA_Outcome <- lookup[metadata$IA_Outcome]
metadata$T1D_Outcome <- lookup[metadata$T1D_Outcome]
metadata$t1d <- as.numeric(metadata$t1d)
metadata$two_or_more_persistent <- as.numeric(metadata$two_or_more_persistent)

colnames(metadata)=c('sampleID','subjectID','age','seroconverted_ever','seroconverted_at_sampling','t1d_ever','diabetes_at_sampling')

cag_abundances=read.csv('teddy_cags.csv',sep='\t')
cag_abundances=setNames(data.frame(t(cag_abundances[,-1])), cag_abundances[,1])
cag_abundances=cag_abundances %>% rownames_to_column()

full_dataset=merge(metadata,cag_abundances,by.x='sampleID',by.y='rowname') %>% filter(t1d_ever==1)

cagset=colnames(cag_abundances)

outputList_full=list()
outputList_subset=list()

for(cag in cagset){
  tryCatch({
    output=lmer(full_dataset[,cag]~diabetes_at_sampling+age+(1|subjectID),data=full_dataset)
    sum=summary(output)$optinfo$conv$lme4$messages
    if(is.null(sum)==TRUE){
      sum=''
    }
    output=tidy(output)
    output$CAG=cag
    output$message=sum
    outputList_full[[cag]]=output
  }, 
  error=function(e){
    print('Error')
  },
  warning=function(w){
    print('Warning')
  })  
} 

outputList_full=do.call(rbind,outputList_full)

write.csv(outputList_full,'teddy_association_outputs_pre_post_t1d.csv')
