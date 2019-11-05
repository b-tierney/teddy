library(lme4)
library(tidyverse)
library(lmerTest)
library(plyr)
library(broom.mixed)
options(warn=2)
###run mixed regression for each CAG

metadata=read.csv('diabimmune_metadata.tsv',sep='\t')
cag_abundances=t(read.csv('diabimmune_cags.tsv',sep='\t'))
metadata=metadata %>% transform(t1d_ever=revalue(t1d_ever,c("TRUE"=1))) 


full_dataset=merge(metadata,cag_abundances,by.x='sampleID',by.y=0)
full_dataset2=full_dataset %>% filter(diabetes_at_sampling!=1)

cagset=colnames(cag_abundances)

outputList_ever_full=list()
outputList_sampling_full=list()
outputList_subset=list()
for(cag in cagset){
  tryCatch({
    output=lmer(full_dataset[,cag]~t1d_ever+age+(1|subjectID),data=full_dataset)
    sum=summary(output)$optinfo$conv$lme4$messages
    if(is.null(sum)==TRUE){
      sum=''
    }
    output=tidy(output)
    output$CAG=cag
    output$message=sum
    outputList_ever_full[[cag]]=output
  }, 
  error=function(e){
    print('Error')
  },
  warning=function(w){
    print('Warning')
  })  
  tryCatch({
    output=lmer(full_dataset[,cag]~diabetes_at_sampling+age+(1|subjectID),data=full_dataset)
    sum=summary(output)$optinfo$conv$lme4$messages
    if(is.null(sum)==TRUE){
      sum=''
    }
    output=tidy(output)
    output$CAG=cag
    output$message=sum
    outputList_sampling_full[[cag]]=output
  }, 
  error=function(e){
    print('Error')
  },
  warning=function(w){
    print('Warning')
  })  
  tryCatch({
    output=lmer(full_dataset2[,cag]~t1d_ever+age+(1|subjectID),data=full_dataset2)
    sum=summary(output)$optinfo$conv$lme4$messages
    if(is.null(sum)==TRUE){
      sum=''
    }
    output=tidy(output)
    output$CAG=cag
    output$message=sum
    outputList_subset[[cag]]=output
  }, 
  error=function(e){
    print('Error')
  },
  warning=function(w){
    print('Warning')
  })  
} 

outputList_ever_full=do.call(rbind,outputList_ever_full)
outputList_sampling_full=do.call(rbind,outputList_sampling_full)
outputList_subset=do.call(rbind,outputList_subset)

write.csv(outputList_ever_full,'diabimmune_association_outputs_t1d_ever_full.csv')
write.csv(outputList_sampling_full,'diabimmune_association_outputs_t1d_sampling_full.csv')
write.csv(outputList_subset,'diabimmune_association_outputs_pret1d.csv')
