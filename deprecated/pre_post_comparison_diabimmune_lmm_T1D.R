#pre and post seroconversion lmms

library(lme4)
library(tidyverse)
library(plyr)
library(lmerTest)
library(broom.mixed)
options(warn=2)
###run mixed regression for each CAG

metadata=read.csv('diabimmune_metadata.tsv',sep='\t')
cag_abundances=t(read.csv('diabimmune_cags.tsv',sep='\t'))
metadata=metadata %>% transform(t1d_ever=revalue(t1d_ever,c("TRUE"=1))) 

full_dataset=merge(metadata,cag_abundances,by.x='sampleID',by.y=0) %>% filter(t1d_ever==1)

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
    print(e)
  },
  warning=function(w){
    print('Warning')
  })  
} 

outputList_full=do.call(rbind,outputList_full)

write.csv(outputList_full,'diabimmune_association_outputs_pre_post_t1d.csv')
