library(lme4)
library(tidyverse)
library(lmerTest)
library(plyr)
library(broom.mixed)
options(warn=2)
###run mixed regression for each CAG

run_regression <- function(cag,antibody,outList,dataset){
  tryCatch({
    output=lmer(dataset[,cag]~dataset[,antibody]+age+(1|subjectID),data=dataset)
    sum=summary(output)$optinfo$conv$lme4$messages
    if(is.null(sum)==TRUE){
      sum=''
    }
    output=tidy(output)
    output$CAG=cag
    output$message=sum
    outList[[cag]]=output
  }, 
  error=function(e){
    print('Error')
  },
  warning=function(w){
    print('Warning')
  })
  return(outList)
}

metadata=read.csv('diabimmune_metadata_with_aa_data.csv',sep=',')
cag_abundances=t(read.csv('diabimmune_cags.tsv',sep='\t'))

full_dataset=merge(metadata,cag_abundances,by.x='sampleID',by.y=0)

gada_subset=full_dataset %>% filter(GADA_at_sampling!=1)
iaa_subset=full_dataset %>% filter(IAA_at_sampling!=1)
ia2a_subset=full_dataset %>% filter(IA2A_at_sampling!=1)

cagset=colnames(cag_abundances)

outputList_gada=list()
outputList_gada_subset=list()
outputList_iaa=list()
outputList_iaa_subset=list()
outputList_ia2a=list()
outputList_ia2a_subset=list()

for(cag in cagset){
  outputList_gada=run_regression(cag,'GADA',outputList_gada,full_dataset)
  outputList_gada_subset=run_regression(cag,'GADA_at_sampling',outputList_gada_subset,gada_subset)
  outputList_iaa=run_regression(cag,'IAA',outputList_iaa,full_dataset)
  outputList_iaa_subset=run_regression(cag,'IAA_at_sampling',outputList_iaa_subset,iaa_subset)
  outputList_ia2a=run_regression(cag,'IA2A',outputList_ia2a,full_dataset)
  outputList_ia2a_subset=run_regression(cag,'IA2A_at_sampling',outputList_ia2a_subset,ia2a_subset)
} 

outputList_gada=do.call(rbind,outputList_gada)
outputList_gada_subset=do.call(rbind,outputList_gada_subset)
outputList_iaa=do.call(rbind,outputList_iaa)
outputList_iaa_subset=do.call(rbind,outputList_iaa_subset)
outputList_ia2a=do.call(rbind,outputList_ia2a)
outputList_ia2a_subset=do.call(rbind,outputList_ia2a_subset)

write.csv(outputList_gada,'diabimmune_association_outputs_gada.csv')
write.csv(outputList_gada_subset,'diabimmune_association_outputs_gada_subset.csv')
write.csv(outputList_iaa,'diabimmune_association_outputs_iaa.csv')
write.csv(outputList_iaa_subset,'diabimmune_association_outputs_iaa_subset.csv')
write.csv(outputList_ia2a,'diabimmune_association_outputs_ia2a.csv')
write.csv(outputList_ia2a_subset,'diabimmune_association_outputs_ia2a_subset.csv')




