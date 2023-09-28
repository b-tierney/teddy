#process teddy metadata

library(lme4)
library(tidyverse)
library(lmerTest)
library(broom.mixed)
options(warn=2)

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

metadata = read.csv('teddy_metadata_20190821.csv',stringsAsFactors = FALSE)

metadata = metadata %>% select(Run,maskid,age_at_collection,two_or_more_persistent,IA_Outcome,t1d,T1D_Outcome,GAD_pos,MIAA_pos,IA2A_pos,age_first_GAD,age_first_MIAA,age_first_IA2A)

metadata$GADA=0
metadata$GADA[metadata$age_at_collection>=metadata$age_first_GAD & metadata$age_first_GAD>0]=1
metadata$IAA=0
metadata$IAA[metadata$age_at_collection>=metadata$age_first_MIAA & metadata$age_first_MIAA>0]=1
metadata$IA2A=0
metadata$IA2A[metadata$age_at_collection>=metadata$age_first_IA2A & metadata$age_first_IA2A>0]=1

metadata$GAD_pos = as.factor(as.integer(as.logical(metadata$GAD_pos)))
metadata$MIAA_pos = as.factor(as.integer(as.logical(metadata$MIAA_pos)))
metadata$IA2A_pos = as.factor(as.integer(as.logical(metadata$IA2A_pos)))

metadata = metadata %>% select(-c(age_first_GAD,age_first_MIAA,age_first_IA2A))

colnames(metadata)=c('sampleID','subjectID','age','seroconverted_ever','seroconverted_at_sampling','t1d_ever','diabetes_at_sampling','GADA','IAA','IA2A','GADA_at_sampling','IAA_at_sampling','IA2A_at_sampling')

cag_abundances=read.csv('teddy_cags.csv',sep='\t')
cag_abundances=setNames(data.frame(t(cag_abundances[,-1])), cag_abundances[,1])
cag_abundances=cag_abundances %>% rownames_to_column()

full_dataset=merge(metadata,cag_abundances,by.x='sampleID',by.y='rowname')

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

write.csv(outputList_gada,'teddy_association_outputs_gada.csv')
write.csv(outputList_gada_subset,'teddy_association_outputs_gada_subset.csv')
write.csv(outputList_iaa,'teddy_association_outputs_iaa.csv')
write.csv(outputList_iaa_subset,'teddy_association_outputs_iaa_subset.csv')
write.csv(outputList_ia2a,'teddy_association_outputs_ia2a.csv')
write.csv(outputList_ia2a_subset,'teddy_association_outputs_ia2a_subset.csv')
