library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)
library(broom)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/22417/for_validation/')

metadata=read.csv('teddy_metadata_20190821.csv',stringsAsFactors = FALSE)

#subset metadata to relevant samples
metadata_presero=metadata[metadata$IA_Outcome %in% c('Before','Never'),]
metadata_presero$IA_Outcome[metadata_presero$IA_Outcome == 'Before'] = 'Seroconverter'
metadata_presero$IA_Outcome[metadata_presero$IA_Outcome == 'Never'] = 'Healthy'

#load abundance data
relative_abundances=read.csv('cag_relative_abundance_for_validation.tsv',sep='\t',row.names=1)
relative_abundances=t(relative_abundances)
combined=merge(metadata_presero,relative_abundances,by.x='Run',by.y='row.names')

#perform t test over averaged sample data
output=list()
cags=colnames(combined[135:ncol(combined)])
tosubset=c(c('IA_Outcome','dbgap_maskid'),cags)
combined=combined[,tosubset]
combined=aggregate(combined[,3:ncol(combined)],list(combined$IA_Outcome,combined$dbgap_maskid),mean)
colnames(combined)[1]='IA_Outcome'

count=0
for(c in cags){
  count=count+1
  x=combined[,c][which(combined$IA_Outcome=='Healthy')]
  y=combined[,c][which(combined$IA_Outcome=='Seroconverter')]
  output[[count]]=tidy(t.test(x,y))
  ggplot(combined,aes(x=as.factor(combined$IA_Outcome),y=combined[,c]))+geom_boxplot()+geom_jitter()+ggtitle(paste('Validation for',c))+ylab('ln(Relative Abundance)')+xlab('')
  ggsave(paste('./validation_plots/',c,'.pdf',sep=''),width = 5,height=5,units = 'in')
}

out=ldply(output,rbind)
out$FDR=p.adjust(out$p.value)
out$CAGID=cags
write.csv(out,'validate_t_tests_after_sero.csv')
