#merge teddy metadata in pandas

import pandas as pd

#krisch1=pd.read_csv('m_45_jkrischer_niddk_30apr2016_1.csv')
krisch=pd.read_csv('m_45_jkrischer_niddk_30apr2016_2.csv')
srarundata=pd.read_csv('SraRunTable.txt',sep='\t')
niddk_data=pd.read_csv('niddk_age_sample_mapping_data',sep='\t')

#merge the krischner data and take relevant columns
#krisch=krisch1.merge(krisch2,on='MASKID')
krisch=krisch[['VISIT_AGE_MOS','IA_ANY_LONG','IA_GAD_ONLY_LONG','IA_IAA_ONLY_LONG','IA_MULTIPLE_LONG','T1D_LONG','MASKID']]
krisch.MASKID=[int(x) for x in krisch.MASKID]

#filter sra columns and remove unnecessary data
srarundata=srarundata.loc[srarundata['LibrarySource']=='METAGENOMIC']
srarundata=srarundata[['Sample_Name','submitted_subject_id','sex']]

#merge the sra data and the krischner data
sra_krisch=krisch.merge(srarundata,left_on='MASKID',right_on='submitted_subject_id')

#solve the time from seroconversion problem 