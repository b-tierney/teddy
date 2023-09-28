import pandas as pd
import sys

tsv_data=sys.argv[1]
genecat=sys.argv[2]

tsv_data_for_df=[]
with open(tsv_data) as f:
		for line in f:
				if 'CDS' in line:
					line=line.rstrip().split('\t')
					tsv_data_for_df.append([line[0],line[2],line[3]])

tsv_data=pd.DataFrame(tsv_data_for_df)
tsv_data_annotated = tsv_data[tsv_data.loc[:,2]!='']
tsv_data_annotated = tsv_data_annotated[tsv_data_annotated.groupby([2],sort=False)[1].transform(max)==tsv_data_annotated[1]].drop_duplicates(2).drop(1,axis=1).set_index(2)
tsv_data_annotated.to_csv('temp_tsv_annotation.csv')
tsv_data=tsv_data.set_index(0).drop(1,axis=1)
tsv_data.to_csv('temp_tsv_data.csv')

with open(genecat) as f:
	with open('%s_collapsed.tsv'%genecat.replace('.tsv',''),'w') as w:
		for line in f:
			line=line.rstrip().split('\t')
			rawgene=line[1]
			congene=line[0]
			annotation=tsv_data.loc[congene,2][0]
			if annotation=='':
				annotation_congene=congene
			else:
				annotation_congene=tsv_data_annotated.loc[annotation].values[0]
			output=[rawgene,congene,annotation,annotation_congene]
			w.write('\t'.join(output)+'\n')