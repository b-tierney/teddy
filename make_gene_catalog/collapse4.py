
import pandas as pd
import sys
import math

genecat='all_seqs_clu_90_teddy.tsv'

tsv_data_annotated=pd.read_csv('temp_tsv_annotation.csv',index_col=0)
tsv_data_annotated=tsv_data_annotated.to_dict()

tsv_data=pd.read_csv('temp_tsv_data.csv',index_col=0)
tsv_data=tsv_data.to_dict()

with open(genecat) as f:
        with open('%s_collapsed.tsv'%genecat.replace('.tsv',''),'w') as w:
                for line in f:
                        line=line.rstrip().split('\t')
                        rawgene=line[1]
                        congene=line[0]
                        annotation=tsv_data['2'][congene]
                        if str(annotation)=='nan':
                                annotation_congene=congene
                                annotation='None'
                        else:
                                annotation_congene=tsv_data_annotated['0'][annotation]
                        output=[rawgene,congene,annotation,annotation_congene]
                        w.write('\t'.join(output)+'\n')
                                                           
