import os
import pandas as pd
import sys

abundance_locs = sys.argv[1]
gene_name_locs = os.listdir('gene-names_batched')
suffix = abundance_locs.split('_')[0]

gene_names = [x for x in gene_name_locs if x[-2:]==suffix][0]
gene_names = 'gene-names_batched/'+gene_names

norm_data = []
norm_data.append(pd.read_csv(gene_names,header=None))

colnames = []
colnames.append('genename')
with open(abundance_locs) as f:
    for line in f:
        colnames.append(line.rstrip().split('/')[1].split('_')[0])
        norm_data.append(pd.read_csv(line.rstrip(),header=None))

norm_data_merged = pd.concat(norm_data,axis=1)
norm_data_merged.columns = colnames

norm_data_merged.to_csv('parsed_data/%s.csv'%suffix)
