import pandas as pd
import os
import sys

locationfile = sys.argv[1]

with open(locationfile) as f:
	for line in f:
		location = line.rstrip()
		name = location.split('/')[-1][:-7].split('_')[0]
		outputloc = 'normalized_data/' + location.split('/')[-1][:-7] + '_normalized.tsv'
		print('Loading data for %s'%location)
		os.system('pigz -d %s'%location)
		data = pd.read_csv(location[:-3],sep='\t',index_col=0,header=None)
		os.system('pigz %s'%location[:-3])
		print('	Normalizing data')
		data[4] = data[2]/data[1]
		output = data[4]/sum(data[4])
		output.name=name
		print('	Writing to file')
		output.to_csv(outputloc,sep='\t',index=False,columns=[name])
		print('	Done')
