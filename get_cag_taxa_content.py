###get taxonomic mapping of CAGs

import os
from collections import Counter
import pandas as pd
import sys

def load_taxa_data(ncbi_tax_id_file):
	taxaNameMap={}
	with open(ncbi_tax_id_file) as f:
	  for line in f:
		taxaNameMap[line.split('\t')[0]]=line.rstrip().split('\t')[1]
	return taxaNameMap

def run_diamond(input_file):
	os.system('./diamond blastx -d ./nr.dmnd -q %s -o %s_taxa_mapping.tsv --taxonmap ./prot.accession2taxid.gz --taxonnodes ./nodes.dmp --outfmt 102 -k 1'%(input_file,input_file))

def find_species_abundances(input_file,ncbi_tax_id_file):
	taxaNameMap=load_taxa_data(ncbi_tax_id_file)
	mapped_data=[]
	mapped_data.append('GENE/CAG ID\tNCBI TAXA ID\tE-VALUE\tNAME\n')
	with open("%s_taxa_mapping.tsv"%input_file) as f:
		for line in f:
			line=line.rstrip().split('\t')[:3]
			if line[1]=='0' or line[1]=='1' or line[1]=='2':
				line.append('No annotation')
				mapped_data.append('\t'.join(line)+'\n')
				continue
			try:
				line.append(taxaNameMap[line[1]])
				mapped_data.append('\t'.join(line)+'\n')
			except:
				continue

	with open('%s_taxa_mapping.tsv'%input_file,'w') as w:
		for line in mapped_data:
			w.write(line)

	mapped_data_counts=pd.DataFrame.from_dict(Counter([x.rstrip().split('\t')[-1] for x in mapped_data[1:]]),orient='index')
	mapped_data_counts.to_csv("%s_taxa_mapping_counts.tsv"%input_file,sep='\t')

if __name__ == '__main__':
	input_file=sys.argv[1]
	ncbi_tax_id_file=sys.argv[2]
	run_diamond(input_file)
	find_species_abundances(input_file,ncbi_tax_id_file)
