#get collapsed consensus sequences

import os
import sys
from Bio import SeqIO


#sequencefile=sys.argv[1]
sequencefile='all_seqs_rep_90_teddy.fasta'
congenetsv='all_seqs_clu_90_teddy_collapsed.tsv'

#os.system("""awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%%10000000==0){file=sprintf("tempseq%%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < %s"""%sequencefile)

#files=os.listdir('.')
#files=[x for x in files if 'tempseq' in x]

files=[sequencefile]
print('loading congene tsv')
congeneids=[]
with open(congenetsv) as f:
	for line in f:
		congeneids.append(line.rstrip().split('\t')[-1])

congeneids = set(congeneids)
print(len(congeneids))
for f in files:
	print('loading fasta file')
	record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
	with open('subsampled_seqs_'+f,'w') as w:
		for c in congeneids:
			try:
				record = record_dict[c]
				w.write('>'+c+'\n')
				w.write(str(record.seq)+'\n')
			except:
				continue

#os.system("cat subsampled_seqs* > %s_collapsed.fasta"%(sequencefile.split('.fasta')[0]))
#os.system("rm subsampled_seqs_*")
#os.system("rm tempseq_*")
