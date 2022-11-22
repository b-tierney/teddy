

import os
import sys
from Bio import SeqIO


sequencefile=sys.argv[1]
congenetsv='missing'

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

