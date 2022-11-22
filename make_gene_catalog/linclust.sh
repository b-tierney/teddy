#!/bin/bash
mmseqs/bin/mmseqs createdb $1 all_seqs_db_"$2"   

mmseqs/bin/mmseqs cluster -c 0.9 --min-seq-id 0.3 all_seqs_db_"$2" all_seqs_clu_30_"$2" tmpfolder_"$2"

mmseqs/bin/mmseqs createtsv all_seqs_db_"$2" all_seqs_db_"${2}" all_seqs_clu_30_"$2" all_seqs_clu_30_"$2".tsv

mmseqs/bin/mmseqs result2repseq all_seqs_db_"$2" all_seqs_clu_30_"$2" all_seqs_rep_30_"$2"

mmseqs/bin/mmseqs result2flat all_seqs_db_"$2" all_seqs_db_"$2" all_seqs_rep_30_"$2" all_seqs_rep_30_"$2".fasta --use-fasta-header
