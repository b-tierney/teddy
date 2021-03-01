while read p; do

sbatch -p short -c 5 --mem=80G -t 0-11:59 ./align_db.sh $p /n/scratch3/users/b/btt5/diabimmune/count_files /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/adapters.fa /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.bitmask /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.srprism /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/tmp 5 /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs_db.dmnd /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs.fasta

done<all

while read p; do

sbatch -p short -c 4 --mem=50G -t 0-11:59 ./align_db.sh $p /n/scratch3/users/b/btt5/diabimmune/count_files /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/adapters.fa /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.bitmask /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.srprism /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/tmp 4 /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs_db.dmnd /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs.fasta

done<subsetab

while read p; do

sbatch -p short -c 4 --mem=50G -t 0-8:59 ./align_db.sh $p /n/scratch3/users/b/btt5/diabimmune/count_files /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/adapters.fa /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.bitmask /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.srprism /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/tmp 4 /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs_db.dmnd /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs.fasta

done<subsetac

while read p; do

sbatch -p short -c 4 --mem=50G -t 0-8:59 ./align_db.sh $p /n/scratch3/users/b/btt5/diabimmune/count_files /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/adapters.fa /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.bitmask /home/sez10/kostic_lab/gene_catalogue/files_for_gene_catalogue/GRCh38/GRCh38.primary_assembly.genome.srprism /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/tmp 4 /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs_db.dmnd /n/scratch3/users/a/adk9/TEDDY/all_seqs_rep_30_db_ted_merged_seqs.fasta

done<subsetad