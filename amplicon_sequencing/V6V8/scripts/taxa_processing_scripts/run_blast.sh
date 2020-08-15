#!bin/bash/

#Kristen R. Hunter-Cevera, Marine Biological Laboratory, Spring 2020

# bash command to make reference database
# makeblastdb -in organism_reference.fasta -input_type fasta -dbtype nucl -out referencedatabase -title "referencedatabase"

# bash command for blast

database="/home/kristen/Documents/V6V8_analysis/marine_syn_16S_db/syn_blast_database/Syn_full_length_16S_db_plus_new_long_mvco_syn.fasta"

#query_file="/home/kristen/Documents/V6V8_analysis/data/vamps_downloads/fasta_files/mvco_env_seqs.fasta"
#query_file="/home/kristen/Documents/V6V8_analysis/analysis_products/mvco_timeseries_uniq.fasta"
query_file=$1
outfile=$2
echo "running blast on $query_file"
echo "saving to $outfile"

#blastn -db ${database} -query ${query_file} -out ${outfile} -outfmt 6 -num_threads 7 -max_target_seqs 1 -evalue 0.0000000001
