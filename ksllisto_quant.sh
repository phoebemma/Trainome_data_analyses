set -evx

input_dir=/datastore/Chidimma/Trainome_data/Training_data/Trimmed_seqs/F2-FP12w2preR_S6_L001_R1_001_1P.fastq.gzF
kallisto_index=/datastore/Chidimma/Trainome_data/Genome_index/GRCh38_kallisto_index
dir=$input_dir/*1P.fastq.gz
for file in $dir
do
file2=${file/_1P/_2P}
kallisto quant -i $kallisto_index -o $input_dir $file $file2 -t 16 --plaintext
done


