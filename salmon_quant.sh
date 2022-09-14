set -evx
input_dir=/datastore/Chidimma/Trainome_data/Training_data/Trimmed_seqs/*1P.fastq.gz
salmon_index=/datastore/Chidimma/Trainome_data/Genome_index/salmon_index
output_dir=/datastore/Chidimma/Trainome_data/Training_data/Salmon_output

for file in $input_dir
do
file2=${file/_1P/_2P}
/datastore/Chidimma/salmon-1.9.0_linux_x86_64/bin/salmon quant -i $salmon_index -l IU --writeUnmappedNames -1 $file -2 $file2 -z -p 16 -o $output_dir/${file##*/}_dir

cd $output_dir/${file##*/}_dir
mv quant.sf ${file##*/}.sf

cp ${file##*/}.sf /datastore/Chidimma/Pipeline/Salmon_output/1Salmon_output
done
#Passing the --writeUnmappedNames flag to Salmon will tell Salmon to write out the names of reads (or mates in paired-end reads) that do not map to the transcriptome.
