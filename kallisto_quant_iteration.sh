input_dir=/datastore/Chidimma/Trainome_data/Training_data/Trimmed_seqs
dir=$input_dir/*1P.fastq.gz
kallisto_index=/datastore/Chidimma/Trainome_data/Genome_index/GRCh38_kallisto_index
for file in $dir
do
file2=${file/_1P/_2P}
#Create a folder in the input directory with which to identify the input files
#I like using "1" in front of folder names to keep them on top of the directory

mkdir -p $input_dir/1kallisto_output/1A${file##*/}_dir

#get the newly created folder path and assign it to "Results"
Results=$(readlink -e $input_dir/1kallisto_output/1A${file##*/}_dir)

#copy the paired sequences to the folder. These would serve as markers to identify the sample ID of the "abundance.tsv" files 
#cp $file $Results
#cp $file2 $Results

#Run kallisto quantification
kallisto quant -i $kallisto_index -o $Results $file $file2 -t 16 --plaintext

cd $Results
mv abundance.tsv ${file##*/}.tsv
done

