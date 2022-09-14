set -a
#input_dir=/datastore/Chidimma/Trainome_data/COPD_data/190815_J00146./subject_151
#files=$input_dir//sample_151_trimmed_seqs/*P.fastq.gz

#cd $input_dir
#mkdir subject_151_fastqc_posttrimming

for file in $files
do
 fastqc $file --outdir ./results --threads 16
done
