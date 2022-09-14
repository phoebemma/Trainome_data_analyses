set -evx
input_dir=/datastore/Stian_RNA_Seq/rawfiles_training_data/rawfiles/*R1_001.fastq.gz;
trim_path=/datastore/Chidimma/tools/Trimmomatic-0.39
output_dir=/datastore/Chidimma/Trainome_data/Training_data/Trimmed_seqs

for f in $input_dir
do
java -jar $trim_path/trimmomatic-0.39.jar PE -threads 24 -basein $f -baseout $output_dir/${f##*/} ILLUMINACLIP:$trim_path/adapters/TruSeq3-PE-2.fa:2:30:15:8:true MAXINFO:40:0.7 MINLEN:10 LEADING:3 TRAILING:3
done

#${f##*/}

