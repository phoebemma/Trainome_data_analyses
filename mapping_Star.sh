set -evx
input_dir=/datastore/Chidimma/Trainome_data/Training_data/Trimmed_seqs/*1P.fastq.gz
genome_dir=/datastore/Chidimma/Trainome_data/Genome_index/STAR 
annotated_dir=/datastore/Chidimma/Trainome_data/Genome_index/gencode.v40.primary_assembly.annotation.gtf 
output_dir=/datastore/Chidimma/Trainome_data/Training_data/STAR_outputs_for_stringtie
for file in $input_dir
do
file2=${file/_1P/_2P}
STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN 16 --sjdbGTFfile $annotated_dir --readFilesCommand zcat --winAnchorMultimapNmax 150 --outFilterMultimapNmax 80 --outReadsUnmapped Fastx --outMultimapperOrder Random --outWigType wiggle --outFilterMismatchNoverLmax 0.1 --genomeDir $genome_dir --outSAMstrandField intronMotif --twopassMode Basic --readFilesIn $file $file2 --outFileNamePrefix $output_dir/${file##*/};

done
#add --genomeLoad LoadAndKeep for one_pass mapping
#{file##*/} extracts the filename from a path. This can also be done using "basename /path to file/"

 #--twopassMode Basic 
