input_dir=/datastore/Chidimma/rawfiles_COPD_data/Trimmed_reads/sample_103/*_1P.fastq.gz
genome_dir=/datastore/Chidimma/Genome/STAR 
annotated_dir=/datastore/Chidimma/Genome/gencode.v40.annotation.gtf 
output_dir=/datastore/Chidimma/Genome/Annotation_output/
for file in $input_dir
file2=${file/_1P/_2P}
do
STAR --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN 16 --sjdbGTFfile $annotated_dir --readFilesCommand zcat --outReadsUnmapped Fastx --outMultimapperOrder Random --outWigType wiggle --genomeDir $genome_dir --readFilesIn $file $file2 --outFileNamePrefix $output_dir/${file##*/};

#STAR  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN 16 --sjdbGTFfile /datastore/Chidimma/Genome/gencode.v40.primary_assembly.annotation.gtf --readFilesCommand zcat --outReadsUnmapped Fastx --outMultimapperOrder Random --outWigType wiggle --genomeDir /datastore/Chidimma/Genome/STAR --readFilesIn /datastore/Chidimma/rawfiles_COPD_data/Trimmed_reads/sample_103/8-103PreSuppVLR21_S14_L001_R1_001trimmed_2P.fastq.gz --outFileNamePrefix /datastore/Chidimma/Genome/star_output/sample_103_PreSupp
done

#{file##*/} extracts the filename from a path. This can also be done using "basename /path to file/"
