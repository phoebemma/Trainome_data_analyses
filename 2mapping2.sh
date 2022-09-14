set -evx

input_dir=/datastore/Chidimma/old_vs_young/young/*_1P.fastq.gz
genome_dir=/datastore/Chidimma/Trainome_data/Genome_index/STAR 
annotated_dir=/datastore/Chidimma/Trainome_data/Genome_index/gencode.v40.primary_assembly.annotation.gtf 
output_dir=/datastore/Chidimma/old_vs_young/Star_outputs/young
#one_pass STAR mapping
#for file in $input_dir
#do
#file2=${file/_1P/_2P}
#STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN 16 --sjdbGTFfile $annotated_dir --readFilesCommand zcat --winAnchorMultimapNmax 150 --outFilterMultimapNmax 80 --outReadsUnmapped Fastx --outMultimapperOrder Random --outWigType wiggle --outFilterMismatchNoverLmax 0.1 --genomeDir $genome_dir --readFilesIn $file $file2 --outFileNamePrefix $output_dir/${file##*/};
#done
#echo "Mapping completed"

#Indexing using samtools. 
for sorted_file in $output_dir/*sortedByCoord.out.bam
do
samtools index $sorted_file
done

echo "Indexing completed"
 
#counting using kallisto


#echo "Transcript quantification completed"

#Splicing analyses using Spladder
#mkdir /datastore/Chidimma/old_vs_young/Spladder_outputs

Spladder_outputs=/datastore/Chidimma/old_vs_young/Spladder_outputs
for sorted_file in $output_dir/*sortedByCoord.out.bam
do

#spladder build -o $Spladder_outputs -a $annotated_dir -b $sorted_file --merge-strat single --no-extract-ase --ase-edge-limit 1500 --sparse-bam --parallel 24
#done

#spladder build -o $Spladder_outputs -a $annotated_dir -b /datastore/Chidimma/old_vs_young/young_bam.txt --merge-strat merge_graphs --no-extract-ase --reference $genome_dir --sparse-bam --parallel 16
spladder build -o $Spladder_outputs -a $annotated_dir -b $sorted_file --merge-strat merge_graphs --no-extract-ase --ase-edge-limit 1500 --sparse-bam --quantify-graph --qmode single --parallel 24
done
