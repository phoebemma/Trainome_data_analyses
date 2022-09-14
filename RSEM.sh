set -evx
input_dir=/datastore/Chidimma/Trainome_data/Training_data/STAR_outputs/*Aligned.toTranscriptome.out.bam
refseq=/datastore/Chidimma/Trainome_data/Genome_index/RSEM/Rsem_human_index
outdir=/datastore/Chidimma/Trainome_data/Training_data/RSEM_+_STAR_outputs
for file in $input_dir
do
rsem-calculate-expression --paired-end --alignments  --append-names -p 16 $file $refseq $outdir/${file##*/}
done

#THE --APPEND NAMES OPTION IS RECOMMENDED BY THE MAKERS
