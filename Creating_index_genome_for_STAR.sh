DIR=/datastore/Chidimma/Trainome_data/Genome_index

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $DIR/STAR --genomeFastaFiles $DIR/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $DIR/gencode.v40.primary_assembly.annotation.gtf --sjdbOverhang 150 



