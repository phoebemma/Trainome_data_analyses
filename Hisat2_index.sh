
#Directions from http://daehwankimlab.github.io/hisat2/howto/
#Extract splice sites into a file called "Hisat2_genome.ss"
hisat2_extract_splice_sites.py /datastore/Chidimma/Trainome_data/Genome_index/gencode.v40.primary_assembly.annotation.gtf > /datastore/Chidimma/Trainome_data/Genome_index/genome.ss

#Extract exons into a file called "Hisat2_genome.exon
hisat2_extract_exons.py /datastore/Chidimma/Trainome_data/Genome_index/gencode.v40.primary_assembly.annotation.gtf > /datastore/Chidimma/Trainome_data/Genome_index/Hisat2_genome.exon

#Building the index
hisat2-build -p 8 --exon /datastore/Chidimma/Trainome_data/Genome_index/HISAT2_index/Hisat2_genome.exon --ss /datastore/Chidimma/Trainome_data/Genome_index/HISAT2_index/Hisat2_genome.ss /datastore/Chidimma/Trainome_data/Genome_index/GRCh38.primary_assembly.genome.fa /datastore/Chidimma/Trainome_data/Genome_index/HISAT2_index/Human_genome

