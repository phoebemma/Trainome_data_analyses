set -evx
#
#I RAN THIS FIRST TO  GET THE LIST OF GENOME TARGETS
grep "^>" < GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt


#Concatenate the transcript and primary assemply files
cat gencode.v40.transcripts.fa GRCh38.primary_assembly.genome.fa > gentrome.fa

/datastore/Chidimma/salmon-1.9.0_linux_x86_64/bin/salmon index -t salmon_gentrome.fa -d decoys.txt -p 12 -i salmon_index
