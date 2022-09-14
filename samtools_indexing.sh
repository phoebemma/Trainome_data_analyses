set -evx
directory=/datastore/Chidimma/Trainome_data/Training_data/STAR_outputs/*sortedByCoord.out.bam

for file in $directory
do
samtools index $file
done

