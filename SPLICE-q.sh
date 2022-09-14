set -evx
input_dir=/datastore/Chidimma/Trainome_data/Phase_1_analyses/COPD_dataset/STAR_outputs/*sortedByCoord.out.bam
output_dir=/datastore/Chidimma/Trainome_data/Phase_1_analyses/COPD_dataset/SPLICEq_outputs
annotation=/datastore/Chidimma/Trainome_data/Genome_index/gencode.v40.primary_assembly.annotation.gtf
for file in $input_dir
do
SPLICE-q.py -b $file -g $annotation -p 16 -o $output_dir/${file##*/}.tsv
done
