set -evx

directory=/datastore/Chidimma/Trainome_data/Training_data/STAR_outputs
input=$directory/*sortedByCoord.out.bam
annotated_dir=/datastore/Chidimma/Trainome_data/Genome_index/gencode.v40.primary_assembly.annotation.gtf 
output_dir=$directory/1Spladder
cd $directory

ls -1 $input > bamfiles.txt
#EACH COMMENDED ONE IS A STEP RUN BEFORE THE ONE THAT FOLLOWS IT

spladder build -o $output_dir -a $annotated_dir -b bamfiles.txt --merge-strat single --no-extract-ase --sparse-bam --parallel 24 --confidence 2 --output-txt


spladder build -o $output_dir -a $annotated_dir -b bamfiles.txt --merge-strat merge_graphs --no-extract-ase --sparse-bam --parallel 24 --confidence 2 --output-txt --quantify-graph

for file in $input
do
spladder build -o $output_dir -b $file -a $annotated_dir --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode single --confidence 2 --parallel 24
done 

spladder build -o $output_dir -a $annotated_dir -b bamfiles.txt --merge-strat merge_graphs --no-extract-ase --sparse-bam --parallel 24 --confidence 2 --output-txt --quantify-graph --qmode collect

spladder build -o $output_dir -a $annotated_dir -b bamfiles.txt --merge-strat merge_graphs --ase-edge-limit 1500 --sparse-bam --parallel 24 --confidence 2 --output-txt --no-compress-text

