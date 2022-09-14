set -evx

directory=/datastore/Chidimma/Trainome_data/Phase_1_analyses/COPD_dataset/SplAdder_workspace/vitd_group
input=$directory/*.bam 
output_dir=$directory/1SplAdder
cd $directory

spladder test -o $output_dir --conditionA 4-102PreExcVLR12_S4_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,37-113PreExcVLR74_S48_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 5-102PostExcVLR13_S5_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,39-113PostExcVLR76_S11_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA Pre_exercise_Right_10 --labelB Post_exercise_Right_10 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12

spladder test -o $output_dir --conditionA 4-102PreExcVLR12_S4_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,37-113PreExcVLR74_S48_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 1-102PreSuppVLR9_S1_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,34-113PreSuppVLR71_S30_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA Pre_exercise_Right_10 --labelB Pre_supplementation_Right_10 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12

spladder test -o $output_dir --conditionA 4-102PreExcVLR12_S4_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,37-113PreExcVLR74_S48_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 2-1023WVLR10_S2_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,40-1133WVLR77_S18_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA Pre_exercise_Right_10 --labelB Week_3_Right_10 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12

spladder test -o $output_dir --conditionA 1-102PreSuppVLR9_S1_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,34-113PreSuppVLR71_S30_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 2-1023WVLR10_S2_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,40-1133WVLR77_S18_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA Pre_supplemenatation_Right_10 --labelB Week_3_Right_10 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12

spladder test -o $output_dir --conditionA 3-1023WVLL11_S3_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,36-1133WVLL73_S42_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 6-102PostExcVLL14_S6_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,38-113PostExcVLL75_S54_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA week_3_left_30 --labelB Postexercise_left_30 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12

spladder test -o $output_dir --conditionA 2-1023WVLR10_S2_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,40-1133WVLR77_S18_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 5-102PostExcVLR13_S5_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,39-113PostExcVLR76_S11_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA Week_3_Right_10 --labelB Post_exercise_Right_10 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12


spladder test -o $output_dir --conditionA 6-102PostExcVLL14_S6_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,38-113PostExcVLL75_S54_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --conditionB 5-102PostExcVLR13_S5_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam,39-113PostExcVLR76_S11_L001_R1_001_1P.fastq.gzAligned.sortedByCoord.out.bam --labelA Post_Exercise_Left_30 --labelB Post_exercise_Right_10 --merge-strat merge_graphs --confidence 2 --diagnose-plots --plot-format pdf --cap-outliers --parallel 12

