 $ samtools view sample_103_repeatAligned.toTranscriptome.out.bam | less

samtools view /path to Bam or sam_file/ 
#less allows us scroll through the file

samtools view -q 30 -c sample_103_repeatAligned.toTranscriptome.out.bam 
#counts "-c" the number of reads with quality score 30 or higher "-q 30"

The flags for samtools are described below
Flag	Description
1	read is mapped
2	read is mapped as part of a pair
4	read is unmapped
8	mate is unmapped
16	read reverse strand
32	mate reverse strand
64	first in pair
128	second in pair
256	not primary alignment
512	read fails platform/vendor quality checks
1024	read is PCR or optical duplicate

samtools view -f 512 -c sample_103_repeatAligned.toTranscriptome.out.bam 
#checks the number of reads that failed quality checks

