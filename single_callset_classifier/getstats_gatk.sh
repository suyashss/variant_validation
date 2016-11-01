testsample_bam=$1
outfile=$2
bedfile=$3


basedir="~/bustamante/progs/gatk"
java -Xmx6G -jar /home/suyashs/bustamante/progs/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
-R ~/bustamante/genomes/genomes/hg19/gatk_resources/ucsc.hg19.fasta \
-T UnifiedGenotyper \
-I $testsample_bam \
-o $outfile \
--output_mode EMIT_ALL_SITES \
-L $bedfile \
-A AlleleBalanceBySample -A FisherStrand -A GCContent -A SpanningDeletions
