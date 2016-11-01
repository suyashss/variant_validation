hwe_illumina_vcf=../../2016_08_01/hwe_illumina.vcf.gz
hwe_rtg_vcf=../../2014_07_29/hwe_rtg.vcf.gz
hwe_gatk_vcf=../../2016_08_01/hwe_gatk.vcf.gz
plink2=~/progs/plink_4_26_14/plink
rtgjar=~/bustamante/progs/rtg-core-3.3.1/RTG.jar 
vcftobed=~/suyashs_newprojects/from_gs1_to_organize/scripts/vcf_to_bed.sh 
illuminapval=1e-6
rtgpval=1e-6
gatkpval=1e-6
sampleid=21-999-056_MAXGT
#indelfile=$(sampleid).indels.pos # File with start position of indels, one location per line
indelfile=../../../data/Indels/21-999-056_MAXGT.pos
#bamfile=$(sampleid).bam
bamfile=/home/suyashs/suyashs_newprojects/caapa_bams/LP6005045-DNA_E04.gatk.bam
bamlist=../../2016_08_01/gatkbams.list 
omnivcf=../../../data/chr22_Omni.vcf.gz 

all:
	echo "Dummy target. Targets are listed below"
	grep '^[a-z]*:' makefile
	echo "make getcoverage" 
	echo "Wait for clusterjobs to finish"
	echo "make allroc_combined.png"

.SECONDARY:

# Generate VCF files for the test sample for each callset
testsample_%hwe.recode.vcf: 
	vcftools --gzvcf $(hwe_$*_vcf) --remove-indels --indv $(sampleid) --non-ref-ac 1 --out testsample_$*hwe --min-alleles 2 --max-alleles 2 --recode

# Convert sample VCFs to bed files for variant positions
testsample_%.bed: testsample_%hwe.recode.vcf
	$(vcftobed) testsample_$*hwe.recode.vcf testsample_$*.bed

# Create union bed file
testsample_illumina_rtg_gatk.bed: testsample_illumina.bed testsample_rtg.bed testsample_gatk.bed
	bedtools multiinter -i $^ | sed 's/^22/chr22/g' > $@

# Create list of genomic positions in the union
testsample_illumina_rtg_gatk.pos: testsample_illumina_rtg_gatk.bed
	cut -f-3 $< > $@

# Get annotations at the union positions from GATK emit-all mode
testsample_union_stats.raw.vcf: testsample_illumina_rtg_gatk.bed
	./getstats_gatk.sh $(bamfile) $@ $<

# Parse GATK annotations to get allele balance from GT field
testsample_union.ab: testsample_union_stats.raw.vcf
	python parse_anno.py $< > $@

# Get other site-specific annotations from the INFO field
testsample_union.sitefeatures: testsample_union_stats.raw.vcf
	python parse_anno_sitefeatures.py $< > $@

# Find coverage for all samples at the union positions
# These jobs run on a cluster and dependent jobs must be run after these jobs are done
getcoverage: $(bamlist) testsample_illumina_rtg_gatk.pos
	./coverage.sh $< $(word 2, $^)
	@touch $@

# Obtain total coverage by summing per-sample coverages
testsample_union.coverage: getcoverage	
	Rscript overallcoverage.r $@

#frequencies: Used to get site frequencies for illumina,rtg,gatk
%.frq:	
	vcftools --gzvcf $(hwe_$*_vcf) --freq2 --out $* --remove-indv $(sampleid) 

# Get Indel proximity annotations (1 if there is an indel with 10 bp of SNP, 0 otherwise)
testsample_union.indelprox: testsample_illumina_rtg_gatk.pos $(indelfile)
	Rscript indel_proximity.r  $< $(word 2, $^) $@

# Combine callset membership information, allele balance, coverage information, indel proximity, simite features and frequency information
testsample_union.trainingdata.noNA: testsample_illumina_rtg_gatk.bed testsample_union.ab testsample_union.coverage testsample_union.indelprox testsample_union.sitefeatures illumina.frq rtg.frq gatk.frq
	Rscript maketrainingdata.r $^ $(basename $@).withNA $@

testsample_omni.recode.vcf:
	vcftools --gzvcf $(omnivcf) --remove-indels --indv $(sampleid) --non-ref-ac 1 --max-alleles 2 --recode --out testsample_omni

testsample_omni_notvar.recode.vcf:
	vcftools --gzvcf $(omnivcf) --remove-indels --indv $(sampleid) --max-non-ref-ac 0 --max-alleles 2 --recode --out testsample_omni_notvar

testsample_omni.bed: testsample_omni.recode.vcf
	$(vcftobed) $< $@

testsample_omni_notvar.bed: testsample_omni_notvar.recode.vcf
	$(vcftobed) $< $@

%_rocinfo.txt: testsample_union.trainingdata.noNA testsample_omni.bed testsample_omni_notvar.bed
	Rscript trainrf_smote.r $* $< $(word 2, $^) $(word 3, $^)

allroc_combined.png: illumina_rocinfo.txt rtg_rocinfo.txt gatk_rocinfo.txt
	Rscript plot_allrocs.r $@
