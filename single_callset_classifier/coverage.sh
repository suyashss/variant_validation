#!/bin/bash

bamlist=$1
position_file=$2

mkdir -p coverage_stats
cat $bamlist | while read file ;
do 
	fname=`basename $file`; 
	qsub -b y -m n -cwd -l h_vmem=2G -e coverr -o coverage_stats/$fname.coverage samtools bedcov $position_file $file ;
done
	
