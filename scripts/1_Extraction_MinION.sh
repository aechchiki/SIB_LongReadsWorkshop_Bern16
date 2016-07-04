#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J Extraction_MinION
#BSUB -q priority 

## prepare directory
# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
mkdir -p $work_dir
# create working directory for MinION extraction
mkdir -p $work_dir/lambda_minion/fast5/
# create links to the raw sequencing files
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/MinION_lambda_reads/*.fast5 $work_dir/lambda_minion/fast5/

## load poretools for read extraction
module add UHTS/Analysis/poretools/0.5.1 

## extract reads 
# convert all fast5 reads to fastq
poretools fastq --type 2D $work_dir/lambda_minion/fast5/*.fast5 > $work_dir/lambda_minion/all_reads.fastq
# get some stats on the raw data using poretools
poretools stats $work_dir/lambda_minion/fast5/*.fast5 > $work_dir/lambda_minion/minion_all_reads_stats.txt

# last tested on: Mon 04 Jul 2016 11:33:07 PM CEST 
# AE
