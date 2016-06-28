#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 1_Extraction_MinION.err
#BSUB -o 1_Extraction_MinION.out
#BSUB -J Extraction_MinION
#BSUB -q priority 

## prepare directory

# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
# create working directory for MinION extraction
mkdir -p $work_dir/lambda_MinION/fast5/
# create links to the raw sequencing files
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/MinION_lambda_reads/*.fast5 $work_dir/lambda_MinION/fast5/

## load poretools for read extraction
module add UHTS/Analysis/poretools/0.5.1 

## extract reads 

# convert all fast5 reads to fastq
poretools fastq --type 2D $work_dir/lambda_MinION/fast5/*.fast5 > $work_dir/lambda_MinION/MinION_Lambda2D.fastq

# get some stats on the raw data using poretools
poretools stats $work_dir/lambda_MinION/fast5/*.fast5 > $work_dir/lambda_MinION/MinION_Lambda2D_stats.txt

1_Extraction_MinION.sh (END)
