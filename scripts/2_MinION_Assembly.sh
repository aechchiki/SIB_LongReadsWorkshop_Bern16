#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J Assembly_all
#BSUB -q priority 
#BSUB -n 4
#BSUB -R "rusage[mem=4096]"
#BSUB -M 4194304

## prepare directory

# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
# create working directory for MinION extraction
mkdir -p $work_dir/lambda_MinION/miniasm

## load canu and miniasm
module add UHTS/Assembler/canu/1.3;
module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;

## Canu assembly

canu -p canu_MinION -d $work_dir/lambda_MinION/asm_canu/ genomeSize=48k -nanopore-raw $work_dir/lambda_MinION/MinION_Lambda2D.fastq -useGrid=false -maxThreads=4 &> canu_minion_lambda.log

## miniasm

minimap -Sw5 -L100 -m0 -t4 $work_dir/lambda_MinION/MinION_Lambda2D.fastq $work_dir/lambda_MinION/MinION_Lambda2D.fastq | gzip -1 > $work_dir/lambda_MinION/miniasm/minion_overlaps.paf.gz

miniasm -f $work_dir/lambda_MinION/MinION_Lambda2D.fastq $work_dir/lambda_MinION/miniasm/minion_overlaps.paf.gz > $work_dir/lambda_MinION/miniasm/lambda_miniasm_contigs.gfa

awk '/^S/{print ">"$2"\n"$3}' $work_dir/lambda_MinION/miniasm/lambda_miniasm_contigs.gfa | fold > $work_dir/lambda_MinION/miniasm/miniasm_MinION_contigs.fa
