#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J RSII_asm
#BSUB -e RSII-assembly.err
#BSUB -o RSII-assembly.out
#BSUB -q dee-hugemem
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "rusage[swp=20000]"
#BSUB -v 20000000
#BSUB -M 8000000

## prepare directory
# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
# create working directory for MinION extraction
mkdir -p $work_dir/lambda_RSII/asm_miniasm

module add UHTS/Assembler/canu/1.3;

# take a subset of reads
head -12000 $work_dir/RSII_LambdaSubreads.fastq > $work_dir/RSII_reads_subset.fastq

# canu assemblu
canu -p canu_RSII -d $work_dir/lambda_RSII/asm_canu/ genomeSize=48k -pacbio-raw $work_dir/lambda_RSII/RSII_reads_subset.fastq -useGrid=false -maxThreads=4

cp $work_dir/lambda_RSII/asm_canu/canu_RSII.contigs.fasta $work_dir/lambda_RSII/canu_RSII_contigs.fa

## miniasm

module purge
module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;

minimap -Sw5 -L100 -m0 -t4 $work_dir/lambda_RSII/RSII_reads_subset.fastq $work_dir/lambda_RSII/RSII_reads_subset.fastq | gzip -1 > $work_dir/lambda_RSII/asm_miniasm/RSII_overlaps.paf.gz

miniasm -f $work_dir/lambda_RSII/RSII_reads_subset.fastq $work_dir/lambda_RSII/asm_miniasm/RSII_overlaps.paf.gz > $work_dir/lambda_RSII/asm_miniasm/lambda_miniasm_contigs.gfa

awk '/^S/{print ">"$2"\n"$3}' $work_dir/lambda_RSII/asm_miniasm/lambda_miniasm_contigs.gfa | fold > $work_dir/lambda_RSII/miniasm_RSII_contigs.fa


