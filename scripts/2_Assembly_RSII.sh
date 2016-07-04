#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J RSII_assembly
#BSUB -e RSII-assembly.err
#BSUB -o RSII-assembly.out
#BSUB -q dee-hugemem
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "rusage[swp=20000]"
#BSUB -v 20000000
#BSUB -M 8000000

### canu

## prepare directory
# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
# prepare assembly directory
mkdir -p $work_dir/lambda_RSII/assembly_canu/

## pre-assembly
# take a subset of reads
head -12000 $work_dir/lambda_RSII/RSII_LambdaSubreads.fastq > $work_dir/RSII_reads_subset.fastq

## assembly
# load canu
module add UHTS/Assembler/canu/1.3
# run assembly
canu -p canu_RSII -d $work_dir/lambda_RSII/assembly_canu/ genomeSize=48k -pacbio-raw $work_dir/lambda_RSII/RSII_reads_subset.fastq -useGrid=false -maxThreads=4

## post-assembly
# copy out the assembled contigs
cp $work_dir/lambda_RSII/assembly_canu/canu_RSII.contigs.fasta $work_dir/lambda_RSII/canu_RSII_contigs.fa
# cleaning modules 
module purge


### miniasm

## prepare directory
#create working directory for miniasm assembly
mkdir -p $work_dir/lambda_RSII/assembly_miniasm

## pre-assembly
# load minimap
module add UHTS/Analysis/minimap/0.2.r124.dirty
# compute overlaps
minimap -Sw5 -L100 -m0 -t4 $work_dir/lambda_RSII/RSII_reads_subset.fastq $work_dir/lambda_RSII/RSII_reads_subset.fastq | gzip -1 > $work_dir/lambda_RSII/asm_miniasm/RSII_overlaps.paf.gz

## assembly
# load miniasm 
module add UHTS/Analysis/miniasm/0.2.r137.dirty
# run assembly
miniasm -f $work_dir/lambda_RSII/RSII_reads_subset.fastq $work_dir/lambda_RSII/asm_miniasm/RSII_overlaps.paf.gz > $work_dir/lambda_RSII/assembly_miniasm/lambda_miniasm_contigs.gfa

## post-assembly
# conversion to fasta
awk '/^S/{print ">"$2"\n"$3}' $work_dir/lambda_RSII/assembly_miniasm/lambda_miniasm_contigs.gfa | fold > $work_dir/lambda_RSII/miniasm_RSII_contigs.fa

# last tested on: Mon 04 Jul 2016 11:33:07 PM CEST 
# AE
