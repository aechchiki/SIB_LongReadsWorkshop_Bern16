#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J minion_assembly
#BSUB -e minion-assembly.err
#BSUB -o minion-assembly.out
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
mkdir -p $work_dir

## run assembly
# load canu 
module add UHTS/Assembler/canu/1.3
# prepare assembly directory
mkdir -p $work_dir/lambda_minion/assembly_canu/
# run assembly
canu -p canu_minion -d $work_dir/lambda_minion/assembly_canu/ genomeSize=48k -nanopore-raw $work_dir/lambda_minion/all_reads.fastq -useGrid=false -maxThreads=4 

## post-assembly processing
# copy out the assembled contigs
cp $work_dir/lambda_minion/assembly_canu/canu_minion.contigs.fasta $work_dir/lambda_minion/canu_minion_contigs.fa

## cleaning
# purge loaded modules
module purge


### miniasm

## prepare directory
#create working directory for miniasm assemly
mkdir -p $work_dir/lambda_minion/assembly_miniasm

## run assembly
# load minimap and miniasm
module add UHTS/Analysis/minimap/0.2.r124.dirty
module add UHTS/Analysis/miniasm/0.2.r137.dirty
# compute overlaps
minimap -Sw5 -L100 -m0 -t4 $work_dir/lambda_minion/all_reads.fastq $work_dir/lambda_minion/all_reads.fastq | gzip -1 > $work_dir/lambda_minion/assembly_miniasm/minion_overlaps.paf.gz
# run assembly
miniasm -f $work_dir/lambda_minion/all_reads.fastq $work_dir/lambda_minion/assembly_miniasm/minion_overlaps.paf.gz > $work_dir/lambda_minion/assembly_miniasm/lambda_miniasm_contigs.gfa

## post-assembly processing
# conversion to fasta
awk '/^S/{print ">"$2"\n"$3}' $work_dir/lambda_minion/assembly_miniasm/lambda_miniasm_contigs.gfa | fold > $work_dir/lambda_minion/miniasm_minion_contigs.fa

# last tested on: Mon 04 Jul 2016 11:33:07 PM CEST 
# AE
