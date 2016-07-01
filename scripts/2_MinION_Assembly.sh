#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J MinION_asm
#BSUB -e MinION-assembly.err
#BSUB -o MinION-assembly.out
#BSUB -q dee-hugemem
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -R "rusage[swp=20000]"
#BSUB -v 20000000
#BSUB -M 8000000

## prepare directory

# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
module add UHTS/Assembler/canu/1.3;

# canu assembly
canu -p canu_MinION -d $work_dir/lambda_MinION/asm_canu/ genomeSize=48k -nanopore-raw $work_dir/lambda_MinION/MinION_Lambda2D.fastq -useGrid=false -maxThreads=4 

# copy out the assembled contigs
cp $work_dir/lambda_MinION/asm_canu/canu_MinION.contigs.fasta $work_dir/lambda_MinION/canu_MinION_contigs.fa

#create working directory for miniasm assemly
mkdir -p $work_dir/lambda_MinION/asm_miniasm

# module purge will remove canu modules
module purge
module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;

# overlaps
minimap -Sw5 -L100 -m0 -t4 $work_dir/lambda_MinION/MinION_Lambda2D.fastq $work_dir/lambda_MinION/MinION_Lambda2D.fastq | gzip -1 > $work_dir/lambda_MinION/asm_miniasm/minion_overlaps.paf.gz

# assembly
miniasm -f $work_dir/lambda_MinION/MinION_Lambda2D.fastq $work_dir/lambda_MinION/asm_miniasm/minion_overlaps.paf.gz > $work_dir/lambda_MinION/asm_miniasm/lambda_miniasm_contigs.gfa

# conversion to fasta
awk '/^S/{print ">"$2"\n"$3}' $work_dir/lambda_MinION/asm_miniasm/lambda_miniasm_contigs.gfa | fold > $work_dir/lambda_MinION/miniasm_MinION_contigs.fa
