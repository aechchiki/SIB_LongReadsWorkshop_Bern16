#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J Mapping

## prepare directory
# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
# create working directory for mapping
mkdir -p $work_dir/mapping_minion
# create working directory for minion fasta
mkdir -p $work_dir/lambda_minion/fasta/
mkdir -p $work_dir/lambda_minion/fasta/2D
mkdir -p $work_dir/lambda_minion/fasta/Complement
mkdir -p $work_dir/lambda_minion/fasta/Template

## pre-assembly
# load LAST aligner
module add SequenceAnalysis/SequenceAlignment/last/531
# load NanoOk
module add UHTS/Analysis/NanoOK/0.72
# index the reference 
lastdb -Q 0 $work_dir/reference/lambda_ref_genome $work_dir/reference/lambda_ref_genome.fa
# extract the read headers and sequences only
cat $work_dir/lambda_minion/all_reads.fastq | awk 'NR%4==1||NR%4==2' > $work_dir/lambda_minion/all_reads.fasta 
# convert header according to fasta specifications
perl -p -i -e 's/^@/>/g' $work_dir/lambda_minion/all_reads.fasta 
# split each read in a different file
nanook_split_fasta -i $work_dir/lambda_minion/all_reads.fasta -o $work_dir/lambda_minion/fasta/
# move reads to folders according to their type: 2D, template, complement
find $work_dir/lambda_minion/fasta/ -maxdepth 1 -name \*_2d.fasta -exec mv {} $work_dir/lambda_minion/fasta/2D/  \;
find $work_dir/lambda_minion/fasta/ -maxdepth 1 -name \*_template.fasta -exec mv {} $work_dir/lambda_minion/fasta/Template/  \;
find $work_dir/lambda_minion/fasta/ -maxdepth 1 -name \*_complement.fasta -exec mv {} $work_dir/lambda_minion/fasta/Complement/ \;

# run alignment
nanook align -s $work_dir/lambda_minion -r $work_dir/reference/lambda_ref_genome.fa

# get report
nanook analyse -s $work_dir/lambda_minion -r $work_dir/reference/lambda_ref_genome.fa
pdflatex lambda_minion/latex_last_passfail/lambda_minion.tex


# last tested on: Mon 04 Jul 2016 11:33:07 PM CEST 
# AE
