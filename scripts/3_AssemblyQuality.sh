#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J AssemblyQuality
#BSUB -q priority 

## prepare directory
# define working directory for user
work_dir=/scratch/beegfs/weekly/$USER
# create working directory for reference download
mkdir -p $work_dir/reference
# create working directory for quast report
mkdir -p $work_dir/assembly_quality/quast

## download reference genome
# cd to reference directory  
cd $work_dir/reference
# get reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
# rename reference
mv NC_001416.fna lambda_ref_genome.fa 

## quast
# load module 
module add UHTS/Quality_control/quast/4.1
# generate report
cd $work_dir/assembly_quality/quast
quast.py -R $work_dir/reference/lambda_ref_genome.fa $work_dir/lambda_minion/canu_minion_contigs.fa $work_dir/lambda_RSII/canu_RSII_contigs.fa $work_dir/lambda_RSII/miniasm_RSII_contigs.fa

# last tested on: Mon 04 Jul 2016 11:33:07 PM CEST 
# AE
