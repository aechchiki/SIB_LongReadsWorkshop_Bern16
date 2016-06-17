# SIB workshop "Bioinformatics of long read sequencing" hands-on

July 5th, 2016. University of Bern. 

https://www.isb-sib.ch/training/upcoming-training-events/training/2016-07-longreads

-

## Contributors: 
- Amina Echchiki: Evolutionary Bioinformatics Group, UniL and SIB
- Walid Gharib: Interfaculty Bioinformatics Unit, UniBe and SIB
- Kamil Jaron: Evolutionary Bioinformatics Group, UniL and SIB
- Julien Roux: Evolutionary Bioinformatics Group, UniL and SIB

## Introduction
- biological source: DNA from phage lambda
  describe the protocol. Was there fragmentation?

- sequencing technologies: MinION, PacBio
- aims of this hands-on: reads extraction, quality control, mapping to a reference genome, genome assembly


## NanoOK

Short intro about this tool/wrapper

TO DO: generate full report with nanoOK. See `results/Emmanuel/sample.pdf` as an example.


NanoOK v0.72 (10_Nanook)
complete Analysis (fast5->fasta extraction->aln->bam->stats)
https://github.com/TGAC/NanoOK
A tool for extracting reads, aligning and producing a PDF containing a
comprehensive analysis report

Documentation: 
https://documentation.tgac.ac.uk/display/NANOOK/NanoOK

## TO DO add modules
This should work:
    Type java -h to see the help text for Java.
    Type R -h to see if R is installed.
    Type pdflatex -v to see if LaTeX is installed.
    Type h5dump -h to see if HDF5 tools are installed.
    Type lastal -h to see help text for LAST.


# load needed modules for the session 
```sh 
module add R/3.2.2;
module add SequenceAnalysis/SequenceAlignment/last/531 # load last (alignment)


### module add UHTS/Analysis/poretools/0.5.1 # load poretools (read extraction)
### module add UHTS/Aligner/bwa/0.7.13 # load bwa (alignment)
### module add UHTS/PacBio/blasr/20140829;
### module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
### module add UHTS/Analysis/seqtk/2015.10.15 # fastq to fasta
```

```sh
export NANOOK_DIR=/home/jroux/bin/NanoOK
export PATH=/home/jroux/bin/NanoOK/bin:$PATH

mkdir -p lambda_minion/fast5

tar -xvf /home/jroux/archive/MinION/run_MinION_2015_06_15_Burn-In.tar 
## Takes some time! 
## TO DO: extract in adavance and touch these files!
## TO DO: should we rename the files to remove the HiSeq header? It could be confusing

ln -s /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/data/Downloads_burnin/*fast5 lambda_minion/fast5/
## Be sure to check the read-only permissions

nanook extract -q -s lambda_minion
```
What folder were created by nanoOK

question: look at fastq file
What does it mean? How many lines per read? How quality encoded? Does it look good?

What are Template / Complement / 2D reads?


The reference genome first needs to be indexed with the LAST aligner:
```sh
mkdir -p lambda_minion/reference/
cd lambda_minion/reference/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
mv NC_001416.fna lambda_ref_genome.fa
lastdb -Q 0 lambda_ref_genome lambda_ref_genome.fa
```

Now launch the alignment:
```
cd ../../
nanook align -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
```
Whta directory has been created? Whta does it contain? Look at a maf file
## TO DO: try with other aligners? First need to index the reference
for BWA-MEM:
```sh
module add UHTS/Aligner/bwa/0.7.13
bwa index referencename.fasta
nanook align -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa -aligner bwa
```

Note: it is possible to specify the alignment parameters. default are:
nanook align -s SampleDir -r referencename.fasta -alignerparams "-s 2 -T 0 -Q 0 -a 1"

TO DO? refer to LAST documentation?

```sh
nanook analyse -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
```
## Long! Problem?
nanoOK created a PDF directory. Open it.
