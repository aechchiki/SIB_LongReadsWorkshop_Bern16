# SIB workshop "Bioinformatics of long read sequencing" hands-on

July 5th, 2016. University of Bern. 

https://www.isb-sib.ch/training/upcoming-training-events/training/2016-07-longreads

-

## Contributors: 
- Amina Echchiki: Evolutionary Bioinformatics Group, UNIL and SIB
- Walid Gharib: Interfaculty Bioinformatics Unit, UNIBe and SIB
- Kamil Jaron: Evolutionary Bioinformatics Group, UNIL and SIB
- Julien Roux: Evolutionary Bioinformatics Group, UNIL and SIB

## Introduction
The aim of this practicals session is for you to get your hands on real long reads sequencing data generated from two different technologies:
* Oxford Nanopore MinION
* Pacific Bioscience RS II 

The biological material that was sequenced using these two platforms is some DNA from the lambda phage (http://en.wikipedia.org/wiki/Lambda_phage). This is not a particularly interesting genomic material for a long read sequencing study, since such a small genome can be assembled easily with short Illumina reads (see for example https://peerj.com/articles/2055/). However, it is small (48kb), which makes it feasible to run an analysis yourself during the limited time of this practicals session.

You will go through different steps, which include the extraction of reads from their native encoding formats (HDF5 formats, see http://en.wikipedia.org/wiki/Hierarchical_Data_Format), their quality control, their mapping to a reference genome, and a **de-novo** genome assembly. Most of these steps will be performed on MinION data only, but the assembly step will also be performed on PacBio data.

*TO DO Julien: describe quickly the sequencing protocol. Check how big was the fragmentation step for MinION. For pacBio, 10kb fragmentation or selection?*

## How to connect to the vital-it cluster?
*TO DO Walid. Mention that they can use their vital-IT account if they have one*

## Tools
Although the long reads sequencing technologies are quite recent, there is already a variety of tools available for their analysis. In the practicals, we will make use of the following tools, which are convenient, fast and well-performing. But we strongly encourage you to try other tools as well for your own analyses!
### NanoOK
*TO DO: short intro about this tool. Very helpful because allows numerous analyses. fast5 -> fasta/fastq extraction -> alignment -> stats and QC analysis). Note: made primarily for MinION, but easy to hack to use PacBio reads.*

* Paper: http://bioinformatics.oxfordjournals.org/content/32/1/142.full
* GitHub page: http://github.com/TGAC/NanoOK
* Documentation: http://documentation.tgac.ac.uk/display/NANOOK/NanoOK

The ```NanoOK``` module is not installed on vital-IT, but you can access the software at this path: ```/home/jroux/bin/NanoOK/```. 

![To do](wrench-and-hammer.png)
To make it work, you need to set up the following environment variables:
```sh
export NANOOK_DIR=/home/jroux/bin/NanoOK
export PATH=/home/jroux/bin/NanoOK/bin:$PATH
```
You also need to load the modules for the softwares required by ```NanoOK```:
```sh 
# Load R
module add R/3.2.2;
# Load the LAST aligner
module add SequenceAnalysis/SequenceAlignment/last/531
```
To test if everything is fine for ```NanoOK```, these commands should work:
* Type ```java -h``` to see the help text for Java.
* Type ```R -h``` to see if R is installed.
* Type ```pdflatex -v``` to see if LaTeX is installed.
* Type ```h5dump -h``` to see if HDF5 tools are installed.
* Type ```lastal -h``` to see help text for LAST.

### Minimap
* Paper: http://bioinformatics.oxfordjournals.org/content/early/2016/05/01/bioinformatics.btw152.full

### Miniasm

## Read extraction
You will convert the MinION reads from their ```.fast5``` to a more readable ```.fastq``` format.

*TO DO Julien: copy the files from /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/Downloads_burnin/ to my home directory. Rename them to remove the HiSeq header which is confusing... Make the directory readble only*
 
![Question](round-help-button.png)
How many reads?

*TO DO: take a small file (HiSeq_lambda_5830_1_ch12_file4_strand.fast5) and use hdf5 unix commands to visualize/summarize it*


![To do](wrench-and-hammer.png)
Create a working directory
```sh
mkdir -p lambda_minion/fast5/
ln -s /home/jroux/.../*.fast5 lambda_minion/fast5/
```

![To do](wrench-and-hammer.png)
Extract the reads:
```sh
nanook extract -q -s lambda_minion
```

What folder were created by nanoOK?

question: look at fastq file
What does it mean? How many lines per read? How quality encoded? Does it look good?

What are Template / Complement / 2D reads folders?


## Mapping to a reference genome
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
What directory has been created? Whta does it contain? Look at a maf file
*TO DO: try with other aligners? First need to index the reference*
for BWA-MEM:
```sh
module add UHTS/Aligner/bwa/0.7.13
bwa index referencename.fasta
nanook align -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa -aligner bwa
```

Note: it is possible to specify the alignment parameters. default are:
nanook align -s SampleDir -r referencename.fasta -alignerparams "-s 2 -T 0 -Q 0 -a 1"

TO DO? refer to LAST documentation?

## Statistics and QC report
```sh
nanook analyse -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
```
*Long! Problem? Missing module? Problem because only 1 read was analyzed?*
nanoOK created a PDF directory. Open it.

TO DO: questions to ask:
What is most common error type (indels, homopolymer). Does it make sense?
Is there a systematic error trend?
(Some GC bias and repeated errors could be due to PCR)



## Assembly
TO DO: refer to minimap/miniasm paper



TO DO: check this paper https://dx.doi.org/10.7554/eLife.14258
TO DO: check the porecamp material (analysis part): http://porecamp.github.io/timetable.html

TO DO: needed?
```sh
module add UHTS/Analysis/poretools/0.5.1 # load poretools (read extraction)
module add UHTS/Aligner/bwa/0.7.13 # load bwa (alignment)
module add UHTS/PacBio/blasr/20140829;
module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
module add UHTS/Analysis/seqtk/2015.10.15 # fastq to fasta
```

![Question](round-help-button.png)
![Tip](elemental-tip.png)
![To do](wrench-and-hammer.png)
![Warning](warning.png)

<!--
-->