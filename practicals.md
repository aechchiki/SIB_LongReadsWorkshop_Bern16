# SIB workshop "Bioinformatics of long read sequencing" hands-on

July 5th, 2016. University of Bern. 

https://www.isb-sib.ch/training/upcoming-training-events/training/2016-07-longreads

-

## Contributors: 
- Amina Echchiki: Evolutionary Bioinformatics Group, UNIL and SIB
- Walid Gharib: Interfaculty Bioinformatics Unit, UNIBE and SIB
- Kamil Jaron: Evolutionary Bioinformatics Group, UNIL and SIB
- Julien Roux: Evolutionary Bioinformatics Group, UNIL and SIB

## Introduction
The aim of this practicals session is for you to get your hands on real long reads sequencing data generated from two different technologies:
* Oxford Nanopore MinION
* Pacific Bioscience RS II 

The biological material that was sequenced using these two platforms is some DNA from the lambda phage (http://en.wikipedia.org/wiki/Lambda_phage). This is not a particularly interesting genomic material for a long read sequencing study, since such a small genome can be assembled easily with short Illumina reads (see for example https://peerj.com/articles/2055/). However, it is small (48kb), which makes it feasible to run an analysis yourself during the limited time of this practicals session.

You will go through different steps, which include the extraction of reads from their native encoding formats (HDF5 formats, see http://en.wikipedia.org/wiki/Hierarchical_Data_Format), their quality control, their mapping to a reference genome, and a **de-novo** genome assembly. Most of these steps will be performed on MinION data only, but the assembly step will also be performed on PacBio data.

*TO DO Julien: describe quickly the sequencing protocol. Fragmentation with Covaris sonication to yield ~8kb fragments. Used 1ug for library preparation, but protocols exist for smaller amounts. For pacBio, 10kb fragmentation or selection? How long was the run? Desrcibe that the run produces ~500 raw fast5 files, which are uploaded to the cloud for basecalling. Then, new fast5 files (1 per read) downloaded back*

<!--
Details of protocol: https://community.nanoporetech.com/protocols/experiment-companion-for-control-dna/v/cde_1001_v1_revm_18may2016-374
-->

## How to connect to the vital-it cluster?
*TO DO Walid. Probably they should not work on the home directory? In that case, tell them to create a directory at their name in ```/scratch/beegfs/weekly/```? Mention that they can use their vital-IT account if they have one.*

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
...

### Miniasm
...

## Read extraction
<!--
DONE (Julien): 
cd /scratch/beegfs/monthly/jroux/tp_long_reads/MinION_lambda_reads
cp /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/Downloads_burnin/* .
rename 'HiSeq_' '' *.fast5
-->

![To do](wrench-and-hammer.png)
First, create a working directory and create links to the raw sequencing files.

```sh
mkdir -p lambda_minion/fast5/
ln -s /scratch/beegfs/monthly/jroux/tp_long_reads/MinION_lambda_reads/*.fast5 lambda_minion/fast5/
```
![Question](round-help-button.png)
Can you guess what each file corresponds to? How many reads has this sequencing experiment produced?

![Tip](elemental-tip.png)
You can get a glimpse of the structure and organization of a ```.fast5/``` file using the HDF5 command-line tools:
```sh
h5dump file.fast5 | less
h5stat file.fast5
```

![To do](wrench-and-hammer.png)
You will need to convert the MinION raw reads from their ```.fast5``` to a more classical and readable ```.fasta```. This can be done with the ```nanook extract``` utility, and takes around 10 minutes.
<!--
nanook extract -fasta -s lambda
-->

![Question](round-help-button.png)
What folders were created by ```nanoOK```? What do the ```2D```, ```Template``` and ```Complement``` folders represent?

### Bonus (or at home)
The ```.fasta``` format is nice and simple, but does not include any information on the quality of the sequences produced. The ```.fastq``` format however has this information. Launch ```nanook extract``` with the option to extract to ```.fastq``` format (you don't need to let ```nanook extract``` finish the extraction of all reads, a few dozens should be enough).  Choose one file at random in the ```fastq/2D``` folder and compare the quality scores with the corresponding file in the ```fastq/Template``` folder. You can refer to this page for help on the PHRED quality scores in ```.fastq``` files: http://en.wikipedia.org/wiki/FASTQ_format#Encoding.
<!--
nanook extract -fastq -s lambda
-->

![Question](round-help-button.png)
Do the quality scores seem to be improved in 2D reads? 

## Mapping to a reference genome
You will download the lambda phage reference genome and map the reads to it using the ```LAST``` aligner. To map to a reference genome, the reference sequence first needs to be indexed.

```sh
mkdir -p lambda_minion/reference/
cd lambda_minion/reference/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
mv NC_001416.fna lambda_ref_genome.fa # we rename the file to a clearer name
lastdb -Q 0 lambda_ref_genome lambda_ref_genome.fa
cd ../../
```
![Tip](elemental-tip.png)
The size of the indexed sequences can be visualized in the ```.sizes``` file.

![To do](wrench-and-hammer.png)
Now launch the alignment with ```LAST``` (default aligner) using the ```nanook align``` utility.

<!--
nanook align -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
-->

![Question](round-help-button.png)
What directories have been created? What do they contain? Look at one randomly chosen alignment file. What is striking?

![Tip](elemental-tip.png)
It is possible to specify to ```nanook align``` the alignment parameters to be used with ```LAST```, with the ```-alignerparams``` option. By default the parameters are ```-s 2 -T 0 -Q 0 -a 1```.

### Bonus (or at home)
You can try to align reads with other aligners. For example to use ```BWA-MEM```, you first need to load the corresponding module on vital-it (```module add UHTS/Aligner/bwa/0.7.13```), create the index of the reference sequence (```bwa index reference.fasta```). Then you can relaunch ```nanook align``` with the ```-aligner bwa``` option.

## Statistics and QC report
![To do](wrench-and-hammer.png)
Launch the generation of the final report including QC and alignment statistics using the ```nanook analyse utility```. If no PDF file is present in the ```latex_last_passfail/``` folder, you can generate it yourself with the ```pdflatex file.tex``` command (press enter everytime the program prompt you with a question).

<!--
nanook analyse -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
pdflatex lambda_minion/latex_last_passfail/lambda_minion.tex
-->

*TO DO: questions on the report:*
*What is most common error type (indels, homopolymer). Does it make sense? Is there a systematic error trend? (Some GC bias and repeated errors could be due to PCR?)*

## Assembly
* TO DO Amina and Kamil*

## TO DOs
* check this paper https://dx.doi.org/10.7554/eLife.14258
* check the porecamp material (analysis part): http://porecamp.github.io/timetable.html
* Will these be needed?
```sh
module add UHTS/Analysis/poretools/0.5.1 # load poretools (read extraction)
module add UHTS/Aligner/bwa/0.7.13 # load bwa (alignment)
module add UHTS/PacBio/blasr/20140829;
module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
module add UHTS/Analysis/seqtk/2015.10.15 # fastq to fasta
```

<!--
![Question](round-help-button.png)
![Tip](elemental-tip.png)
![To do](wrench-and-hammer.png)
![Warning](warning.png)
-->