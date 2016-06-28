---
output: html_document
---
# SIB workshop "Bioinformatics of long read sequencing" hands-on

July 5th, 2016. University of Bern. 

https://www.isb-sib.ch/training/upcoming-training-events/training/2016-07-longreads


## Contributors

- Amina Echchiki: Evolutionary Bioinformatics Group, UNIL and SIB
- Walid Gharib: Interfaculty Bioinformatics Unit, UNIBE and SIB
- Kamil Jaron: Evolutionary Bioinformatics Group, UNIL and SIB
- Julien Roux: Evolutionary Bioinformatics Group, UNIL and SIB

## Introduction

The aim of this practicals session is for you to get your hands on real long reads sequencing data generated from two different technologies:

- Oxford Nanopore [MinION](https://www.nanoporetech.com/products-services/minion-mki)
- Pacific Biosciences [RSII](http://www.pacb.com/products-and-services/pacbio-systems/rsii/)

The biological material sequenced using these two platforms is DNA from the [lambda phage](http://en.wikipedia.org/wiki/Lambda_phage). This is not a particularly interesting genomic material for a long read sequencing study, since such a small genome can be assembled easily with short Illumina reads (check [this](https://peerj.com/articles/2055/) publication for example ). However, its genome is small (48kb), which makes it feasible to run an analysis yourself during the limited time of this practicals session.

You will go through different steps, which include the extraction of reads from their native encoding formats ([HDF5](http://en.wikipedia.org/wiki/Hierarchical_Data_Format)), their quality control, their mapping to a reference genome, and a *de-novo* genome assembly.

The MinION library preparation protocol starts by DNA fragmentation, then the fragmented DNA is end-repaired and dA-tailed. Adapters are ligated to the dsDNA fragments. These adapters are in two flavors: Y-form and hairpin-form, allowing the generation of 2D reads. Both adapters are ligated to each end of the dsDNA fragments. The adapters are conjugated with motor proteins that help control the translocation speed of DNA through the pore. An illustration of the protocol is provided in this [figure](http://www.genetics.org/content/202/1/37):

![protocol](img/protocol.jpg)

In the MinION experiment you are going to analyze, the lambda phage DNA was fragmented using sonication (Covaris) to yield ~8kb-long fragments. One microgram of fragmented DNA material was then used for library preparation. Other protocols exist for smaller amounts of starting material. 

The sequencing run produced a single ```.fast5``` file per pore. All files were then uploaded to the cloud for basecalling. The base-called files were downloaded back to vital-it. They are also in the ```.fast5``` format, but there is one file per read.

<!--
Details of protocol: https://community.nanoporetech.com/protocols/experiment-companion-for-control-dna/v/cde_1001_v1_revm_18may2016-374
-->

***

## How to connect to the Vital-IT cluster?

By now, you should have received a username and password to access the high performance computing cluster of the SIB (Vital-IT).

Note 1: Even if you have access to another cluster, you will not be able to access the data which is stored on Vital-IT so all participants should connect to the latter.

Note 2: Obviously, you can use your own Vital-IT account if you have one.

In order to connect to the cluster and set up your working directory follow the next steps:

### For Linux/OSX users

* In a terminal, type:

```ssh <username>@prd.vital-it.ch```

* You will be prompted to input your user password: 

```<username>@prd.vital-it.ch's password:```

type in your password (you will not see what you are typing) and press Enter.

* You are in! Jump to: "Setting up your working directory".

### For Windows users

You should first have a ssh client installed e.g. PuTTY.

* If you already have it, we assume you know how to use it: connect to Vital-it and jump to point 3.

* If you don't have an ssh client, follow the following link steps. When done, proceed to point 3. https://github.com/aechchiki/SIB_LongReadsWorkshop_Bern16/blob/master/vital-it_connect_Putty.pdf

### Setting up your working directory

On Vital-IT, it is highly recommended yet mandatory to read and write in the ```/scratch``` directory:

* move to `/scratch`: 
    
```sh
cd /scratch/beegfs/weekly/
```

* create your own directory:
    
```sh
mkdir <username> ; cd <username>
```

* You will always be working from this directory. Before launching commands please be sure that you are located in the right directory by typing:

```sh
pwd
```
Expected output: ```/scratch/beegfs/weekly/<username>```

    
***

## 1. Read extraction

The output of MinION and PacBio RSII is stored in Hierarchical Data Format, what is basically an archive format (like `.zip` or `.tar`), but with very quick access to its content (You can find details on [wikipedia](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)). In those files are more details about reads and basecalling. What, we need is just a simple `.fastq` for downstream analysis.

**fast5**

MinION basecaller produces one file per read. In the file there is a time series of a current used for basecalling, basecalled template and complement reads (also called 1D reads) and consensus of 1D reads called 2D reads with comparably higher accuracy.

**bas.h5 & bax.h5**

RSII system produces three `bax.h5` files and a `bas.h5` per SMRT cell. The three `bax.h5` files correspond to the first, second and third part of the movie capturing the SMRT cell, `bas.h5` contains metainformations. PacBio announced a change of data format to specialised `.bam` for platform Sequel, however all data produced by RSII will be still in `.h5` formats we are going to work with.


### Tools

**poretools** is a toolkit to analyze MinION data, including the extraction of read sequences from native ```.fast5``` files. Source code is available on [GitHub](https://github.com/arq5x/poretools) and software usage is detailed in [documentation](http://poretools.readthedocs.io/en/latest/). **poretools** is available on Vital-IT: 

```sh
module add UHTS/Analysis/poretools/0.5.1
```

**pbh5tools** is a set of python scripts to extract `.fasta` and `.fastq` from `bas.h5` and `bax.h5` files. These scripts allow filtering based on error rate, read length and read type. Source code is available on [GitHub](https://github.com/PacificBiosciences/pbh5tools) and software usage is detailed in [documentation](https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst). **pbh5tools** is available on Vital-IT: 

```sh
module add UHTS/PacBio/pbh5tools/0.8.0;
```

Read types of SMRT cell are: 

<!--
syntax of scripts to tools?

text or picture??

- polymerase reads: The full read
- subreads: Sebreads separated by cutting out adapters.
- circular consensus read: The consensus seqeunce of subreads. Note, that `ccs` seqeunce does not have to be computed in the `h5` files.
-->

![protocol](img/pb_reads.png)

### MinION

![To do](img/wrench-and-hammer.png)
First, create a working directory for MinION raw reads and create links to the raw sequencing files.

```sh
mkdir -p lambda_minion/fast5/
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/minion_lambda_reads/*.fast5 lambda_minion/fast5/
```
![Question](img/round-help-button.png)
Can you guess what each file corresponds to? How many reads has this sequencing experiment produced? [TO ADD]

```sh
ls lambda_minion/fast5/ | wc -l
```

![Tip](img/elemental-tip.png)
You can get a glimpse of the structure and organization of a ```.fast5/``` file using the HDF5 command-line tools:

```sh
h5dump file.fast5 | less
h5stat file.fast5
```

![To do](img/wrench-and-hammer.png)
You will need to convert the MinION raw reads from their ```.fast5``` to a more classical and readable ```.fastq```. This can be done with the ```poretools fastq``` utility, and takes some seconds.

```sh
bsub -q priority "poretools fastq --type 2D lambda_minion/fast5/*.fast5 >  lambda_minion/MinION_Lambda2D.fastq"
```

### Bonus (or at home)
The ```.fasta``` format is nice and simple, but does not include any information on the quality of the sequences produced. The ```.fastq``` format however has this information. Launch ```nanook extract``` with the option to extract to ```.fastq``` format (you don't need to let ```nanook extract``` finish the extraction of all reads, a few dozens should be enough).  Choose one file at random in the ```fastq/2D``` folder and compare the quality scores with the corresponding file in the ```fastq/Template``` folder. You can refer to this page for help on the PHRED quality scores in ```.fastq``` files: http://en.wikipedia.org/wiki/FASTQ_format#Encoding.

```sh
nanook extract -fastq -s lambda
```

![Question](img/round-help-button.png)
Do the quality scores seem to be improved in 2D reads? 

![Tip](img/elemental-tip.png)
To access quality values, you can use `fastqc` known from Illumina data. Reported stats are correct, just keep in mind, that warning flags in the report are for assembly by short reads and therefore not very informative.

![help](img/help.png) If you are lost, you can get extracted MinION reads by executing
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/1_Extraction_MinION.sh
```

### RSII
```sh
module add UHTS/PacBio/pbh5tools/0.8.0;
```

![To do](img/wrench-and-hammer.png) Extract also RSII reads using 

```sh
mkdir -p lambda_RSII/raw_reads/
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/RSII_lambda_reads/*.h5 lambda_RSII/raw_reads/
```

![Question](img/round-help-button.png)
How reads of many SMRT cells we have? [1]

![To do](img/wrench-and-hammer.png) Extract RSII reads using `bash5tools.py` from `pbh5tools`.

```sh
bsub -q priority bash5tools.py <file.bas.h5> --outFilePrefix <prefix_for_extracted_reads> --readType subreads --outType fastq
```

![Question](img/round-help-button.png)
How many reads has was produced by the SMRT cell? [132269]


```sh
wc -l <prefix_for_extracted_reads>.fastq   # reads in fastq = number of lines / 4
```

![help](img/help.png) If you are lost, you can get extracted RSII reads by executing

```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/1_Extraction_RSII.sh
```

***

## 2. Genome assembly

Long reads genome assembly is just like established short-reads assembly based on graph theory, but string graphs are used instead of De Brujin. Instead of kmer decomposition, the graph of reads (nodes) and their overlaps (edges) is constructed.

There are two leading assemblers for long reads Canu and Falcon. Canu is a fork of Celera with added modules for error correction and trimming of reads. Falcon is an assembler aware of ploidy developed from scratch for PacBio data.

### Tools

We will use Canu and Miniasm for both MinION and RSII assembly, then we check the assembly quality using the report provided by Quast. 

**Canu** is an assembler for noisy long-reads sequences. Canu will correct the reads, trim suspicious regions, then assemble the corrected and cleaned reads into unitigs. Note: this software was designed for both MinION and RSII data.

* Link to the paper: http://www.nature.com/nbt/journal/v33/n6/abs/nbt.3238.html
* GitHub page: https://github.com/marbl/canu
* Canu manual: http://canu.readthedocs.io/en/stable/

**Miniasm** is lightweight assembler. It very fast, but for a cost of simplicity and low parametrization. It can used for a first proxy of the data content, but for the final assembly, another assembler should be considered. The first step, the overlap of reads is computed in a separated step using a standalone program called `minimap`.

### MinION

We assemble lambda phage using the 3,068 2D reads (cca 500x). 1D reads are in substantially worse, quality, so if we have enough 2D reads, we can use just them. The assembly should be trivial: The genome is 48.8kbp and we have some reads of 20kb! We expect only one contig!

**Canu**

```sh
module add UHTS/Assembler/canu/1.3;
``` 

![To do](img/wrench-and-hammer.png) Chose appropriate parameters:

```
parameters:	
				-p : use specified prefix to name canu output files
				-d : use specified prefix for output directories 
				-errorRate : (optional) specifies expected error of reads.
				-genomeSize : expected genome size
				-nanopore-raw : specifies ONT MinION data 
				-useGrid=false : disables automatic submission to cluster
```

![To do](img/wrench-and-hammer.png) and run the assembly in the directory with extracted MinION 2D reads.

```sh
bsub -n 4 -q priority 'canu -p lambda -d <name_of_folder_for_output> errorRate=<error_rate> genomeSize=<genome_size> -nanopore-raw <reads_name> -useGrid=false -maxThreads=4 &> canu_minion_lambda.log'
```

The output is a directory containing several files. The most interesting ones are:
- *.correctedReads.fasta.gz : file containing the input sequences after correction, trim and split based on consensus evidence. 
- *.trimmedReads.fastq : file containing the sequences after correction and final trimming
- *.layout : file containing informations about read inclusion in the final assembly
- *.gfa : file containing the assembly graph by Canu
- *.contigs.fasta: file containing everything that could be assembled and is part of the primary assembly

![Question](img/round-help-button.png) How many contigs were produced? [1]

**Miniasm**

```sh
module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;
```

![To do](img/wrench-and-hammer.png) Compute overlaps of reads using `minimap`

```sh
bsub -q priority -n 4 'minimap -Sw5 -L100 -m0 -t4 <reads.fq> <reads.fq> | gzip -1 > <overlaps.paf.gz>'
```

![To do](img/wrench-and-hammer.png) Perform an assembly using `miniasm`.

```sh
bsub -q priority 'miniasm -f <reads.fq> <overlaps.paf.gz> > Lambda_contigs.gfa'
```

![To do](img/wrench-and-hammer.png) Convert `.gfa` to `.fasta`.

```sh
awk '/^S/{print ">"$2"\n"$3}' Lambda_contigs.gfa | fold > Lambda_contigs.fa
```

![Question](img/round-help-button.png) How is it different, the assembly computed by `Canu` and `Miniasm`? [not much]

### RSII

One SMRT cell produces yield between 0.5-1gbp. What roughly correspond to 10,000 to 20,000x coverage of lambda phage genome. As you can imagine, it is a bit overkill. We can a decrease a computational load very much by assembling a subset of a few thousand reads only.

![To do](img/wrench-and-hammer.png) Go to directory with extracted RSII reads and take a subset of the first 3,000 reads.

```sh
head -12000 RSII_reads.fastq > RSII_reads_subset.fastq
```

The assembly is now analogical to the assembly of MinION data

![To do](img/wrench-and-hammer.png) Modify in the command of `Canu` the type of input reads to `-pacbio-raw`, name of the input file and the name of output name and perform the assembly.


![To do](img/wrench-and-hammer.png) Perform assembly of RSII data using miniasm as well.

![Question](img/round-help-button.png) What is are the lengths of all four assemblies? Is there anything weird going on? [yes, RSII reads are not assembled well.]


![help](img/help.png) If you are lost, you can get both assemblies of MinION and RSII reads by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/2_assembly.sh
```

***

## 3. Quality of assembly

### Tools

**Quast** is a quality assessment tool for genome assemblies. Running Quast on an assembly provides a detailed report about the statistics of the assembly, in pdf format.

* Link to the paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/
* GitHub page: https://github.com/ablab/quast
* Quast manual: http://quast.bioinf.spbau.ru/manual

The `Quast` module is installed on vital-IT. You can access the software by loading the module: `module add UHTS/Quality_control/quast/4.1`. 

### Usage?

We will use the information in *.contigs.fasta to generate the assembly report.

```sh
# load the quality assessment module 
module add UHTS/Quality_control/quast/4.1
# generate the report 
quast.py -R <path/to/reference_genome> *.contigs.fasta
```

The output is a directory, typically `quast_results/<date_of_launch>` containing several files. The most interesting for us is the report in pdf format. We can copy it to the local machine using `scp`. 

TODO: interesting questions about the report?


TODO: questions on comparison of the assemblies 


***TO DO: Emmanuel has done a nanook analysis on these reads. transfer the PDF to your laptop using `scp`. Couldn't do it because permissions not set up correctly (email sent to Emmanuel).***

***TO DO: Emmanuel told me the DNA fragments for PacBio were 10kb long. Is that made through sonication or size-selection?***

![Question](img/round-help-button.png)
Have a look at the `nanoOK` report for these reads. How does it compare to the report made on MinION reads?


<!--
TO DOs
* check this paper https://dx.doi.org/10.7554/eLife.14258
* check the porecamp material (analysis part): http://porecamp.github.io/timetable.html
* Check Pradervand talk: http://edu.isb-sib.ch/pluginfile.php/3390/course/section/1574/Pradervand_Sequencing_CUSO2015.pdf
* Will these modules be needed?
```sh
module add UHTS/Analysis/poretools/0.5.1 # load poretools (read extraction)
module add UHTS/PacBio/blasr/20140829;
module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
module add UHTS/Analysis/seqtk/2015.10.15 # fastq to fasta
```
* Idea: make them download maf files and reference genome, and open in IGV... Or is it possible to generate a picture on vital-it from IGV?

![Question](img/round-help-button.png)
![Tip](img/elemental-tip.png)
![To do](img/wrench-and-hammer.png)
![Warning](img/warning.png)
-->

![help](img/help.png) If you are lost, you can get all quality stats of all assemblies by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/3_assembly_quality.sh
```

***

## 4. Mapping to a assembled genome

### Tools

You will download the lambda phage reference genome and map the reads to it using the `LAST` aligner. To map to a reference genome, the reference sequence first needs to be indexed.

### MinION

```sh
mkdir -p lambda_minion/reference/
cd lambda_minion/reference/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
mv NC_001416.fna lambda_ref_genome.fa # we rename the file to a clearer name
lastdb -Q 0 lambda_ref_genome lambda_ref_genome.fa
cd ../../
```
![Tip](img/elemental-tip.png)
The size of the indexed sequences can be visualized in the `.sizes` file.

![To do](img/wrench-and-hammer.png)
Now launch the alignment with `LAST` (default aligner) using the `nanook align` utility.

```sh
nanook align -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
```

![Question](img/round-help-button.png)
What directories have been created? What do they contain? Look at one randomly chosen alignment file. What is striking?

![Tip](img/elemental-tip.png)
It is possible to specify to `nanook align` the alignment parameters to be used with `LAST`, with the `-alignerparams` option. By default the parameters are `-s 2 -T 0 -Q 0 -a 1`.

### Bonus (or at home)
You can try to align reads with other aligners. For example to use `BWA-MEM`, you first need to load the corresponding module on vital-it (`module add UHTS/Aligner/bwa/0.7.13`), create the index of the reference sequence (`bwa index reference.fasta`). Then you can relaunch `nanook align` with the `-aligner bwa` option.

![help](img/help.png) If you are lost, you can perform mapping of RSII reads by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/4_RSII_mapping.sh
```

### RSII

![help](img/help.png) If you are lost, you can perform mapping of RSII reads by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/4_RSII_mapping.sh
```

### Statistics and QC report
![To do](img/wrench-and-hammer.png)
Launch the generation of the final report including QC and alignment statistics using the ```nanook analyse``` utility. A PDF file should be created in the ```latex_last_passfail/``` folder. If it's not the case, you can generate it yourself with the ```pdflatex file.tex``` command (press enter every time the program prompt you with a question).

<!--
nanook analyse -s lambda_minion -r lambda_minion/reference/lambda_ref_genome.fa
pdflatex lambda_minion/latex_last_passfail/lambda_minion.tex
-->

To download the report to your laptop:
```sh
scp username@prd.vital-it.ch:[path_to_file_on_vital_it] [path_to_file_on_laptop]
```

![Question](img/round-help-button.png)
Take some time to read and understand the report. Here are a few questions that will guide you:
* What size is the longest template read? Is that surprising?
* What does the N50 values indicate us? Is that consistent with the size of the DNA fragments used to create the library?
* What is most common error type within reads?
* In comparison, the typical error rate reported for the Illumina sequencing technology is around 0.1%. For Sanger sequencing, it can be as low as 0.001%... An interesting comparison of error rates across platforms was recently published: http://bib.oxfordjournals.org/content/17/1/154
* Why is the alignment rate of 2D reads higher than those of template and complement reads?
* Which are the most accurate: shortest or longest reads?
* Have a look at the coverage plot. There is a "bump" in coverage around 45kb. This corresponds to some control spike-in DNA that was added during library preparation (more precisely around 3.6kb of a region of the lambda phage genome, with a single mutation G45352A). Can additional variation be explained by GC content differences?

*TO DO: tell them to have a look at alignments in this region to see the mutation?*

* Have a look at the k-mer over and under-representation analysis. What sort of k-mers are under-represented in 2D reads? Is that expected given how the technology works?

### RSII

![help](img/help.png) If you are lost, you can perform a quality control report of mapping by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/4_mapping_qc_report.sh
```
