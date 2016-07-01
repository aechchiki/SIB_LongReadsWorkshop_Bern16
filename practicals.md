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

The aim of this practicals session is to get your hands on real long reads sequencing data generated from two different technologies:

- Oxford Nanopore [MinION](https://www.nanoporetech.com/products-services/minion-mki)
- Pacific Biosciences [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/)

The biological material that was sequenced using these two platforms is DNA from the [lambda phage](http://en.wikipedia.org/wiki/Lambda_phage). This is not a particularly interesting genomic material for a long read sequencing study, since such a small genome can be assembled easily with short Illumina reads (check [this](https://peerj.com/articles/2055/) publication for example). However, its genome is small (48kb), which makes it feasible to run an analysis yourself during the limited time of this practicals session.

You will go through different steps, which include the extraction of reads from their native encoding formats ([HDF5](http://en.wikipedia.org/wiki/Hierarchical_Data_Format)), quality control, *de novo* genome assembly, and a mapping of reads to a reference genome.

After DNA fragmentation, the **MinION** library is prepared by ligating adapters to the double-stranded DNA fragments. One side receives a Y-form adapter, and the other side a hairpin-form adaptor. The adapters are conjugated with motor proteins that help control the translocation speed of DNA through the pore. When the Y-form adapter approaches the pore, one strand starts to be sequenced. This sequence is called *template*. After the hairpin-form adapter is sequenced, the *complement* sequence is read. The template and complement sequences are sometimes called *1D reads*. The consensus of template and complement sequences is called a *2D read*. An illustration of the protocol is provided in this [figure](http://www.genetics.org/content/202/1/37):

![protocol](img/protocol.jpg)
<!--
TO DO: find a better picture of the read going through the pore?
-->

For **PacBio** sequencing (or **SMRT** sequencing), hairpin-form adapters (SMRTbell) are ligated on both sides of double-stranded DNA fragments, forming a circular sequence. Circular sequences attach to the bottom of theSMRT cell wells, and are continuously read, and a movie is recorded. Nowadays, the total length of movie allows to record up to 75,000 nucleotide basecalls. This implies that short DNA fragments will be read many times, while long fragments may be read just one or twice. The complete sequence is called a *polymerase* read. After the adapter sequences are removed, the sequence is split into *subreads*. The consensus of subreads is called a *circular consensus read* or *read of insert*.

![protocol](img/pb_reads.png)

<!--
In the MinION experiment you are going to analyze, the lambda phage DNA was fragmented using sonication (Covaris) to yield ~8kb-long fragments. One microgram of fragmented DNA material was then used for library preparation. Other protocols exist for smaller amounts of starting material. 

Details of protocol: https://community.nanoporetech.com/protocols/experiment-companion-for-control-dna/v/cde_1001_v1_revm_18may2016-374

PacBio: Note, that `ccs` sequence does not have to be computed in the `h5` files.
-->

## How to connect to the Vital-IT cluster?

By now, you should have received a username and password to access the high performance computing cluster of the SIB (Vital-IT). Obviously, you can use your own Vital-IT account if you have one.

![Warning](img/warning.png)
Even if you have access to another cluster, the data required for this practical is stored on Vital-IT so all participants should connect to vital-IT.

In order to connect to the cluster and set up your working directory follow the next steps:

### For Linux / Mac OSX users

In a terminal, type ```ssh <username>@prd.vital-it.ch```. You will be prompted to input your user password (```<username>@prd.vital-it.ch's password:```). Type in your password (you will not see what you are typing) and press Enter.

You are in! Jump to [Setting up your working directory](#setting-up-your-working-directory).

### For Windows users

You should first install a ssh client (e.g., `PuTTY`). If you already have one, we assume you know how to use it. Connect to Vital-IT and jump to [Setting up your working directory](#setting-up-your-working-directory).

If you do not have a ssh client, follow [these steps](https://github.com/aechchiki/SIB_LongReadsWorkshop_Bern16/blob/master/vital-it_connect_Putty.pdf). 

### Setting up your working directory

* On Vital-IT, you must read and write files in the "scratch" directory. Move to `/scratch` weekly directory:
```sh
cd /scratch/beegfs/weekly/
```

![Warning](img/warning.png)
Be careful, files in this folder are erased after a week. At the end of the practicals, if you want to keep your results, you need to back-up the data, for example in a compressed tarball that you move to your home or archive folder, or to another computer.

* Create your own directory:
```sh
mkdir <username>; cd <username>
```

* You will always be working from this directory. Before launching commands please be sure that you are located in the right directory by typing: ```pwd```. The expected output is: ```/scratch/beegfs/weekly/<username>```

### Submitting commands to the cluster

It is not allowed to launch any calculation on the frontal machine where you are connected (```prd.vital-it.ch```). You need to submit each job for batched execution through a job scheduler that will dispatch it on the cluster nodes. For example:
```sh
bsub '[my command line here]'
```
Or, better, write your commands in a script and submit it with:
```sh
bsub < script.sh
```
Please have a look at [this short tutorial](https://github.com/aechchiki/SIB_LongReadsWorkshop_Bern16/blob/master/vital-it-usage.md) to help you write such a script yourself, with the `Nano` text editor. To save you some typing, a skeleton of a submission script can be found [here](script/skeleton.sh) ;)  

<!--
TO DO?
* Memory usage? `-M / -R`
* Make them download a skeleton of a submission script, to save time in writting it?
-->

## 1. Read extraction

### Formats
The output of MinION and PacBio RS II are stored in Hierarchical Data Format, that is basically an archive format (like `.zip` or `.tar`), but allowing a very quick access to its content (you can find details on [wikipedia](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)). Those files include all details about reads and basecalling. What we need will need for this practical are just the sequences and their qualities, that can be stored in a simple `.fastq`.

![Tip](img/elemental-tip.png)
The files saved in Hierarchical Data Format can be explored using HDF5 command-line tools:
```sh
h5dump <HDF_file> | less
h5stat <HDF_file> | less
```

#### fast5
The MinION basecaller produces one file per read, which includes the time series of the measured ionic current used for basecalling, the basecalled template and complement read sequences, and the consensus 2D reads.

#### bas.h5 and bax.h5
The PacBio RS II system produces three `bax.h5` files and a `bas.h5` per SMRT cell. The three `bax.h5` files correspond to the first, second and third part of the movie capturing the SMRT cell. The `bas.h5` file contains metadata. PacBio announced a change of data format to more classcial `.bam` files for the next platform (called *Sequel*), however all data produced by the RS II platform will be still in the `.h5` formats we are going to work with.

### Tools
All used tools are already installed on Vital-IT. Before the first use, you will just need to load the package.

**Poretools** is a toolkit for working with sequencing data from MinION, for the purposes of quality control and downstream analysis. It operates directly on the native `.fast5` format and provides format conversion utilities and data exploration tools. Software usage is detailed in the [documentation](http://poretools.readthedocs.io/en/latest/).

**pbh5tools** is a set of python scripts to extract `.fasta` and `.fastq` from `bas.h5` and `bax.h5` files. These scripts allow filtering based on error rate, read length and read type. Source code is available on [GitHub](https://github.com/PacificBiosciences/pbh5tools) and software usage is detailed in [documentation](https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst). 

### Extraction of MinION reads
First, create a working directory for MinION raw reads and create links to the raw sequencing files.

```sh
mkdir -p lambda_minion/fast5/
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/MinION_lambda_reads/*.fast5 lambda_minion/fast5/
```

![Question](img/round-help-button.png)
Can you guess what each file corresponds to? How many reads has this sequencing experiment produced?

<!--
```sh
ls lambda_minion/fast5/ | wc -l
```
Answer: 5559 reads
-->

<!--
![To do](img/wrench-and-hammer.png)
You will convert the MinION raw reads from their `.fast5` to a more classical and readable `.fastq`. This can be done with the `nanook extract` utility, and takes around 10 minutes. Please refer to the documentation to find out how to run the program and which options to use.

```sh
module add UHTS/Analysis/NanoOK/0.72
bsub -q priority -o minion_extract.out -e minion_extract.err -R "rusage[mem=4096]" 'nanook extract [...]'
bsub -q priority -o minion_extract.out -e minion_extract.err -J minion_extract -R "rusage[mem=4096]" 'nanook extract -fastq -s lambda_minion'
```

![Question](round-help-button.png)
What folders are created by ```nanoOK```? What do they include?

![To do](img/wrench-and-hammer.png)
Choose one file at random in the ```fastq/2D``` folder and compare the quality scores with the corresponding file in the ```fastq/Template``` folder. You can refer to this page for help on the PHRED quality scores in ```.fastq``` files: http://en.wikipedia.org/wiki/FASTQ_format#Encoding.
-->

![To do](img/wrench-and-hammer.png)
You will convert the MinION raw reads from their `.fast5` to a more classical and readable `.fastq`. This can be done with the `poretools fastq` utility. To find out which options to use, please refer to the [documentation](http://poretools.readthedocs.io/en/latest/content/examples.html#poretools-fastq), or type `poretools fastq -h`. For now, we are interested in extracting all types of reads (template, complement and 2D), so check out the option `--type`. Be careful, the output `.fastq` file is written in the standard output, so you need to save the standard output to a file using `>`. For example: 
``` sh
poretools fastq [...] > lambda_minion/all_reads.fastq
## And even better you can add a step to compress the output file:
gzip -9 lambda_minion/all_reads.fastq
```

![Warning](img/warning.png)
**Remember, that it is not allowed to launch your command on the frontal machine!** You need to write the above commands in a bash script (called for example `minion_extract.sh`) that you will submit to the cluster with `bsub` (see [Submitting commands to the cluster](#submitting-commands-to-the-cluster)):
```sh
module add UHTS/Analysis/poretools/0.5.1
bsub < minion_extract.sh
```

<!--
bsub -q priority -o minion_extract.out -e minion_extract.err -J minion_extract 'poretools fastq lambda_minion/fast5 | gzip -9 > lambda_minion/all_reads.fastq.gz'
For 2D reads only add: --type 2D 
-->

![Question](img/round-help-button.png)
Have a look at the `.fastq` file created, using `less` or `nano`. Looking at header lines (starting with a `@`), do you identify which reads are 2D, template or complement reads? Do you identify which lines correspond to the quality scores (their header is a `+` sign).

Focus on a 2D read chosen at random. Look for the corresponding template and complement reads, which should be the next ones in the file (or you can search for the unique hash at the beginning of the header, which is common between 2D, template and complement reads).

![Question](img/round-help-button.png)
Do the quality scores seem to be improved in 2D reads? You can refer to [this wikipedia page](http://en.wikipedia.org/wiki/FASTQ_format#Encoding) for help on the PHRED quality scores in `.fastq` files. 

<!--
TO DO? 
Poretools stats
Would be nice to know what we are dealing with
-->

![Tip](img/elemental-tip.png)
We will not do this today, but for basic quality control of the reads, you can also launch the `fastqc` software (widely used for Illumina data) on this `fastq` file. The reported statistics are correct, just keep in mind that warning flags in the report are meaningful for short reads, and sometimes not very informative for long reads.

![help](img/help.png) If you are lost, you can get extracted MinION reads by executing the following script. While the script is running, have a look at it to understand what is done.
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/1_Extraction_MinION.sh
```
<!--
TO DO? 
- compress fastq file
- extract all types of reads
-->

### Extraction of PacBio RS II reads
First, create a working directory for raw reads and create links to the raw sequencing files.
```sh
mkdir -p lambda_RSII/raw_reads/
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/RSII_lambda_reads/*.h5 lambda_RSII/raw_reads/
```

![Question](img/round-help-button.png)
How many SMRT cells do we have? 
<!--
Answer: 1
-->

![To do](img/wrench-and-hammer.png)
Extract PacBio subreads using `bash5tools.py` from `pbh5tools`. Again, refer to the [documentation](https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst). The commands will look like this:
```sh
bash5tools.py <file.bas.h5> --outFilePrefix lambda_RSII/<prefix_for_extracted_reads> --readType subreads --outType fastq
## compress the output file
gzip -9 lambda_RSII/[prefix.fastq]
```
Put the full command in a script and submit it using bsub:
```sh
module add UHTS/PacBio/pbh5tools/0.8.0
bsub < pacbio_extract.sh
```
<!--
bsub -q priority -o RSII_extract.out -e RSII_extract.err -J RSII_extract 'bash5tools.py lambda_RSII/raw_reads/m140715_200214_42182_c100661412550000001823125411271451_s1_p0.bas.h5 --outFilePrefix pacbio --readType subreads --outType fastq'
-->

![Question](img/round-help-button.png)
How many subreads were produced? 
<!--
wc -l <prefix_for_extracted_reads>.fastq   # reads in fastq = number of lines / 4
Answer: 132269
-->

![help](img/help.png) If you are lost, you can get extracted PacBio reads by executing
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/1_Extraction_RSII.sh
```

## 2. Genome assembly
Genome assembly with long reads is, just like the established assembly with short-reads, based on graph theory. But string graphs are used instead of De Brujin graphs, and instead of kmer decomposition, the graph of reads (nodes) and their overlaps (edges) is constructed.

### Tools

#### Canu
Canu is an assembler for noisy long-reads sequences. It is a fork of Celera with added modules for error correction and trimming of reads. Canu was designed for both MinION and PacBio data, and is one of the two leading assemblers for longs reads (together with Falcon). Canu has many parameters, which are not discussed on this workshop, but many details are given in the [manual](http://canu.readthedocs.io/en/stable/).

#### Miniasm
Miniam is lightweight assembler. It very fast, but for a cost of simplicity and low parametrization. It can used as a first proxy of the data content, but for a final assembly, another assembler should be considered. As a first step, the overlap of reads is computed in a separated step using a standalone program called `minimap`.

### MinION
We will assemble the lambda phage genome using only the 2D reads. 1D reads are of substantially lower quality, so we would like to avoid using them if we have enough 2D reads. Since we do (~3,000 reads, i.e., coverage of ~500X), the assembly should be trivial: the genome is 48.8 kb and we have some reads of 20kb. We expect only one contig!

#### Canu
![To do](img/wrench-and-hammer.png) Choose the appropriate parameters to run [Canu](http://canu.readthedocs.io/en/latest/commands/canu.html):
```
usage: canu -p <assembly-prefix>    # use a specified prefix to name canu output files
            -d <assembly-directory> # use a specified prefix for output directories
            genomeSize=<number>     # expected genome size
   				  useGrid=false           # disables automatic submission to cluster
            -nanopore-raw           # specifies that ONT MinION data are used
            <data in fastq format>
```

Put the full command in a script and submit it using bsub:
```sh
module add UHTS/Assembler/canu/1.3;
bsub < minion_assembly_canu.sh
```

<!--
bsub -q priority -o minion_canu.out -e minion_canu.err -J minion_canu 'canu -p lambda -d lambda_minion/canu/ genomeSize=49k -nanopore-raw lambda_minion/all_reads.fastq.gz useGrid=false'
-->

The output directory contains many files. The most interesting ones are:
- *.correctedReads.fasta.gz : file containing the input sequences after correction, trim and split based on consensus evidence. 
- *.trimmedReads.fastq : file containing the sequences after correction and final trimming
- *.layout : file containing informations about read inclusion in the final assembly
- *.gfa : file containing the assembly graph by Canu
- *.contigs.fasta: file containing everything that could be assembled and is part of the primary assembly

![Question](img/round-help-button.png) How many contigs were produced? Does the total size seem to match your expectations?
<!--
Answer: 1, yes
->

![Warning](img/warning.png)
During our tests, the assembly did not always work. Sometimes the job was killed by the cluster, please check carefully the standard output and error files. We think that the memory requirement of Canu might be too big for some machines of the cluster. If this happens to you, rerun the job using the queue `dee-hugemem` instead of the queue `priority`. It night be useful to increase the amount of memory requested with the options `-v 20000000 -R "rusage[swp=20000]"` and `-M 10000000 -R "rusage[mem=10000]"`. If this still does not work, you can consult a successfull run that we pre-computed in the folder [...] 

**TO DO (see `/scratch/beegfs/monthly/jroux/tp_long_reads/lambda_minion/canu/`) + check that submission on this queue works with student accounts**.

<!--
bsub -q dee-hugemem -o minion_canu_hugemem.out -e minion_canu_hugemem.err -J minion_canu_hugemem -v 20000000 -R "rusage[swp=20000]" -M 10000000 -R "rusage[mem=10000]" 'canu -p lambda -d lambda_minion/canu/ genomeSize=49k -nanopore-raw lambda_minion/all_reads.fastq.gz useGrid=false'
-->

<!-- 
For now I removed this part since it is done for pacbio below

#### Bonus: Miniasm
![To do](img/wrench-and-hammer.png) You first need to compute the overlaps of reads using `minimap`. As recommended in the [documentation](https://github.com/lh3/miniasm), put the following command and parameters in your submission script:

```sh
minimap -S -w 5 -L 100 -m 0 <reads.fq> <reads.fq> | gzip -9 > <overlaps.gz>
```

<!--
bsub -q priority -o minion_miniasm.out -e minion_miniasm.err -J minion_miniasm 'mkdir lambda_minion/miniasm/; minimap -S -w 5 -L 100 -m 0 lambda_minion/all_reads.fastq.gz lambda_minion/all_reads.fastq.gz | gzip -9 > lambda_minion/miniasm/overlaps.gz'
-->

Then, to perform the assembly using `miniasm`:
```sh
miniasm -f <reads.fq> <overlaps.gz> > <contigs.gfa>
## convert .gfa to .fasta
awk '/^S/{print ">"$2"\n"$3}' <contigs.gfa> | fold > <contigs.fa>
```

<!--
bsub -q priority -o minion_miniasm.out -e minion_miniasm.err -J minion_miniasm 'miniasm -f lambda_minion/all_reads.fastq.gz lambda_minion/miniasm/overlaps.gz > lambda_minion/miniasm/contigs.gfa'
awk '/^S/{print ">"$2"\n"$3}' lambda_minion/miniasm/contigs.gfa | fold > lambda_minion/miniasm/contigs.fa
-->

Put these commands in a script and submit it using bsub:
```sh
module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;
bsub < minion_assembly_miniasm.sh
```

![Tip](img/elemental-tip.png) The size of the `.fasta` file in bytes corresponds to the number of nucleotides, newline characters and characters in headers. This is useful to roughly estimate the length of an assembly.

![help](img/help.png) If you are lost, you can get the assembly of MinION reads by executing:
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/2_MinION_Assembly.sh
```

**Miniasm assembly of MinION reads is only 29kb -> only use miniasm with pacbio (since canu doesn't give an assembly of the rigth size)**
-->

### PacBio RS II
One PacBio RS II SMRT cell produces between 0.5 and 1 gb. This roughly correspond to a 10,000x coverage of the lambda phage genome which is clearly an overkill. We can decrease the computational load by using a subset of only a few thousand reads for the assembly.

![To do](img/wrench-and-hammer.png) 
Go to the directory including the extracted PacBio reads and extract a subset of the first 3,000 reads.
```sh
zcat RSII_reads.fastq.gz | head -12000 | gzip -9 > RSII_reads_subset.fastq
```

The assembly step is now analogous to the assembly of MinION data

![To do](img/wrench-and-hammer.png) Modify the command performing `Canu` assembly for PacBio data. Change the type of input reads to `-pacbio-raw`, the name of the input file and the name of output folder, and launch the assembly.

<!--
module add UHTS/Assembler/canu/1.3;
bsub -q dee-hugemem -o RSII_canu_hugemem.out -e RSII_canu_hugemem.err -J RSII_canu_hugemem -v 20000000 -R "rusage[swp=20000]" -M 10000000 -R "rusage[mem=10000]" 'canu -p lambda -d lambda_RSII/canu/ genomeSize=49k -pacbio-raw lambda_RSII/RSII_reads_subset.fastq useGrid=false'
Gives one contig of length 58.7kb

But 3000 subreads do not correspond to 3000 unique sequences... Maybe we should take more subreads as input?
awk 'NR%4==1' RSII_reads_subset.fastq | cut -f2 -d"/" | sort | uniq | wc -l ## 1735
head -n 22000 lambda_RSII/pacbio.fastq | awk 'NR%4==1' | cut -f2 -d"/" | sort | uniq | wc -l ## 3156
Still contig of 59kb :(

head -40000 lambda_RSII/pacbio.fastq | gzip -9 > lambda_RSII/RSII_reads_subset.fastq.gz
Still one contig, but 67.4kb :(
-->

![Question](img/round-help-button.png)
What is the length of the resulting assembly? Is there anything weird going on? 

<!--
yes, PacBio reads are not assembled well. Maybe too much coverage? Quality of the data?
-->

#### Bonus: Miniasm
Since the assembly with Canu does not seem to be successful, we will try a different assembler, `miniasm`.

![To do](img/wrench-and-hammer.png)
You first need to compute the overlaps of reads using `minimap`. As recommended in the [documentation](https://github.com/lh3/miniasm), put the following command and parameters in your submission script:

```sh
minimap -S -w 5 -L 100 -m 0 <reads.fq> <reads.fq> | gzip -9 > <overlaps.gz>
```
<!--
module add UHTS/Analysis/minimap/0.2.r124.dirty;
bsub -q priority -o RSII_miniasm.out -e RSII_miniasm.err -J RSII_miniasm 'mkdir lambda_RSII/miniasm/; minimap -S -w 5 -L 100 -m 0 lambda_RSII/RSII_reads_subset.fastq lambda_RSII/RSII_reads_subset.fastq | gzip -9 > lambda_RSII/miniasm/overlaps.gz'
-->

![To do](img/wrench-and-hammer.png)
Then, to perform the assembly using `miniasm`:
```sh
miniasm -f <reads.fq> <overlaps.gz> > <contigs.gfa>
## convert .gfa to .fasta
awk '/^S/{print ">"$2"\n"$3}' <contigs.gfa> | fold > <contigs.fa>
```

<!--
module add UHTS/Analysis/miniasm/0.2.r137.dirty;
bsub -q priority -o RSII_miniasm.out -e RSII_miniasm.err -J RSII_miniasm 'miniasm -f lambda_RSII/RSII_reads_subset.fastq lambda_RSII/miniasm/overlaps.gz > lambda_RSII/miniasm/contigs.gfa'
awk '/^S/{print ">"$2"\n"$3}' lambda_RSII/miniasm/contigs.gfa | fold > lambda_RSII/miniasm/contigs.fa
-->

Put these commands in a script and submit it using bsub:
```sh
module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;
bsub < RSII_assembly_miniasm.sh
```
![Tip](img/elemental-tip.png) The size of the `.fasta` file in bytes corresponds to the number of nucleotides, newline characters and characters in headers. This is useful to roughly estimate the length of an assembly. Ortherwise you can remove the header with tail and count the number of nucleotides: `tail -n+2 lambda_RSII/miniasm/contigs.fa | wc -m`

![Question](img/round-help-button.png)
What do you think about the length of this assembly?

![help](img/help.png) If you are lost, you can get the assembly based on PacBio reads by executing:
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/2_RSII_Assembly.sh
```

## 3. Quality of assemblies

`Canu` produces a `.html` report for every performed step. You can find there the amount of reads after correction, filtering, etc. Miniasm, however does not provide such a report, therefore we will use `Quast` to estimate basic statistics of the all assemblies in a comparable format. Then we will map the assemblies to the lambda phage reference genome to inspect their quality.

### Tools

#### Quast
Quast is a tool to evaluate genome assemblies. It provides a PDF report providing metrics on contigs. It can also be used for comparison between assemblies, or to compare a *de novo* assembly to a reference genome. Here are the [paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/), the [GitHub page](https://github.com/ablab/quast) and the [manual](http://quast.bioinf.spbau.ru/manual). 

**TO DO: add MUMmer**
[Documentation](http://mummer.sourceforge.net/manual/)

### Downloading the lambda phage reference genome
```sh
mkdir reference/
cd reference/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
ln -s NC_001416.fna lambda_ref_genome.fa # a clearer name
```

### Running Quast
The command in your submission script should look something like this:
```sh
quast.py -R <path to reference genome> <path to Canu MinION assembly> <path to Canu PacBio assembly> <path to Miniasm PacBio assembly> 
```

![To do](img/wrench-and-hammer.png)
You can submit your script after loading the appropriate module: 
```sh
module add UHTS/Quality_control/quast/4.1
bsub < quast.sh
```

<!--
bsub -q priority -o quast.out -e quast.err -J quast 'quast.py -R reference/lambda_ref_genome.fa lambda_minion/canu/lambda.contigs.fasta lambda_RSII/canu/lambda.contigs.fasta lambda_RSII/miniasm/contigs.fa'

bsub -q priority -o quast.out -e quast.err -J quast 'quast.py -R reference/lambda_ref_genome.fa lambda_minion/canu/lambda.contigs.fasta lambda_minion/miniasm/contigs.fa lambda_RSII/canu/lambda.contigs.fasta lambda_RSII/miniasm/contigs.fa'
-->

The output is a directory, typically `quast_results/<date_of_launch>`, containing several files including a report in PDF format. Copy it to your local machine using `scp`: `scp username@prd.vital-it.ch:[path_to_file_on_vital_it] [path_to_file_on_laptop]`

![Question](img/round-help-button.png)
What do you notice in this report? What did go wrong with the PacBio Canu assembly? Is the Miniasm assembly better? Compare the Canu assemblies made wiht both MinION and PacBio data: hwo does the error and indel rates compare?

<!--
Answers: 
pacbio canu assembly has a duplicated part (duplication ratio = 1.2)
Miniasm assembly worse: not even mapping to ref !!!
Minion has a lot more errors and indels!!!
-->

### Bonus: getting further
It would be interesting to know more precisely what is going on with PacBio assemblies. 
* The software `Mauve` (see [website](http://darlinglab.org/mauve/user-guide/introduction.html)) could be interesting, but it is GUI only, so you need to install it on your laptop.  
* It is probably a good idea to generate a dot plot to visualize and compare assemblies. This can be done with the `MUMmer` software. The command below will find all maximal unique matches (-mum) between the assemblies on both the forward and reverse strands (-b) and report all the match positions relative to the forward strand (-c). The output lists all of the matches between the two input sequences.

```sh
mkdir mummer/
mummer -mum -b -c <assembly_1.fa> <assembly_2.fa> > mummer/<assembly_1_vs_2.mums>
```

<!--
mkdir mummer
bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummer -mum -b -c reference/lambda_ref_genome.fa lambda_minion/canu/lambda.contigs.fasta > mummer/lambda_minion_canu.mums'

bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummer -mum -b -c reference/lambda_ref_genome.fa lambda_minion/miniasm/contigs.fa > mummer/lambda_minion_miniasm.mums'

bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummer -mum -b -c reference/lambda_ref_genome.fa lambda_RSII/canu/lambda.contigs.fasta > mummer/lambda_RSII_canu.mums'

bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummer -mum -b -c reference/lambda_ref_genome.fa lambda_RSII/miniasm/contigs.fa > mummer/lambda_RSII_miniasm.mums'
-->

To generate a dotplot of all the matches between two sequences:
```sh
mummerplot -postscript -p <assembly_1_vs_2> <assembly_1_vs_2.mums>
```
<!--
bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummerplot -postscript -p mummer/lambda_minion_canu mummer/lambda_minion_canu.mums'
bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummerplot -postscript -p mummer/lambda_minion_miniasm mummer/lambda_minion_miniasm.mums'
bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummerplot -postscript -p mummer/lambda_RSII_canu mummer/lambda_RSII_canu.mums'
bsub -q priority -o mummer.out -e mummer.err -J mummer 'mummerplot -postscript -p mummer/lambda_RSII_miniasm mummer/lambda_RSII_miniasm.mums'
-->

![To do](img/wrench-and-hammer.png)
Write a script to perform the comparison of all assemblies to the reference genome. You can submit your script after loading the appropriate module: 
```sh
module add UHTS/Analysis/MUMmer/3.23;
bsub < mummer.sh
```

![Question](img/round-help-button.png)
Download the postscript dot plots to your computer with `scp`. What do you think of the dot plots? Check this [cheat sheet](http://mummer.sourceforge.net/manual/AlignmentTypes.pdf) for interpreting dot plots. Is the problem with the PacBio Canu assembly really important? Can you guess from this plot what went wrong when passing the Miniasm assembly through Quast?

<!--
Answers: Canu pacbio: overlap between beginning and end of assembly: not so bad!
Miniasm quite bad at correcting errors: probably why did not map with Quast

TO DO?
![help](img/help.png)
If you are lost, you can get all quality stats of all assemblies by executing:
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/_Assembly_quality.sh
```
-->

## 4. Bonus: mapping reads to a reference genome

### Tools

#### NanoOK
This tool allows to perform multiple analyses over MinION data, including the extraction of read sequences from `.fast5` files, their alignment to a reference, and the generation of a summary report of QC and mapping statistics. This software was designed for MinION data, but it is easy to "hack" it to use PacBio reads. Source code can be found on [GitHub](http://github.com/TGAC/NanoOK) and software usage is detailed in the [documentation](http://documentation.tgac.ac.uk/display/NANOOK/NanoOK).

#### LAST
LAST finds similar regions between sequences (local alignment). It is faster than BLAST and can handle largeer datasets. See the project [website](http://last.cbrc.jp/).

### MinION

To map the reads to the reference genome of the lambda phage using the `LAST` aligner, the reference sequence first needs to be indexed. The `-Q 0` tells LAST to expect a fasta file.

```sh
cd reference/
module add UHTS/Analysis/NanoOK/0.72
lastdb -Q 0 lambda_ref_genome lambda_ref_genome.fa
cd ..
```
<!--
This is very fast, so probably no need to bsub...
-->

![To do](img/wrench-and-hammer.png)
Now launch the alignment with `LAST` (default aligner) using the `nanook align` utility. NanoOK expects a fasta file for each read, in a folder called `fasta/`, so we need a bit of file manipulation first. The commands in your submission script should look like this:

```sh
mkdir lambda_minion/fasta/

# Extract the read headers and sequences only
zcat lambda_minion/all_reads.fastq.gz | awk 'NR%4==1||NR%4==2' \
  > lambda_minion/all_reads.fasta 

# Convert header according to fasta specifications
perl -p -i -e 's/^@/>/g' lambda_minion/all_reads.fasta 

# split each read in a different file
nanook_split_fasta -i lambda_minion/all_reads.fasta -o lambda_minion/fasta/

# move reads to folders reflecting their type: 2D, template, complement
mkdir lambda_minion/fasta/2D
mkdir lambda_minion/fasta/Complement
mkdir lambda_minion/fasta/Template
find lambda_minion/fasta/ -maxdepth 1 -name \*_2d.fasta \
  -exec mv {} lambda_minion/fasta/2D/ \;
find lambda_minion/fasta/ -maxdepth 1 -name \*_template.fasta \
  -exec mv {} lambda_minion/fasta/Template/ \;
find lambda_minion/fasta/ -maxdepth 1 -name \*_complement.fasta \
  -exec mv {} lambda_minion/fasta/Complement/ \;

# Launch the alignment
nanook align -s lambda_minion -r reference/lambda_ref_genome.fa
```

![Question](img/round-help-button.png)
What directories have been created? What do they contain? Look at one randomly chosen alignment file. What is striking?

![Tip](img/elemental-tip.png)
It is possible to specify to `nanook align` the alignment parameters to be used with `LAST`, with the `-alignerparams` option. By default the parameters are `-s 2 -T 0 -Q 0 -a 1`.

![Tip](img/elemental-tip.png)
You can try to align reads through `NanoOK` with other aligners. For example to use `BWA-MEM`, you first need to create the index of the reference sequence (`bwa index reference.fasta`). Then you can relaunch `nanook align` with the `-aligner bwa` option.

### Statistics and QC report
![To do](img/wrench-and-hammer.png)
Launch the generation of the NanoOK alignment report including QC and alignment statistics using the ```nanook analyse``` utility. A PDF file should be created in the ```latex_last_passfail/``` folder. If it's not the case, you can generate it yourself with the ```pdflatex file.tex``` command (press enter every time the program prompt you with a question).

```sh 
nanook analyse -s lambda_minion -r reference/lambda_ref_genome.fa
pdflatex lambda_minion/latex_last_passfail/lambda_minion.tex
```

![Question](img/round-help-button.png)
Download the report to your laptop, and take some time to read and understand it. Here are a few questions that will guide you:
* What size is the longest template read? Is that surprising?
* What does the N50 values indicate us? Is that consistent with the size of the DNA fragments used to create the MinION library (~8kb)?
* What is most common error type within reads?
* In comparison, the typical error rate reported for the Illumina sequencing technology is around 0.1%. For Sanger sequencing, it can be as low as 0.001%... An interesting comparison of error rates across platforms was recently [published](http://bib.oxfordjournals.org/content/17/1/154).
* Why is the alignment rate of 2D reads higher than those of template and complement reads?
* Which are the most accurate: shortest or longest reads?
* Have a look at the coverage plot. There is a "bump" in coverage around 45kb. This corresponds to some control spike-in DNA that was added during library preparation (more precisely around 3.6kb of a region of the lambda phage genome, with a single mutation G45352A). Can additional variation be explained by GC content differences?
* Have a look at the k-mer over and under-representation analysis. What sort of k-mers are under-represented in 2D reads? Is that expected given how the technology works?

<!-- 
TO DO: tell them to have a look at alignments in this region to see the mutation?
-->

### PacBio RS II
It is quite easy to redo the above steps with PacBio data. Do it if you have time ;-) 

<!--
TO DOs
* check this paper https://dx.doi.org/10.7554/eLife.14258
* check the porecamp material (analysis part): http://porecamp.github.io/timetable.html
* Check Pradervand talk: http://edu.isb-sib.ch/pluginfile.php/3390/course/section/1574/Pradervand_Sequencing_CUSO2015.pdf
* Will these modules be needed?
```sh
module add UHTS/PacBio/blasr/20140829;
module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
module add UHTS/Analysis/seqtk/2015.10.15 # fastq to fasta
```
* Idea: make them download maf files and reference genome, and open in IGV... Or is it possible to generate a picture on vital-IT from IGV?

* TO DO: export final HTML and send to Patricia for testing
Fork to SIB github
https://github.com/sib-swiss/2016-07-05-longreads-bern 

* TO DO: cut lines of code that are too long

* TO DO: git pull of scratch/beegfs/monthly/SIB_long_read_workshop/

* TO DO: remove gzip commands?

* TO DO: change "we to "you"

* TO DO: test dee-hugemem with student account

![Question](img/round-help-button.png)
![Tip](img/elemental-tip.png)
![To do](img/wrench-and-hammer.png)
![Warning](img/warning.png)
-->

