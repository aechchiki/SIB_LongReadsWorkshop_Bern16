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

You will go through different steps, which include the extraction of reads from their native encoding formats ([HDF5](http://en.wikipedia.org/wiki/Hierarchical_Data_Format)), quality control, *de-novo* genome assembly, and a mapping of reads to a reference genome.

After DNA fragmentation, the **MinION** library is prepared by ligating adapters to the double-stranded DNA fragments. One side receives a Y-form adapter, and the other side a hairpin-form adaptor. The adapters are conjugated with motor proteins that help control the translocation speed of DNA through the pore. When the Y-form adapter approaches the pore, one strand starts to be sequenced. This sequence is called *template*. After the hairpin-form adapter is sequenced, the *complement* sequence is read. The consensus of template and complement sequences is called a *2D read*. An illustration of the protocol is provided in this [figure](http://www.genetics.org/content/202/1/37):

![protocol](img/protocol.jpg)

For **PacBio** sequencing (or **SMRT** sequencing), hairpin-form adapters (SMRTbell) are ligated on both sides of double-stranded DNA fragments, forming a circular sequence. Circular sequences attach to the bottom of theSMRT cell wells, and are continuously read, and a movie is recorded. Nowadays, the total length of movie allows to record up to 75000 nucleotide basecalls. This implies that short DNA fragments will be read many times, while long fragments may be read just one or twice. The complete sequence is called a *polymerase* read. After the adapter sequences are removed, the sequence is split into *subreads*. The consensus of subreads is called a *circular consensus read* or *read of insert*.

![protocol](img/pb_reads.png)

<!--
In the MinION experiment you are going to analyze, the lambda phage DNA was fragmented using sonication (Covaris) to yield ~8kb-long fragments. One microgram of fragmented DNA material was then used for library preparation. Other protocols exist for smaller amounts of starting material. 

The sequencing run produced a single ```.fast5``` file per pore. All files were then uploaded to the cloud for basecalling. The base-called files were downloaded back to vital-IT. They are also in the ```.fast5``` format, but there is one file per read.

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

You are in! Jump to: "Setting up your working directory".

### For Windows users

You should first install a ssh client (e.g., `PuTTY`). If you already have one, we assume you know how to use it. Connect to Vital-IT and jump to: "Setting up your working directory".

If you do not have a ssh client, follow [these steps](https://github.com/aechchiki/SIB_LongReadsWorkshop_Bern16/blob/master/vital-it_connect_Putty.pdf). 

### Setting up your working directory

* On Vital-IT, you must read and write files in the "scratch" directory. Move to `/scratch` weekly directory:
```sh
cd /scratch/beegfs/weekly/
```

![Warning](img/warning.png)
Be careful, files in this folder are erased after a week. At the end of the practicals, if you want to keep your results, you need to back-up the data, for example in a compressed tarball that you move to your home or archive folder.

* Create your own directory:
```sh
mkdir <username> ; cd <username>
```

* You will always be working from this directory. Before launching commands please be sure that you are located in the right directory by typing: ```pwd```. The expected output is: ```/scratch/beegfs/weekly/<username>```

### Submitting commands to the cluster

It is not allowed to launch any calculation on the frontal machine where you are connected (```prd.vital-it.ch```). You need to submit each job for batched execution through a job scheduler that will dispatch it on the cluster nodes. For example:
```sh
bsub '[my command line here]'
bsub < script.sh
```

TO DO (walid): list useful bsub commands options:
* priority queue `-q priority`
* Memory usage? `-M / -R`
* How to get output and error in files instead of email: `-o / -e`

## 1. Read extraction

The output of MinION and PacBio RS II is stored in Hierarchical Data Format, what is basically an archive format (like `.zip` or `.tar`), but with very quick access to its content (You can find details on [wikipedia](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)). In those files are more details about reads and basecalling. What, we need is just a simple `.fastq` for downstream analysis.

![Tip](img/elemental-tip.png)
The files saved in Hierarchical Data Format can be explored using HDF5 command-line tools:

```sh
h5dump <HDF_file> | less
h5stat <HDF_file> | less
```

**fast5**

MinION basecaller produces one file per read. In the file there is a time series of a current used for basecalling, basecalled template and complement reads (also called 1D reads) and consensus of 1D reads called 2D reads with comparably higher accuracy.

**bas.h5 & bax.h5**

The RS II system produces three `bax.h5` files and a `bas.h5` per SMRT cell. The three `bax.h5` files correspond to the first, second and third part of the movie capturing the SMRT cell, `bas.h5` contains metainformations. PacBio announced a change of data format to specialised `.bam` for platform Sequel, however all data produced by RS II will be still in `.h5` formats we are going to work with.

### Tools

All used tools are installed already on Vital-IT. Just before the first usage, you will have a command which loades the package.

**NanoOK** allows to perform multiple analyses over MinION data, including the extraction of read sequences from `.fast5` files, their alignment to a reference, and the generation of a summary report of QC and mapping statistics. Note: this software was designed for MinION data, but it is easy to hack to use PacBio reads. Source code can be found on [GitHub](http://github.com/TGAC/NanoOK) and software usage is detailed in [documentation](http://documentation.tgac.ac.uk/display/NANOOK/NanoOK).

**pbh5tools** is a set of python scripts to extract `.fasta` and `.fastq` from `bas.h5` and `bax.h5` files. These scripts allow filtering based on error rate, read length and read type. Source code is available on [GitHub](https://github.com/PacificBiosciences/pbh5tools) and software usage is detailed in [documentation](https://github.com/PacificBiosciences/pbh5tools/blob/master/doc/index.rst). 

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

![To do](img/wrench-and-hammer.png)
You will need to convert the MinION raw reads from their `.fast5` to a more classical and readable `.fastq`. This can be done with the `nanook extract` utility, and takes around 10 minutes.

```sh
module add UHTS/Analysis/NanoOK/0.72
bsub -q priority 'nanook extract -fastq -s lambda_minion'
```
TO DO: shouldn't it be rather this?
```sh
bsub -q priority 'module add UHTS/Analysis/NanoOK/0.72; nanook extract -fastq -s lambda_minion'
```

![To do](img/wrench-and-hammer.png) Cat all to one file... (TO ADD)

![Question](img/round-help-button.png)
Do the quality scores seem to be improved in 2D reads? 

![Tip](img/elemental-tip.png)
To access quality values, you can use `fastqc` software, as for short-reads data. Reported stats are correct, just keep in mind that warning flags in the report are for assembly by short reads and therefore not very informative.

![help](img/help.png) If you are lost, you can get extracted MinION reads by executing
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/1_Extraction_MinION.sh
```

### Pacbio RS II
```sh
module add UHTS/PacBio/pbh5tools/0.8.0;
```

![To do](img/wrench-and-hammer.png) Extract also RS II reads using 

```sh
mkdir -p lambda_RSII/raw_reads/
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/RSII_lambda_reads/*.h5 lambda_RSII/raw_reads/
```

![Question](img/round-help-button.png)
How reads of many SMRT cells we have? [1]

![To do](img/wrench-and-hammer.png) Extract RS II reads using `bash5tools.py` from `pbh5tools`.

```sh
bsub -q priority bash5tools.py <file.bas.h5> --outFilePrefix <prefix_for_extracted_reads> --readType subreads --outType fastq
```

![Question](img/round-help-button.png)
How many reads has was produced by the SMRT cell? [132269]


```sh
wc -l <prefix_for_extracted_reads>.fastq   # reads in fastq = number of lines / 4
```

![help](img/help.png) If you are lost, you can get extracted RS II reads by executing

```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/1_Extraction_RSII.sh
```

***

## 2. Genome assembly

Long reads genome assembly is just like established short-reads assembly based on graph theory, but string graphs are used instead of De Brujin. Instead of kmer decomposition, the graph of reads (nodes) and their overlaps (edges) is constructed.

There are two leading assemblers for long reads Canu and Falcon. Canu is a fork of Celera with added modules for error correction and trimming of reads. Falcon is an assembler aware of ploidy developed from scratch for PacBio data.

### Tools

We will use Canu and Miniasm for both MinION and RS II assembly, then we check the assembly quality using the report provided by Quast. 

**Canu** is an assembler for noisy long-reads sequences. Canu will correct the reads, trim suspicious regions, then assemble the corrected and cleaned reads into unitigs. Note: this software was designed for both MinION and RS II data. Canu has many parameters, which are not discussed on this workshop, however all of them can be found in very detailed [manual](http://canu.readthedocs.io/en/stable/) of canu.

**Miniasm** is lightweight assembler. It very fast, but for a cost of simplicity and low parametrization. It can used as a first proxy of the data content, but for the final assembly, another assembler should be considered. The first step, the overlap of reads is computed in a separated step using a standalone program called `minimap`.

### MinION

We assemble lambda phage using the 3,068 2D reads (cca 500x). 1D reads are in substantially worse quality, so we would like to avoid using them if we have enough 2D reads. Since we do, the assembly should be trivial: The genome is 48.8kbp and we have some reads of 20kb! We expect only one contig!

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

![Tip](img/elemental-tip.png) The size of the `.fasta` file in bytes is number of nucleotides, newlines and letters in all headers inside. Therefire you can roughly estimate the length of the assembly from the size of files.

![help](img/help.png) If you are lost, you can get both assemblies of MinION and RS II reads by executing
```sh
bsub < /scratch/beegfs/monthly/SIB_long_read_workshop/scripts/2_MinION_Assembly.sh
```

### Pacbio RS II

One SMRT cell produces yield between 0.5-1gbp. What roughly correspond to 10,000 to 20,000x coverage of lambda phage genome. As you can imagine, it is a bit overkill. We can a decrease a computational load very much by assembling a subset of a few thousand reads only.

![To do](img/wrench-and-hammer.png) Go to directory with extracted RS II reads and take a subset of the first 3,000 reads.

```sh
head -12000 RSII_reads.fastq > RSII_reads_subset.fastq
```

The assembly is now analogical to the assembly of MinION data

![To do](img/wrench-and-hammer.png) Modify the command performing `Canu` assembly for RS II data. Change the type of input reads to `-pacbio-raw`, name of the input file and the name of output name and perform the assembly.


![To do](img/wrench-and-hammer.png) Perform assembly of RS II data using miniasm as well.

![Question](img/round-help-button.png) What is are the lengths of all four assemblies? Is there anything weird going on? [yes, RS II reads are not assembled well.]


![help](img/help.png) If you are lost, you can get both assemblies of MinION and RSII reads by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/2_RSII_Assembly.sh
```

***

## 3. Quality of assembly

`Canu` produces a `.html` report for every performed step. You can find there amounts of reads after correction, filtering, etc. However, miniasm, does not therefore we use `Quast` for stats of all smeblies to get unified format. Then we will map assemblies to reference to find amount of expected differences.

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



![help](img/help.png) If you are lost, you can get all quality stats of all assemblies by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/3_Assembly_quality.sh
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
You can try to align reads with other aligners. For example to use `BWA-MEM`, you first need to load the corresponding module on vital-IT (`module add UHTS/Aligner/bwa/0.7.13`), create the index of the reference sequence (`bwa index reference.fasta`). Then you can relaunch `nanook align` with the `-aligner bwa` option.

![help](img/help.png) If you are lost, you can perform mapping of RS II reads by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/4_RSII_mapping.sh
```

### Pacbio RS II

![help](img/help.png) If you are lost, you can perform mapping of RS II reads by executing
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

### Pacbio RS II

![help](img/help.png) If you are lost, you can perform a quality control report of mapping by executing
```sh
bsub < /scratch/beefaskf/monthly/SIB_long_reads/4_mapping_qc_report.sh
```

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
* Idea: make them download maf files and reference genome, and open in IGV... Or is it possible to generate a picture on vital-IT from IGV?

![Question](img/round-help-button.png)
![Tip](img/elemental-tip.png)
![To do](img/wrench-and-hammer.png)
![Warning](img/warning.png)
-->
