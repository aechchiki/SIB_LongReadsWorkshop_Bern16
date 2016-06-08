# MinION practical

# data location 
ls /archive/dee/robinson/jroux/MinION/Minion_Lambda_Library
# data already in .fastq


### read extraction

# -
# Q: I have .fast5, how to convert to .fastq? 
# A: Use poretools on basecalled reads.
#
# install poretools:
git clone https://github.com/arq5x/poretools
# check for write rights
cd poretools; python setup.py install --root
# example
# cd /your_dir/MinION/basecalled_reads
# /home/user/bin/poretools/poretools/fastq.py input.fast5 > output.fastq
#
# new! poretools is now on Vital-IT (thx Sébastien)
# module add UHTS/Analysis/poretools/0.5.1

# -
# Q: I want one big fastq containing all reads
# A: Use direct globbing
#
# cd /your_dir/MinION/basecalled_reads
module add UHTS/Analysis/poretools/0.5.1 # load poretools
poretools fastq input.fast5 > output.fastq

# -
# Q: I want one fastq per fast5 input
# A: Use bash loop instead
#
# cd /your_dir/MinION/basecalled_reads
a=0; for i in $(ls *.fast5); do echo $i;a=$(echo $i | cut -d'.' -f1); echo $a; b=$((b + 1)); module add UHTS/Analysis/poretools/0.5.1; poretools fastq $i > $a'.fastq'; done

# -
# Q: How many resources does it take to convert fast5 to fastq?
# A: For an archive of 256M (fast5), CPU time: 12.19 sec; Max Memory: 29.20 MB

# -
# Q: How to check the length of each read sequenced by MinION?
# A: extract raw sequence from fastq, then count the number of characters per line
#
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
less LambdaBurnIn.2D.fastq | grep -E '^[ACTGN]+$' | while read rawseq; do echo -n "$rawseq" | wc -c ; done > readlength.txt

# -
# Q: How to check some stats on the length of the reads sequenced by MinION?
# A: use extracted raw sequences character count, then calculate the some stats
#
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
awk '{ total += $1; count++ } END { print total/count }' readlength.txt # mean read length
sort -nk 1 readlength.txt | head -n 1 # min read length
sort -nrk 1 readlength.txt | head -n 1 # max read length
less LambdaBurnIn.2D.fastq | grep ^@ | grep .fast5$ | wc -l # how many reads in the input
#
# OR, more infos with poretools stats on fast5 directly
# cd /your_dir/MinION/basecalled_reads
# module add UHTS/Analysis/poretools/0.5.1; poretools stats input.fast5 > stats_output.txt


### read alignment to reference

# -
# Q: How to align MinION reads to known reference?
# A: Can use BWA (faster but less sensitive), or LAST (slower but more sensitive)

# part 1: BWA
#
# genome indexing 
module add UHTS/Aligner/bwa/0.7.13 # load bwa
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome directory
bwa index Lambda.fasta # create index
# 
# create fasta index for samtools
module add UHTS/Analysis/samtools/1.3 # load samtools
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome directory
samtools faidx Lambda.fasta # create index
# 
# run bwa, pipe to samtools 
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
bwa mem -x ont2d /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda.fasta LambdaBurnIn.2D.fastq | samtools view -T /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda.fasta -bS - | samtools sort -T Lambda.bwa -o Lambda.bwa.bam -
# options: -x define read type

# part 2: LAST
# 
# genome indexing
module add SequenceAnalysis/SequenceAlignment/last/531 # load last 
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome directory
lastdb Lambda Lambda.fasta  # create index
# 
# convert fastq to fasta for last
# note: I tried it on fastq with -Q1 option (def fastq-sanger format) but got error 
module add UHTS/Analysis/seqtk/2015.10.15
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq # move to genome
seqtk seq -a LambdaBurnIn.2D.fastq > Lambda.fasta 
# 
# run lastal on fasta reads
lastal -q 1 -a 1 -b 1 /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda Lambda.fasta > Lambda.maf
# options: -q = mismatch cost; -a = gap cost; -b = gap extension cost; -Q1 = read input in sanger-fastq format 
#
# get nanopore scripts 
# cd /home/aechchik/bin/; git clone https://github.com/arq5x/nanopore-scripts.git
# 
# convert maf to bam with complete CIGAR (matches and mismatches)
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
module add UHTS/Analysis/samtools/1.3 # load samtools
python /home/aechchik/bin/nanopore-scripts/maf-convert.py sam Lambda.maf | samtools view -T /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda.fasta -bS - | samtools sort -T Lambda.last -o Lambda.last.bam -




