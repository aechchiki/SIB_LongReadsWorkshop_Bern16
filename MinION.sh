# MinION practical

# fast5 location 
ls /home/jroux/archive/MinION/run_MinION_2015_06_15_Burn-In.tar
# extract to my /scratch
tar -xf /home/jroux/archive/MinION/run_MinION_2015_06_15_Burn-In.tar -C /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/data/

# load needed modules for the session 
module add UHTS/Analysis/poretools/0.5.1 # load poretools (read extraction)
# note: other nanopore tools are not installed in cluster, should I ask for them?
module add UHTS/Aligner/bwa/0.7.13 # load bwa (alignment)
module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
module add SequenceAnalysis/SequenceAlignment/last/531 # load last (alignment)
module add UHTS/Analysis/seqtk/2015.10.15 # fasatq to fasta


### read extraction

# -
# Q: I have .fast5, how to convert to .fastq? 
# A: Use poretools on basecalled reads.
#
# to install poretools:
# git clone https://github.com/arq5x/poretools
# cd poretools; python setup.py install --root
#
# NOW poretools is now on Vital-IT (thx Sébastien)
# module add UHTS/Analysis/poretools/0.5.1 # load poretools

# -
# Q: I want one big fastq containing all reads
# A: Use direct globbing
#
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/data/Downloads_burnin
poretools fastq *.fast5 > /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq/Lambda.fastq

# -
# Q: I want one fastq per fast5 input
# A: Use bash loop instead
#
# not advised
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/data/Downloads_burnin
a=0; for i in $(ls *.fast5); do echo $i;a=$(echo $i | cut -d'.' -f1); echo $a; b=$((b + 1)); module add UHTS/Analysis/poretools/0.5.1; poretools fastq $i > /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq/$a'.fastq'; done

# -
# Q: How many resources does it take to convert fast5 to fastq?
# A: Current archive contains 13G fast5. 0m43.505s CPU time

# -
# Q: How to check the length of each read sequenced by MinION?
# A: extract raw sequence from fastq, then count the number of characters per line
#
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
less Lambda.fastq | grep -E '^[ACTGN]+$' | while read rawseq; do echo -n "$rawseq" | wc -c ; done > readlength.txt

# -
# Q: How to check some stats on the length of the reads sequenced by MinION?
# A: use extracted raw sequences character count, then calculate the some stats
#
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
awk '{ total += $1; count++ } END { print total/count }' readlength.txt # mean read length
sort -nk 1 readlength.txt | head -n 1 # min read length
sort -nrk 1 readlength.txt | head -n 1 # max read length
less Lambda.fastq | grep ^@ | grep .fast5$ | wc -l # how many reads in the inputq
#
# OR, more infos with poretools stats on fast5 directly


### read alignment to reference

# -
# Q: How to align MinION reads to known reference?
# A: Can use BWA (faster but less sensitive), or LAST (slower but more sensitive)

# part 0: get ref genome 
# 
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
mv NC_001416.fna Lambda_RefGenome.fa

# part 1: BWA
#
# genome indexing for bwa
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome directory
bwa index Lambda_RefGenome.fa # create index
# 
# create fasta index for samtools
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome directory
samtools faidx Lambda_RefGenome.fa # create index
# 
# run bwa, pipe to samtools 
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
bwa mem -x ont2d /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda_RefGenome.fa Lambda.fastq | samtools view -T /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda_RefGenome.fa -bS - | samtools sort -T Lambda.bwa -o Lambda.bwa.bam -
# options: -x define read type

# part 2: LAST
# 
# genome indexing
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome # move to genome directory
lastdb Lambda Lambda_RefGenome.fa  # create index
# 
# convert fastq to fasta for last
# note: I tried it on fastq with -Q1 option (def fastq-sanger format) but got error 
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq # move to genome
seqtk seq -a Lambda.fastq > Lambda.fasta
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
python /home/aechchik/bin/nanopore-scripts/maf-convert.py sam Lambda.maf | samtools view -T /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/ref_genome/Lambda_RefGenome.fa -bS - | samtools sort -T Lambda.last -o Lambda.last.bam -

# part 3: get basic statistics on alignment 
#
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/fastq
# bam sort
samtools sort Lambda.bwa.bam > Lambda.bwa_sorted.bam
samtools sort Lambda.last.bam > Lambda.last_sorted.bam
# index on sorted bam 
samtools index Lambda.bwa_sorted.bam
samtools index Lambda.last_sorted.bam
# get stats of sorted bam
samtools stats Lambda.bwa_sorted.bam > Lambda.bwa.stats
samtools stats Lambda.last_sorted.bam > Lambda.last.stats
# get coverage from stats
grep ^COV Lambda.bwa.stats > Lambda.bwa.coverage
grep ^COV Lambda.last.stats > Lambda.last.coverage


