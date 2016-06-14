# MinION practicals

### get the data

# fast5 location 
ls /home/jroux/archive/MinION/ #run_MinION_2015_06_15_Burn-In.tar
# extract to my /scratch
tar -xf /home/jroux/archive/MinION/run_MinION_2015_06_15_Burn-In.tar -C /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/
# touch files to keep them alive
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/
find . -exec touch {} \;

# load needed modules for the session 
module add UHTS/Analysis/poretools/0.5.1 # load poretools (read extraction)
module add UHTS/Aligner/bwa/0.7.13 # load bwa (alignment)
module add UHTS/Analysis/samtools/1.3 # load samtools (alignment)
module add SequenceAnalysis/SequenceAlignment/last/531 # load last (alignment)
module add UHTS/Analysis/seqtk/2015.10.15 # fasatq to fasta
module add UHTS/Quality_control/quast/4.1 # post-alignment control 
# other software in my /home:
ls /home/aechchik/bin
# TODO: ask Sébastien to install them on Vital-IT?


### read extraction

cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/wdir

# convert 2D reads to fast5
poretools fastq --type 2D *.fast5 > Lambda2D.fastq
# 41M     Lambda2D.fastq

# check read length
less Lambda2D.fastq | grep -E '^[ACTGN]+$' | while read rawseq; do echo -n "$rawseq" | wc -c ; done > Lambda2D_readlength.txt

# some stats
less Lambda2D.fastq | grep ^@ | grep .fast5$ | wc -l # how many reads in the input
awk '{ total += $1; count++ } END { print total/count }' readlength.txt # mean read length
sort -nk 1 Lambda2D_readlength.txt | head -n 1 # min read length
sort -nrk 1 Lambda2D_readlength.txt | head -n 1 # max read length
poretools stats ../Downloads_burnin/*.fast5 > Lambda_stats.txt # OR, get infos with poretools stats on fast5 directly


### read alignment to reference

# TODO: bowtie or last?

# part 0: get ref genome 
# 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna
mv NC_001416.fna Lambda_RefGenome.fa

# part 1: BWA
#
bwa index Lambda_RefGenome.fa # # genome indexing for bwa
samtools faidx Lambda_RefGenome.fa # create fasta index for samtools
bwa mem -x ont2d Lambda_RefGenome.fa Lambda2D.fastq | samtools view -T Lambda_RefGenome.fa -bS - | samtools sort -T Lambda.bwa -o Lambda.bwa.bam - # run bwa, pipe to samtools 
# options: -x define read type

# part 2: LAST
# 
lastdb Lambda Lambda_RefGenome.fa  # genome indexing
seqtk seq -a Lambda2D.fastq > Lambda2D.fasta # convert fastq to fasta for last
lastal -q 1 -a 1 -b 1 Lambda Lambda2D.fasta > Lambda.maf # run lastal on fasta reads
# options: -q = mismatch cost; -a = gap cost; -b = gap extension cost; -Q1 = read input in sanger-fastq format 
python /home/aechchik/bin/nanopore-scripts/maf-convert.py sam Lambda.maf | samtools view -T Lambda_RefGenome.fa -bS - | samtools sort -T Lambda.last -o Lambda.last.bam - # convert maf to bam with complete CIGAR (matches and mismatches)

# part 3: get basic statistics on alignment 
#
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


### assembly

# overlap
/home/aechchik/bin/minimap/minimap -Sw5 -L100 -m0 -t8 Lambda2D.fastq Lambda2D.fastq | gzip -1 > Lambda2D_reads.paf.gz
# layout
/home/aechchik/bin/miniasm/miniasm -f Lambda2D.fastq Lambda2D_reads.paf.gz > Lambda2D_contigs.gfa
awk '/^S/{print ">"$2"\n"$3}' Lambda2D_contigs.gfa | fold > Lambda2D_contigs.fa # convert contigs to fasta 

# check assembly quality 
quast.py -R  Lambda_RefGenome.fa Lambda2D_contigs.fa # scp pdf report?

# align reads to contigs?
bwa index Lambda2D_contigs.fa
bwa bwasw -t 8 Lambda2D_contigs.fa Lambda2D.fastq > Lambda2D_readmap.sai
bwa bwasw Lambda2D_contigs.fa Lambda2D_readmap.sai Lambda2D.fastq > Lambda2D_readmap.sam
