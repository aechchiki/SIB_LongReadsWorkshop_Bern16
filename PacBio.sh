# PacBio practicals

### get the data 

# reads location
/scratch/beegfs/monthly/SIB_long_read_workshop/pacbio_lambda_reads/
# link them to my scratch
cd /scratch/beegfs/monthly/<mydir>
mkdir pb_reads
ln -s /scratch/beegfs/monthly/SIB_long_read_workshop/pacbio_lambda_reads/* ./pb_reads/

# load module for read extraction
module add UHTS/PacBio/pbh5tools/0.8.0;

# get consensus reads (reads_of_insert.fastq)
cd /scratch/beegfs/monthly/aechchik/SIB_Bern16/minion/wdir

bsub -q priority bash5tools.py <file.bas.h5> --minReadScore 0.85 --minLength 2500 --outFilePrefix pb_reads_fitered --readType subreads --outType fastq

cat reads_of_insert.fastq |  grep @ | wc -l # how many reads in the input
# note: 3068 in minION 
# select same number from PacBio reads:

head -12272 pacbio_reads.fastq > sel_PacBio.fastq



# some stats
less sel_PacBio.fastq | grep -E '^[ACTGN]+$' | while read rawseq; do echo -n "$rawseq" | wc -c ; done > sel_PacBio_readlength.txt
awk '{ total += $1; count++ } END { print total/count }' sel_PacBio_readlength.txt # mean read length
sort -nk 1 sel_PacBio_readlength.txt | head -n 1 # min read length
sort -nrk 1 sel_PacBio_readlength.txt | head -n 1 # max read length

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
bwa mem -x pacbio Lambda_RefGenome.fa sel_PacBio.fastq | samtools view -T Lambda_RefGenome.fa -bS - | samtools sort -T Lambda.bwa -o Lambda.bwa.bam - # run bwa, pipe to samtools 
# options: -x define read type

# part 2: LAST
# 
lastdb Lambda Lambda_RefGenome.fa  # genome indexing
seqtk seq -a sel_PacBio.fastq > sel_PacBio.fasta # convert fastq to fasta for last
lastal -q 1 -a 1 -b 1 Lambda sel_PacBio.fasta > Lambda.maf # run lastal on fasta reads
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

module add UHTS/Analysis/minimap/0.2.r124.dirty;
module add UHTS/Analysis/miniasm/0.2.r137.dirty;

#minimap -Sw5 -L100 -m0 -t8 RSII_selection.fq RSII_selection.fq | gzip -1 > Lambda_reads.paf.gz
#miniasm -f RSII_selection.fq Lambda_reads.paf.gz > Lambda_contigs.gfa
#awk '/^S/{print ">"$2"\n"$3}' Lambda_contigs.gfa | fold > Lambda_contigs.fa

# overlap
/home/aechchik/bin/minimap/minimap -Sw5 -L100 -m0 -t8 sel_PacBio.fastq sel_PacBio.fastq | gzip -1 > Lambda_reads.paf.gz
# layout
/home/aechchik/bin/miniasm/miniasm -f sel_PacBio.fastq Lambda_reads.paf.gz > Lambda_contigs.gfa
awk '/^S/{print ">"$2"\n"$3}' Lambda_contigs.gfa | fold > Lambda_contigs.fa # convert contigs to fasta 

# check assembly quality 
quast.py -R  Lambda_RefGenome.fa Lambda_contigs.fa # scp pdf report?

