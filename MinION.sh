# MinION 

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
module add UHTS/Analysis/poretools/0.5.1; poretools fastq input.fast5 > output.fastq

# -
# Q: I want one fastq per fast5 input
# A: Use bash loop instead
#
# cd /your_dir/MinION/basecalled_reads
a=0; for i in $(ls *.fast5); do echo $i;a=$(echo $i | cut -d'.' -f1); echo $a; b=$((b + 1)); module add UHTS/Analysis/poretools/0.5.1; poretools fastq $i > $a'.fastq'; done
# organize subfolders
mkdir ./fast5; mkdir ./fastq
mv *.fast5 ./fast5; mv *.fastq ./fastq

# -
# Q: How many resources does it take to convert fast5 to fastq?
# A: For an archive of 256M (fast5), CPU time: 12.19 sec; Max Memory: 29.20 MB

# -
# Q: How to check the length of each read sequenced by MinION?
# A: extract raw sequence from fastq, then count the number of characters per line
#
# cd /your_dir/MinION/basecalled_reads
less output.fastq | grep -E '^[ACTGN]+$' | while read rawseq; do echo -n "$rawseq" | wc -c ; done > readlength.txt

# -
# Q: How to check the average length of the reads sequenced by MinION?
# A: use extracted raw sequences character count, then calculate the mean
#
# cd /your_dir/MinION/basecalled_reads
awk '{ total += $1; count++ } END { print total/count }' reads_length.txt > mean_readlength.txt
