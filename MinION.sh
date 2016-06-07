# MinION 

# data location 
ls /archive/dee/robinson/jroux/MinION/Minion_Lambda_Library
# data already in .fastq


### read extraction

# Q: I have .fast5, how to convert to .fastq? 
# A: Use poretools on basecalled reads.

# to install poretools:
git clone https://github.com/arq5x/poretools
cd poretools
python setup.py install --root ./ # must have write rights, ideally /home/user/bin
# example
# cd /your_dir/MinION/basecalled
# /home/aechchik/bin/poretools/poretools fastq input.fast5 > output.fastq

# new! poretools is now on Vital-IT (thx Sébastien)
# example 
# cd /your_dir/MinION/basecalled
# module add UHTS/Analysis/poretools/0.5.1; poretools fastq input.fast5 > output.fastq

# note: 
# - poretools can also be used to convert to .fasta: write to output.fasta
# - can create one big .fastq for all the flowcell by defining *.fast5 as input (default: 1 .fast5 per pore)

# Q: How many resources does it take to convert fast5 to fastq?
# A: For an archive of 256M (fast5), CPU time: 2.19 sec; Max Memory: 29.20 MB

# Q: How to check the length of each read sequenced by MinION?
# A: extract raw sequence from fastq, then count the number of characters per line
less output.fastq | grep -E '^[ACTGN]+$' | while read rawseq; do echo -n "$rawseq" | wc -c ; done > reads_length.txt

