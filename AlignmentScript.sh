#!/bin/sh
# Script to create an index, map reads, and then return results

# Put in my user name
USER=`whoami`

# Create directory on /scratch
mkdir -p  /scratch/$USER

# Copy input files to scratch
cp B.impatiens_genome/B_imp-ncbi.fa /scratch/$USER/.
mv /scratch/$USER/B_imp-ncbi.fa /scratch/$USER/B_imp.fasta
cp TestRename/*.fq /scratch/$USER/.

# Go to scratch
cd /scratch/$USER

# Load BWA module
source /opt/asn/etc/asn-bash-profiles/modules.sh
module load bwa/0.7.12
# Load samtools module
module load samtools/1.2
# Load Picard module
module load picard/1.79

#Create index of reference genome
bwa index B_imp.fasta

for FILE in *.fq
        do
                # From the input file name, cuts the .fq for the creation of new files
                FILENAME = echo ${FILE##*/} | rev | cut -c 4- | rev
                
                # Align the reads from one library to the genome
                # -M: Mark multiply aligned files
                # -t 4: use four threads
                
                bwa mem -M -t 4 B_imp.fasta $FILE > test-$FILENAME.sam

                # First creates a .bam file and then pipes it to the sort tool
                # -bu: -b creates .bam file, -u indicates “don’t compress”- only used because I am $
                # new tool
                samtools view -bu test_-$FILENAME.sam -o test_-$FILENAME.bam

                #contrary to the manual, fomat for samtools is: samtools sort [input] [output prefi$
                samtools sort test_-$FILENAME.bam test_-$FILENAME.sorted
                # Indexes the results
                samtools index test_-$FILENAME.sorted.bam
                # prints index statistics
                samtools idxstats test_-$FILENAME.sorted.bam > test_-$FILENAME.sorted.stats
        done

# Creates a directory on my home directory for the results IF that directory doesn’t exist
# -p option is what gives it this capability
mkdir -p /home/$USER/Test
# Copies the .bam files to that directory
cp ./* /home/$USER/Test/.