#!/bin/sh
# Script to create an index, map reads, and then return results

# Put in my user name
USER=`whoami`

# Create directory on /scratch
mkdir -p  /scratch/$USER

# Copy input files to scratch
cp B.impatiens_genome/B_imp-ncbi.fa /scratch/$USER/.
mv /scratch/$USER/B_imp-ncbi.fa /scratch/$USER/B_imp.fasta
cp -rf TestRename/files_only /scratch/$USER/.

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

for DIR in ./*
do
        if [ -d "$DIR" ]
                then
                echo $(basename "$DIR")
                cd $DIR
                for FILE in *.fq
                        do
                        # From the input file name, cuts the .fq for the creation of new files
                        FILENAME = echo ${FILE##*/} | rev | cut -c 4- | rev
                        
                        READGROUP = @RG\tID:$DIR\tSM:$FILENAME\tPL:illumina\tLB:$FILENAME 
                        
                        # Align the reads from one library to the genome
                        # -M: Mark multiply aligned files
                        # -t 8: use four threads
                
                        bwa mem -M -R '<$READGROUP>' -t 8 -p B_imp.fasta $FILE > aligned-$FILENAME.sam

                        # Convert files from .sam to .bam, sorted
                        java -jar SortSam.jar \
                                INPUT=aligned-$FILENAME.sam \
                                OUTPUT=sorted-$FILENAME.bam \
                                SORT_ORDER=coordinate 

                        # Mark duplicates
                        java -jar MarkDuplicates.jar \
                                INPUT=sorted-$FILENAME.bam \
                                OUTPUT=dedup-$FILENAME.bam \
                                METRICS_FILE=metrics.txt

                        # Create index file from .bam
                        java -jar BuildBamIndex.jar \
                                INPUT=dedup_read.bam

                        # prints index statistics
                        samtools idxstats test_-$FILENAME.sorted.bam > test_-$FILENAME.sorted.stats
                done
        fi
done

# Creates a directory on my home directory for the results IF that directory doesn’t exist
# -p option is what gives it this capability
mkdir -p /home/$USER/Test
# Copies the .bam files to that directory
cp ./* /home/$USER/Test/.