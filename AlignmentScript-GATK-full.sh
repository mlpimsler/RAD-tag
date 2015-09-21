#!/bin/sh
# Script to create an index, map reads, and then return results

# Put in my user name
USER=`whoami`

# Create directory on /scratch
# -p option means to only make it if it doesn't exist already
mkdir -p /scratch/$USER
# Make sure that the directory is empty
rm -rf /scratch/$USER/*

# Copy input files to scratch
cp B.impatiens_genome/B_imp-ncbi.fa /scratch/$USER/.
mv /scratch/$USER/B_imp-ncbi.fa /scratch/$USER/B_imp.fasta
cp -rf C3* /scratch/$USER/.

# Go to scratch
cd /scratch/$USER
ls -alF

# Load BWA module
source /opt/asn/etc/asn-bash-profiles/modules.sh
module load bwa/0.7.12
# Load samtools module
module load samtools/1.2
# Load Picard module
module load picard/1.79
# Load BEDtools module
module load bedtools/2.17.0
# Load GATK module
module load gatk/3.4-46

#Create index of reference genome
bwa index -a bwtsw B_imp.fasta
# Prepare a .dict index (necessary for GATK)
java -Xmx16g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar \
        R= B_imp.fasta \
        O= B_imp.dict

# Prepare a .fai index (necessary for GATK)
samtools faidx B_imp.fasta

# Cycle through each directory
# Each directory should contain the pre-renamed split files, one lane per directory
# Use the renamed_split_files.py python script for this
for DIR in ./*
do
        # if the * wildcard is actually a directory
        if [ -d "$DIR" ]
                then
                # Print the directory
                echo $(basename "$DIR")
                # go "into" the directory
                cd $DIR
                # Cycle through each file in the directory that has a .fq ending
                for FILE in *.fq
                        do
                        # From the input file name, cuts the .fq for the creation of new files
                        FILENAME=$(echo ${FILE##*/} | rev | cut -c 4- | rev)
                        
                        echo '---------------------------------'
                        echo $FILENAME
                        echo
                        # Align the reads from one library to the genome
                        # -M: Mark multiply aligned files
                        # -t 4: use four threads
                
                        bwa mem -M -R '@RG\tID:$DIR\tSM:$FILENAME\tPL:illumina\tLB:$FILENAME' -t 4 -p ../B_imp.fasta $FILE > aligned-$FILENAME.sam

                        # Convert files from .sam to .bam, sorted
                        java -Xmx16g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar \
                                INPUT=aligned-$FILENAME.sam \
                                OUTPUT=sorted-$FILENAME.bam \
                                SORT_ORDER=coordinate 
                        
                        # Indexes the results
                        samtools index sorted-$FILENAME.bam
                        # prints index statistics
                        samtools idxstats sorted-$FILENAME.bam > $FILENAME.bam.stats
                        #                       
                        # Mark duplicates
#                        java -jar MarkDuplicates.jar \
#                                INPUT=sorted-$FILENAME.bam \
#                                OUTPUT=dedup-$FILENAME.bam \
#                                METRICS_FILE=metrics.txt

                        # Get information regarding the coverage
                        genomeCoverageBed -ibam sorted-$FILENAME.bam -g ../B_imp.fasta > $FILENAME-BEDcoverage.txt

                        #Identify regions in need of realignment:
                        java -Xmx16g -jar /mnt/homeapps/apps/dmc/apps/gatk_3.4-46/GenomeAnalysisTK.jar \
                                -T RealignerTargetCreator \
                                -R ../B_imp.fasta \
                                -o sorted-$FILENAME.intervals \
                                -I sorted-$FILENAME.bam \
                                --minReadsAtLocus 3
                        #  defaults for optional parameters:
                        #  --minReadsAtLocus N [the minimum coverage at a locus for the entropy calculation to be enabled; default=4]
                        #  --windowSize N [any two SNP calls and/or high entropy positions are considered clustered when they occur no more than N basepairs apart; default=10]
                        #  --mismatchFraction f [fraction of total sum of base qualities at a position that need to mismatch for the position to be considered to have high entropy; default=0.15; to disable, set to <= 0 or > 1]
                        #  Note that this fraction should be adjusted based on your particular data set. For deep coverage and/or when looking for indels with low allele frequency, this number should be smaller.
                        #  --maxIntervalSize [max size in bp of intervals that we'll pass to the realigner; default=500]

                        #  Run realigner over intervals:
                        java -Xmx16g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar \
                                -I sorted-$FILENAME.bam \
                                -R ../B_imp.fasta \
                                -T IndelRealigner \
                                -targetIntervals sorted-$FILENAME.intervals \
                                -o realigned-$FILENAME.bam \
                                -LOD 3.0 \
                                --maxReadsInMemory 1000000 \
                                --maxReadsForRealignment 100000
                        #Optional parameters:
                        # -compress 0 \
                        #    this argument recommended to speed up the process *if* this is only a temporary file; otherwise, use the default value
                        #    defaults for optional parameters:
                        # -compress, --bam_compression; Compression level to use for output bams; [default:5].
                        # -LOD, --LODThresholdForCleaning; LOD threshold above which the realigner will proceed to realign; default=5.0]
                        #    This term is equivalent to "significance" - i.e. is the improvement significant enough to merit realignment? Note that this number should be adjusted based on your particular data set. For low coverage and/or when looking for indels with low allele frequency, this number should be smaller.
                        # -targetNotSorted, --targetIntervalsAreNotSorted; This tool assumes that the target interval list is sorted; if the list turns out to be unsorted, it will throw an exception. Use this argument when your interval list is not sorted to instruct the Realigner to first sort it in memory.
                        # -knownsOnly, --useOnlyKnownIndels; Don't run 'Smith-Waterman' to generate alternate consenses; use only known indels provided as RODs for constructing the alternate references. 

                        # HaploType caller- SNP variant discovery
                        java -Xmx16g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar \
                                -T HaplotypeCaller \
                                -R ../B_imp.fasta \
                                -I realigned-$FILENAME.bam \
                                --genotyping_mode DISCOVERY \
                                -stand_emit_conf 10 \
                                -stand_call_conf 30 \
                                -drf DuplicateRead \
                                -ERC GVCF \
                                -o raw_variants-$FILENAME.g.vcf
                        echo 
                        echo
                done
                # End of the loop through the files in the directory

                echo "************************************"
                pwd
                echo "ls realigned-*.bam > bam_files.list"
                ## Samtools variant calling
                # First need to create a list of input files
                \ls realigned-*.bam > $(basename "$DIR")_bam_files.list
                more bam_files.list
                echo "samtools mpileup"
                # calling mpileup to make the .vcf file
                samtools mpileup -v -q 10 -f ../B_imp.fasta -t DP --output samtools_$(basename "$DIR").vcf -b $(basename "$DIR")_bam_files.list

                # GATK variant calling 
                # First need to create a list of input files as a .list to pass to GATK
                \ls raw_variants-*.g.vcf > $(basename "$DIR")_gvcf.list      
                #Combine g.vcf files          
                java -Xmx16g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar \
                        -T CombineGVCFs \
                        -R ../B_imp.fasta \
                        -drf DuplicateRead \
                        -V $(basename "$DIR")_gvcf.list \
                        -o $(basename "$DIR")_combined.g.vcf
                
                # genotyping/SNP calling
                java -Xmx16g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar \
                        -T GenotypeGVCFs \
                        -R ../B_imp.fasta \
                        -V $(basename "$DIR")_combined.g.vcf \
                        -o $(basename "$DIR")_combined.vcf
                cd ..
                pwd
        else
                echo
                echo $DIR
                pwd
                echo
        fi
        
done

for f in */realigned-*.bam; do [[ -d "$f" ]] || echo "$f"; done > all_bam_files.list
samtools mpileup -v -q 10 -f B_imp.fasta -t DP --output samtools_AllLibraries.vcf -b all_bam_files.list

for f in */*combined.g.vcf; do [[ -d "$f" ]] || echo "$f"; done > all_gvcf_files.list
java -Xmx16g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R B_imp.fasta \
        -V all_gvcf_files.list \
        -o all_combined.vcf

cd /scratch/$USER
# Creates a directory on my home directory for the results IF that directory doesn’t exist
# -p option is what gives it this capability
mkdir -p /home/$USER/Test/GATK_Align
# Copies the .bam files to that directory
cp -rf ./* /home/$USER/Test/GATK_Align/.
rm -rf /scratch/$USER/*