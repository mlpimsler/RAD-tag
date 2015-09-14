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
# Load BEDtools module
module load bedtools/2.17.0
# Load GATK module
module load gatk/3.4-46

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
                        FILENAME=$(echo ${FILE##*/} | rev | cut -c 4- | rev)
                        
                        # Align the reads from one library to the genome
                        # -M: Mark multiply aligned files
                        # -t 8: use eight threads
                
                        bwa mem -M -R '@RG\tID:$DIR\tSM:$FILENAME\tPL:illumina\tLB:$FILENAME' -t 8 -p B_imp.fasta $FILE > aligned-$FILENAME.sam

                        # First creates a .bam file and then pipes it to the sort tool
                        # -bu: -b creates .bam file, -u indicates “don’t compress”- only used because I am $
                        # new tool
                        samtools view -bu aligned-$FILENAME.sam -o aligned-$FILENAME.bam

                        #contrary to the manual, fomat for samtools is: samtools sort [input] [output prefix]
                        samtools sort aligned-$FILENAME.bam sorted_$FILENAME
                        # Indexes the results
                        samtools index sorted_$FILENAME.bam
                        # prints index statistics
                        samtools idxstats sorted_$FILENAME.bam > $FILENAME.bam.stats
#                       
                        # Mark duplicates
#                        java -jar MarkDuplicates.jar \
#                                INPUT=sorted-$FILENAME.bam \
#                                OUTPUT=dedup-$FILENAME.bam \
#                                METRICS_FILE=metrics.txt


                        # Get information regarding the coverage
                        genomeCoverageBed -ibam sorted_$FILENAME.bam -g B_imp.fasta > $FILENAME-BEDcoverage.txt

                        #Identify regions in need of realignment:
                        java -Xmx2g -jar GenomeAnalysisTK.jar \
                                -T RealignerTargetCreator \
                                -R REFERENCE_ASSEMBLY_NAME.fasta \
                                -o merged_output.intervals \
                                -I merged.bam \
                                --minReadsAtLocus 3
                        #  defaults for optional parameters:
                        #  --minReadsAtLocus N [the minimum coverage at a locus for the entropy calculation to be enabled; default=4]
                        #  --windowSize N [any two SNP calls and/or high entropy positions are considered clustered when they occur no more than N basepairs apart; default=10]
                        #  --mismatchFraction f [fraction of total sum of base qualities at a position that need to mismatch for the position to be considered to have high entropy; default=0.15; to disable, set to <= 0 or > 1]
                        #  Note that this fraction should be adjusted based on your particular data set. For deep coverage and/or when looking for indels with low allele frequency, this number should be smaller.
                        #  --maxIntervalSize [max size in bp of intervals that we'll pass to the realigner; default=500]

                        #  Run realigner over intervals:
                        java -Xmx4g -jar ~/programs/GATK/dist/GenomeAnalysisTK.jar \
                                -I merged.bam \
                                -R REFERENCE_ASSEMBLY_NAME.fasta \
                                -T IndelRealigner \
                                -targetIntervals merged_output.intervals \
                                -o merged_realigned.bam \
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
                        java -jar GenomeAnalysisTK.jar \
                                -T HaplotypeCaller \
                                -R B_imp.fasta \
                                -I sorted-$FILENAME.bam \
                                --genotyping_mode DISCOVERY \
                                -stand_emit_conf 10 \
                                -stand_call_conf 30 \
                                -emitRefConfidence GVCF
                                -o $FILENAME-raw_variants.g.vcf
                done
                

                cd ..

        fi
done

# Creates a directory on my home directory for the results IF that directory doesn’t exist
# -p option is what gives it this capability
mkdir -p /home/$USER/Test
# Copies the .bam files to that directory
cp ./* /home/$USER/Test/.