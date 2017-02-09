#!/bin/bash
#SBATCH --job-name=BWA_alignmment
#SBATCH -n 1
#SBATCH -p owners
#SBATCH --qos jlozier
#SBATCH --error=error.BWA_alignmment.%J.txt
#SBATCH --output=output.BWA_alignmment.%J.txt
#SBATCH --mem-per-cpu=160G
#SBATCH --mail-user=mlpimsler@ua.edu
#SBATCH -t 7-00:00:00

# Creates a variable for your user name
USER=`whoami`

# Create directory for output files
# -p option means to only make it if it doesn't exist already
mkdir -p /grps1/jlozier/$USER/Alignments/BWA
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/SAM
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/BAM
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/BamStat
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/SamStat
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/UnMapped
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/RealignBAM
mkdir -p /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/Intervals

# Go to output directory
cd /grps1/jlozier/$USER/Alignments/BWA

# load dokit
export DK_ROOT=/share/apps/dotkit
. /share/apps/dotkit/bash/.dk_init

# Load your bash preferences
# Especially useful if you have locally installed files and/or paths
source /home/$USER/.bash_profile

# Load Picard module
# Load GATK
use bioinfoJava
# load BWA
# Load BEDtools module
# Load BCFtools module
use bioinfoGCC
# Need to use a specific JAVA in order to use GATK
use java1.8.0

# Create variable that points to the reference genome
REFERENCE_SEQ=/grps1/jlozier/mlpimsler/B_imp_bwa/B_imp.fasta

for DIR in /grps1/jlozier/mlpimsler/SplitRenamed/C*
do
        # if the * wildcard is actually a directory
        if [ -d "$DIR" ]
                then
                # Print the directory
                echo $(basename "$DIR")
                # go "into" the directory
                cd $DIR
                # Cycle through each file in the directory that has a .fq ending
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/SAM/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/BAM/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/BamStat/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/SamStat/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/UnMapped/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/Intervals/$(basename "$DIR")
                for FILE in *.fq
                        do
                        # From the input file name, cuts the .fq for the creation of new files that retain sample information
                        FILENAME=$(echo ${FILE##*/} | rev | cut -c 4- | rev)
                        
                        # Creates a variable for that directory for the readgroup information
                        RUN=$(basename "$DIR")
                        echo '---------------------------------'
                        echo $FILENAME
                        echo
                        
                        cd /grps1/jlozier/$USER/Alignments/BWA

                        ### Alignment using BWA
                        # Align the reads from one library to the genome
                        # -M: Mark multiply aligned files
                        # -t 4: use four threads
                        bwa mem -R "@RG\tID:${RUN}\tSM:${FILENAME}\tPL:illumina\tLB:${FILENAME}" \
                                -t 1 -p $REFERENCE_SEQ $FILE > aligned-$FILENAME.sam
                        
                        # Calculate summary statistics about the alignments
                        samstat aligned-$FILENAME.sam  
                                                                                                
                        # Convert files from .sam to .bam, sorted
                        samtools view -Shu aligned-$FILENAME.sam | samtools sort - sorted-$FILENAME
                        
                        # Pull out the unmapped reads
                        samtools view -f4 aligned-$FILENAME.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > unmapped-$FILENAME.fastq
                        # Compresses the unmapped files
                        gzip unmapped-$FILENAME.fastq
                        
                        # Indexes the alignment results
                        samtools index sorted-$FILENAME.bam
                        # Saves index statistics
                        samtools idxstats sorted-$FILENAME.bam > sorted-$FILENAME.bam.stats

                        ### Realignment using GATK
                        #Identify regions in need of realignment:
                        #  defaults for optional parameters:
                        #  --minReadsAtLocus N [the minimum coverage at a locus for the entropy calculation to be enabled; default=4]
                        #  --windowSize N [any two SNP calls and/or high entropy positions are considered clustered when they occur no 
                        #       more than N basepairs apart; default=10]
                        #  --mismatchFraction f [fraction of total sum of base qualities at a position that need to mismatch for the 
                        #       position to be considered to have high entropy; default=0.15; to disable, set to <= 0 or > 1]
                        #       Note that this fraction should be adjusted based on your particular data set. For deep coverage and/or 
                        #       when looking for indels with low allele frequency, this number should be smaller.
                        #  --maxIntervalSize [max size in bp of intervals that we'll pass to the realigner; default=500]
                        java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
                                -T RealignerTargetCreator \
                                -R $REFERENCE_SEQ \
                                -o /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/Intervals/$(basename "$DIR")/${FILENAME}.intervals \
                                -I sorted-${FILENAME}.bam \
                                --minReadsAtLocus 3

                        #  Run realigner over intervals:
                        # Optional parameters:
                        # -compress 0: this argument recommended to speed up the process *if* this is only a temporary file; 
                        #       otherwise, use the default value defaults for optional parameters:
                        # -compress, --bam_compression; Compression level to use for output bams; [default:5].
                        # -LOD, --LODThresholdForCleaning; LOD threshold above which the realigner will proceed to realign; default=5.0]
                        #       This term is equivalent to "significance" - i.e. is the improvement significant enough to merit realignment? 
                        #       Note that this number should be adjusted based on your particular data set. For low coverage and/or when 
                        #       looking for indels with low allele frequency, this number should be smaller.
                        # -targetNotSorted, --targetIntervalsAreNotSorted; This tool assumes that the target interval list is sorted; 
                        #       if the list turns out to be unsorted, it will throw an exception. Use this argument when your interval 
                        #       list is not sorted to instruct the Realigner to first sort it in memory.
                        # -knownsOnly, --useOnlyKnownIndels; Don't run 'Smith-Waterman' to generate alternate consenses; 
                        #       use only known indels provided as RODs for constructing the alternate references. 

                        java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
                                -I sorted-${FILENAME}.bam \
                                -R $REFERENCE_SEQ \
                                -T IndelRealigner \
                                -targetIntervals /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/Intervals/$(basename "$DIR")/${FILENAME}.intervals \
                                -o /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/$(basename "$DIR")/realigned-${FILENAME}.bam \
                                -LOD 3.0 \
                                --maxReadsInMemory 1000000 \
                                --maxReadsForRealignment 100000
 
                done
                # End of the loop through the files in the directory
                
                ## Move files to their appropriate destinations
                mv *.sam /grps1/jlozier/$USER/Alignments/BWA/SAM/$(basename "$DIR")/.
                mv *.html /grps1/jlozier/$USER/Alignments/BWA/SamStat/$(basename "$DIR")/.
                mv *.bam /grps1/jlozier/$USER/Alignments/BWA/BAM/$(basename "$DIR")/.
                mv *.bai /grps1/jlozier/$USER/Alignments/BWA/BAM/$(basename "$DIR")/.
                mv *.stats /grps1/jlozier/$USER/Alignments/BWA/BamStat/$(basename "$DIR")/.
                mv unmapped* /grps1/jlozier/$USER/Alignments/BWA/UnMapped/$(basename "$DIR")/.
                cd ..
                echo "************************************"

        else
                echo
                echo $DIR
                pwd
                echo
        fi
        
done
