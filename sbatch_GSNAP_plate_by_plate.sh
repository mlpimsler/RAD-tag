#!/bin/bash
#SBATCH --job-name=GSNAP_workflow
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p owners
#SBATCH --qos jlozier
#SBATCH --error=error.GSNAP_workflow.%J.txt
#SBATCH --output=output.GSNAP_workflow.%J.txt
#SBATCH --mem-per-cpu=160G
#SBATCH --mail-user=mlpimsler@ua.edu
#SBATCH -t 7-00:00:00

# Creates a variable for your user name
USER=`whoami`

# Create directory for output files
# -p option means to only make it if it doesn't exist already
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/SAM
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/BAM
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/BamStat
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/SamStat
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/RealignBAM
mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/RealignBAM/Intervals

# Go to output directory
cd /grps1/jlozier/mlpimsler/Alignments/GSNAP

# load dokit
export DK_ROOT=/share/apps/dotkit
. /share/apps/dotkit/bash/.dk_init

# Load your bash preferences
# Especially useful if you have locally installed files and/or paths
source /home/$USER/.bash_profile

# Load Picard module
# Load GATK
use bioinfoJava
# load GSNAP
# Load BEDtools module
# Load BCFtools module
use bioinfoGCC
# Need to use a specific JAVA in order to use GATK
use java1.8.0

# Create variable that points to the reference genome
REFERENCE_DIR=/grps1/jlozier/mlpimsler/B_imp_genome

for DIR in /grps1/jlozier/mlpimsler/SplitRenamed/C*
do
        # if the * wildcard is actually a directory
        if [ -d "$DIR" ]
                then
                # Print the directory
                echo $(basename "$DIR")
                # go "into" the directory
                cd $DIR

                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/SAM/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/SAM/$(basename "$DIR")/AdditionalFiles
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/BAM/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/BamStat/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/SamStat/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/RealignBAM/$(basename "$DIR")
                mkdir -p /grps1/jlozier/$USER/Alignments/GSNAP/RealignBAM/Intervals/$(basename "$DIR")
                # Cycle through each file in the directory that has a .fq ending
                for FILE in *.fq
                        do
                        # From the input file name, cuts the .fq for the creation of new files that retain sample information
                        FILENAME=$(echo ${FILE##*/} | rev | cut -c 4- | rev)
                        
                        # Creates a variable for that directory for the readgroup information
                        RUN=$(basename "$DIR")
                        echo '---------------------------------'
                        echo $FILENAME
                        echo

                        ### Alignment using GSNAP
                        # Align the reads from one library to the genome
                        # -t: number of threads
                        # -n: maximum number of paths to print (number of permitted alignment positions)
                        # -Q, --quiet-if-excessive: If more than maximum number of paths are found, then nothing is printed.
                        # --use-sarray=0: whether or not to use a suffix array; 1 is default, and includes a suffix array
                        #                       use of a suffix array will bias your results against 
                        # --split-output: creates multiple flies; nomapping, halfmapping_uniq, halfmapping_mult, unpaired_uniq, 
                        #               unpaired_mult, paired_uniq, paired_mult, concordant_uniq, and concordant_mult results
                        # --failed-input={string}: creates a .fastq or .fasta file of all of the unmapped reads
                        # --min-coverage={number}: minimum coverage of the read, percentage expressed betwee 0.0 and 1.0
                        # -m, --max-mismatches=FLOAT{int}: maximum number of mismatches
                        # -i, --indel-penalty=INT: Counts against mismatches allowed.  To find indels, make
                        #          indel-penalty less than or equal to max-mismatches. A value < 2 can lead to false positives at read ends
                        gsnap -t 8 -n 1 --quiet-if-excessive --use-sarray=0 \
                                --failed-input=/grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/unmapped-$FILENAME \
                                --split-output=/grps1/jlozier/$USER/Alignments/GSNAP/aligned-${FILENAME} \
                                --min-coverage=1.0 -A sam -m 5 -i 2 \
                                -D $REFERENCE_DIR -d B_imp \
                                --read-group-id="${RUN}" \
                                --read-group-name="${FILENAME}" \
                                --read-group-library="${FILENAME}" \
                                --read-group-platform="illumina" \
                                $FILE
                        
                        ## The output files don't include a .bam at the end, so some programs (samstat) don't recognize the file
                        # We need to rename the file we are interested in to make sure that it is a .sam file
                        mv /grps1/jlozier/$USER/Alignments/GSNAP/aligned-${FILENAME}.unpaired_uniq \
                                /grps1/jlozier/$USER/Alignments/GSNAP/aligned-${FILENAME}.unpaired_uniq.sam

                        ## Move the unneeded files out of the way
                        mv /grps1/jlozier/$USER/Alignments/GSNAP/aligned-${FILENAME}.unpaired_transloc \
                                /grps1/jlozier/$USER/Alignments/GSNAP/SAM/$(basename "$DIR")/AdditionalFiles/aligned-${FILENAME}.unpaired_transloc.sam
                        mv /grps1/jlozier/$USER/Alignments/GSNAP/aligned-${FILENAME}.unpaired_mult \
                                /grps1/jlozier/$USER/Alignments/GSNAP/SAM/$(basename "$DIR")/AdditionalFiles/aligned-${FILENAME}.unpaired_mult.sam
                        mv /grps1/jlozier/$USER/Alignments/GSNAP/aligned-${FILENAME}.nomapping \
                                /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/$(basename "$DIR")/aligned-${FILENAME}.nomapping.sam
                        mv /grps1/jlozier/$USER/Alignments/GSNAP/unmapped-${FILENAME}.1 \
                                /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/$(basename "$DIR")/unmapped-${FILENAME}.1.fastq
                        mv /grps1/jlozier/$USER/Alignments/GSNAP/unmapped-${FILENAME}.2 \
                                /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/$(basename "$DIR")/unmapped-${FILENAME}.2.fastq

                        # Compress the unmapped reads
                        gzip /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/$(basename "$DIR")/unmapped-${FILENAME}.1.fastq
                        gzip /grps1/jlozier/$USER/Alignments/GSNAP/UnMapped/$(basename "$DIR")/unmapped-${FILENAME}.2.fastq

                        cd /grps1/jlozier/$USER/Alignments/GSNAP

                        # Calculate summary statistics about the alignments
                        samstat aligned-$FILENAME.unpaired_uniq.sam  
                                                                                                
                        # Convert files from .sam to .bam, sorted
                        samtools view -Shu aligned-$FILENAME.unpaired_uniq.sam | samtools sort - sorted-$FILENAME
                        
                        # Compresses the unmapped files
                        gzip unmapped-$FILENAME.fastq
                        
                        # Indexes the alignment results
                        samtools index sorted-$FILENAME.bam
                        # Saves index statistics
                        samtools idxstats sorted-$FILENAME.bam > sorted-$FILENAME.bam.stats

                        ### Realignment using GATK
                        #Identify regions in need of realignment:
                        java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
                                -T RealignerTargetCreator \
                                -R ~/B_imp_bwa/B_imp.fasta \
                                -o ${FILENAME}.intervals \
                                -I $FILE \
                                --minReadsAtLocus 3
                        #  defaults for optional parameters:
                        #  --minReadsAtLocus N [the minimum coverage at a locus for the entropy calculation to be enabled; default=4]
                        #  --windowSize N [any two SNP calls and/or high entropy positions are considered clustered when they occur no more than N basepairs apart; default=10]
                        #  --mismatchFraction f [fraction of total sum of base qualities at a position that need to mismatch for the position to be considered to have high entropy; default=0.15; to disable, set to <= 0 or > 1]
                        #  Note that this fraction should be adjusted based on your particular data set. For deep coverage and/or when looking for indels with low allele frequency, this number should be smaller.
                        #  --maxIntervalSize [max size in bp of intervals that we'll pass to the realigner; default=500]

                        #  Run realigner over intervals:
                        java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
                                -I $FILE \
                                -R ~/B_imp_bwa/B_imp.fasta \
                                -T IndelRealigner \
                                -targetIntervals ${FILENAME}.intervals \
                                -o realigned-${FILENAME}.bam \
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

                done
                # End of the loop through the files in the directory
                
                ## Move files to their appropriate destinations
                mv *.sam /grps1/jlozier/$USER/Alignments/GSNAP/SAM/$(basename "$DIR")/.
                mv *.html /grps1/jlozier/$USER/Alignments/GSNAP/SamStat/$(basename "$DIR")/.
                mv realigned* /grps1/jlozier/$USER/Alignments/GSNAP/RealignBAM/$(basename "$DIR")/.
                mv *.intervals /grps1/jlozier/$USER/Alignments/GSNAP/RealignBAM/Intervals/$(basename "$DIR")/.
                mv *.bam /grps1/jlozier/$USER/Alignments/GSNAP/BAM/$(basename "$DIR")/.
                mv *.bai /grps1/jlozier/$USER/Alignments/GSNAP/BAM/$(basename "$DIR")/.
                mv *.stats /grps1/jlozier/$USER/Alignments/GSNAP/BamStat/$(basename "$DIR")/.

                cd ..
                echo "************************************"

        else
                echo
                echo $DIR
                pwd
                echo
        fi
        
done
