#!/bin/bash
#SBATCH --job-name=GATK_BWA_C300
#SBATCH -n 1
#SBATCH -p owners
#SBATCH --qos jlozier
#SBATCH --error=error.GATK_BWA_C300.%J.txt
#SBATCH --output=output.GATK_BWA_C300.%J.txt
#SBATCH --mem-per-cpu=160G
#SBATCH --mail-user=mlpimsler@ua.edu
#SBATCH -t 7-00:00:00

# Put in my user name
USER=`whoami`

# Plate
PLATE=C300

# Create directory on /stutmp
# -p option means to only make it if it doesn't exist already
mkdir -p /grps1/jlozier/$USER/VariantCalls/GATK/BWA/$PLATE
mkdir -p /grps1/jlozier/$USER/VariantCalls/GATK/BWA/$PLATE/gvcf

# Go to output directory
cd /grps1/jlozier/$USER/VariantCalls/GATK/BWA

# load dokit
export DK_ROOT=/share/apps/dotkit
. /share/apps/dotkit/bash/.dk_init

source /home/mlpimsler/.bash_profile

# Load Picard module
# Load GATK
use bioinfoJava
# Load samtools module
# Load BEDtools module
# Load BCFtools module
use bioinfoGCC
# Need to use a specific JAVA in order to use GATK
use java1.8.0

# Create variable that points to the reference genome
REFERENCE_DIR=/grps1/jlozier/mlpimsler/B_imp_bwa

# Create a directory to output the .gvcf files
mkdir -p /grps1/jlozier/$USER/VariantCalls/${PLATE}

cd /grps1/jlozier/$USER/VariantCalls/

# Cycle through each bifarius file and do variant calls for each 
for FILE in /grps1/jlozier/$USER/Alignments/BWA/RealignBAM/$PLATE/realigned-${FILENAME}.bam
        do
        # From the input file name, cuts the .fq for the creation of new files
        echo '---------------------------------'
        FILENAME=$(echo ${FILE##*/} | rev | cut -c 5- | rev | cut -c 11-)
        echo $FILENAME
        echo
        # Align the reads from one library to the genome
        # -M: Mark multiply aligned files
        # -t 4: use four threads

        # HaploType caller- SNP variant discovery- on un-realigned files
        # -T <string> : the GATK module to use
        # -R: reference genome
        # -I: input file (realigned .bam file)
        # -mmq <INT>: map quality cutoff 
        # --genotyping_mode {STRING}: 
        #       DISCOVERY- the program will choose the most likely alleles out of those it sees in the data
        #       GENOTYPE_GIVEN_ALLELES- only use the alleles passed in from a VCF file
        # -stand_emit_conf <INT>: minimum confidence threshold (phred-scaled) at which the program should emit 
        #       sites that appear to be possibly variant
        # -stand_call_conf <INT>: minimum confidence threshold (phred-scaled) at which the program should emit 
        #       variant sites as called. If a site's associated genotype has a confidence score lower than the 
        #       calling threshold, the program will emit the site as filtered and will annotate it as LowQual. 
        #       This threshold separates high confidence calls from low confidence calls.
        # -drf {STRING}: flag to igonore; using 'DuplicateRead' here because we expect RADtags to have duplicates
        # -ERC,--emitRefConfidence {STRING}: Mode for emitting reference confidence scores (NONE|BP_RESOLUTION|GVCF)
        # -o {STRING}: includes information about the mapping quality cutoff
        java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
                -T HaplotypeCaller \
                -R $REFERENCE_DIR/B_imp.fasta \
                -I $FILE \
                -mmq 30 \
                --genotyping_mode DISCOVERY \
                -stand_emit_conf 10 \
                -stand_call_conf 10 \
                -drf DuplicateRead \
                -ERC GVCF \
                -o ${FILENAME}_30.g.vcf

        echo 
done

echo "************************************"
pwd

# Variant calling for bifarius files only
echo "ls spp_coverage.g.vcf> spp-coverage_gvcf.list"

# GATK variant calling 
# First need to create a list of input files as a .list to pass to GATK
\ls /grps1/jlozier/$USER/VariantCalls/${PLATE}/*bif_30.g.vcf > /grps1/jlozier/$USER/VariantCalls/${PLATE}_bif_30_gvcf.list 
\ls /grps1/jlozier/$USER/VariantCalls/${PLATE}/*vos_30.g.vcf > /grps1/jlozier/$USER/VariantCalls/${PLATE}_vos_30_gvcf.list 

## Change directory
cd /grps1/jlozier/$USER/VariantCalls

# bifarius files
# quality 10
# Combine g.vcf files  
java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -R $REFERENCE_DIR/B_imp.fasta  \
        -drf DuplicateRead \
        -V ${PLATE}_bif_30_gvcf.list \
        -o ${PLATE}_bif_30_combined.g.vcf

# Vos files files
# quality 10
# Combine g.vcf files  
java -jar -Xmx160G $JARFILES/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -R $REFERENCE_DIR/B_imp.fasta  \
        -drf DuplicateRead \
        -V ${PLATE}_vos_30_gvcf.list \
        -o ${PLATE}_vos_30_combined.g.vcf

