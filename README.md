# RAD-tag
This will be a workflow for RAD-Tag data
It will include scripts and programs for the bioinformatics of RAD-tag data analysis.
Hopefully this will include generalized programs.

# Rename_split_files.py
This script needs a .csv file with two columns (barcodes in first, sample ID in second) and all of the files created from the 'process_radtags' command and changes the names of all of the files in the directory to the names provided in the 'sample id' column.
If it does not work because of the .csv file, it might mean that you saved it in the wrong format- if you are working on a Mac you should have saved the file as windows version .csv

# sbatch-BWA-plate_by_plate.sh
This script is set up for submission to the University of Alabama High Performance Computing Cluster. It assumes that your RADtag fastq files have already been split, trimmed, and filtered with Stacks and are in plate-specific directories. These plate-specific directories are important for the read group information that GATK and samtools (can) make use of in the variant calling steps. The script will first cycle through the directories, then through each fastq file, and do the following steps:
> bwa mem alignment

> samstat alignment quality summary and assessment (samstat avaialble at: http://samstat.sourceforge.net/)

> pull out the unaligned reads and compress them for storage purposes (for possible downstream applications)

> sort and index the alignments

> Run GATK indel realignment

# sbatch-GSNAP-plate_by_plate.sh
This script is set up for submission to the University of Alabama High Performance Computing Cluster. It assumes that your RADtag fastq files have already been split, trimmed, and filtered with Stacks and are in plate-specific directories. These plate-specific directories are important for the read group information that GATK and samtools (can) make use of in the variant calling steps. The script will first cycle through the directories, then through each fastq file, and do the following steps:
> gsnap alignment using the --split-output and --failed-input options

>     --split-output will create multiple flies; nomapping, halfmapping_uniq, halfmapping_mult, unpaired_uniq, unpaired_mult, paired_uniq, paired_mult, concordant_uniq, and concordant_mult results

>     -failed-input creates a .fastq or .fasta file of all of the unmapped reads

> samstat alignment quality summary and assessment (samstat avaialble at: http://samstat.sourceforge.net/)

> pull out the unaligned reads and compress them for storage purposes (for possible downstream applications)

> sort and index the alignments

> Run GATK indel realignment
