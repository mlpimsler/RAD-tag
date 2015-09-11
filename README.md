# RAD-tag
This will be a workflow for RAD-Tag data
It will include scripts and programs for the bioinformatics of RAD-tag data analysis.
Hopefully this will include generalized programs.

# Rename_split_files.py
This script needs a .csv file with two columns (barcodes in first, sample ID in second) and all of the files created from the 'process_radtags' command and changes the names of all of the files in the directory to the names provided in the 'sample id' column.
If it does not work because of the .csv file, it might mean that you saved it in the wrong format- if you are working on a Mac you should have saved the file as windows version .csv
