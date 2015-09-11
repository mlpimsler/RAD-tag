#!/usr/bin/python

# Created by Meaghan L Pimsler
# 10 September 2015
''' This script needs a .csv file with two columns (barcodes in first, sample ID in second) and all of the files created from the process_radtags command and changes the names of all of the files in the directory to the names provided in the 'sample id' column.
    If it does not work because of the .csv file, it might mean that you saved it in the wrong format- if you are working on a mac you should have saved the file as windows version .csv'''

# Open relevant modules
import csv
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='This is a script that renames files created by process-radtags. Must have two column .csv in Windows format with barcodes in first column and sample ID name in second column. Usage: Rename_split_files.py -i INPUT.csv')
parser.add_argument('-i','--input', help='Input file name',required=True)
args = parser.parse_args()

# Code to read in the appropriate names and change the file names
#
with open(args.input) as csvfile:
    names = csv.reader(csvfile, dialect='excel')
    for row in names:
        if row[0] != 'Barcode':
            short_barcode = str(row[0])[:10]
            old_name = str('sample_' + short_barcode + '.fq')
            new_name = str(row[1] + '.fq')
            print str(row[0]) + ', ' + short_barcode + ', ' + old_name + ', ' + new_name
            os.rename(old_name,new_name)
        else:
            print row[0] + ', ' + '10bp barcode' + ', ' + 'old name' + ', ' + 'new name'

print "The files are: %s"%os.listdir(os.getcwd())
