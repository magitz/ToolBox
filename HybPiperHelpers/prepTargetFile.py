#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import os
import re

# Prepare fasta files of individual genes to be used in a HybPiper target file
#   * Shortens name to first two underscore seperated items in record id (hopefully, genus_species)
#   * Adds a gene number as passed in on command line: >genus_speces-G###
#   * Removes reference sequence as indicated on command line
#   * Resulting files can be concatinated into a single file for the target file.
#
#
# Matt Gitzendanner
# University of Florida
# 
# version: 1.0: Initial release, Jan 7, 2019
#
# Usage: python prepTargetFile.py -i <fasta_file.fa> -n <gene name/number to use>\
#                        -p <prefix to add to gene name> -r <remove taxa with this text> -o <output_fileile.fa>
#
# To run in a loop, do something like:
#   count=1
#   for file in Folder/*
#   do 
#      echo Gene $count
#      python ../prepTargetFile.py -i ${file} -n $count -p G -r AT -o renamed/$file
#      count=$((count+1))
#   done
# 

# Deal with command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input fasta file")
parser.add_argument("-n", help="gene name/number")
parser.add_argument("-p", help="Add prefix to all genes", default=None)
parser.add_argument("-r", help="Remove taxa with this string in their name", default=None)
parser.add_argument("-o", help="Output fasta file with sequences renamed")

args = parser.parse_args()

InFile = args.i
OutFile= args.o
GeneName= args.n
GenePrefix=args.p
RemoveTaxa=args.r

# verify read and write input and output files
try:
	IN=open(InFile, 'r')
except:
	print("Can't open infile for reading: %s" %(InFile))
	exit()
try:
	OUT=open(OutFile, 'w')
except:
	print("Can't open outfile for writing: %s" %(OutFile))
	exit()	


# Parse input fasta file for sequence records.
for Record in SeqIO.parse(IN, "fasta"):	
	
	name_bits=Record.id.split("_")
	
    # Drop sequences specified by the -r option
	if RemoveTaxa in name_bits[0]:
		print("Removing " + name_bits[0])
		continue

	elif len(name_bits) > 1:
		Record.id= name_bits[0] + "_" + name_bits[1] + "-" + GenePrefix + GeneName
	else:
		Record.id= name_bits[0] + "_" + GeneName	

	Record.description="" # Clear anything in the description part of the fasta record.
		
    # Write the sequence to the output file.
	SeqIO.write(Record, OUT, "fasta")


