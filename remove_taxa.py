#!/usr/bin/env python

# remove_taxa.py
# Takes a comma separated list of taxon names to remove from a fasta alignment.
# Sequences with matching names are not included in the output file.
#
# Matt Gitzendanner
# University of Florida

import argparse
from Bio import SeqIO
import os, sys


#Parse commandline options.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input fasta file")
parser.add_argument("-r", help="Comma separated list of taxa to remove")
parser.add_argument("-o", help="Output file")

args = parser.parse_args()
in_file = args.i
taxon_args = args.r
out_file = args.o

#Open the infile
try:
	IN = open (in_file, 'r')
except OSError as exception:
	print("Can't open input file:", in_file)
    sys.exit(1)

#Open the infile
try:
	OUT = open (out_file, 'a')
except OSError as exception:
	print("Can't open output file:", out_file)
    sys.exit(1)
		
#Convert taxon list to python list
taxon_list=  taxon_args.split(",")
	
for record in SeqIO.parse(IN, "fasta") : #Read each record from the in_file
   if record.id in taxon_list:
      continue # Don't include taxa in the taxon_list
   else:
      SeqIO.write(record, OUT, "fasta")
	
OUT.close
