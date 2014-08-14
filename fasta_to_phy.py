#!/usr/bin/env python

import argparse
from Bio import AlignIO

#Convert a fasta file (aligned hopefully!) to phylip format.
# Uses BioPython to do everything. 
#
# Matt Gitzendanner
# University of Florida
#
# fasta_to_phy.py -i file.fa -o file.phy
#
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file")
parser.add_argument("-o", help="output file")

args = parser.parse_args()

infile = args.i
outfile = args.o

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

try:
	OUT=open(outfile, 'a')
except IOError:
		print "Can't open file", outfile
		
alignment = AlignIO.read(IN, "fasta")
AlignIO.write([alignment], OUT, "phylip-relaxed")
