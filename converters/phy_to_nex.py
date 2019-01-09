#!/usr/bin/env python

import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

# A simple converter from Phylip to Nexus format using BioPython.
# Matt Gitzendanner
# University of Florida

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input Phylip formatted file")
parser.add_argument("-o", help="output filename")
parser.add_argument("-a", help="Alphabet: dna or aa, default=dna", default="dna")

args = parser.parse_args()

infile = args.i
outfile = args.o
alphabet = args.a

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

try:
	OUT=open(outfile, 'a')
except IOError:
		print "Can't open file", outfile

if alphabet == "dna":		
	alignment = AlignIO.read(IN, "phylip-relaxed", alphabet=Gapped(IUPAC.ambiguous_dna))
	AlignIO.write([alignment], OUT, "nexus")

elif alphabet == "aa":		
	alignment = AlignIO.read(IN, "phylip-relaxed", alphabet=Gapped(IUPAC.protein))
	AlignIO.write([alignment], OUT, "nexus")
