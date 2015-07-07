#!/usr/bin/env python

# Takes a fasta file with multiple gene sequences for one species and splits it into 
# individual files for each gene. It appends to the gene files, so you can convert a
# set of files with multiple genes for each species into a set with the species in the 
# same file and multiple gene file. 
#
# Matt Gitzendanner
# University of Florida

import argparse
from Bio import SeqIO
import os


#Parse commandline options.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input fasta file")


args = parser.parse_args()
in_file = args.i

#Open the infile
try:
	IN = open (in_file, 'r')
except OSError as exception:
	pass
	
	
species_name=os.path.splitext(os.path.basename(in_file))[0] #Get the name of the input fasta file without the extension.
print species_name
	
for record in SeqIO.parse(IN, "fasta") : #Read each record from the in_file
	gene=record.id  #The gene name is the record.id of the sequence.
	record.id=species_name
	out_file=gene + ".fa"
	
	try:	
		OUT = open( out_file, 'a')	#append to the gene file.
	except OSError as exception:
		print ("Error opening output file for gene %s: %s" %(gene, out_file))
		pass
	
	SeqIO.write(record, OUT, "fasta")
	
	OUT.close



