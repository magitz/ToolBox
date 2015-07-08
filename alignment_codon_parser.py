#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq
import os

# Create datasets by codon and 1st&2nd from an input nucleotide alignment.
# The output file is in fasta format, but that can be easily converted to phylip with MFAtoPHY.pl
#
# Matt Gitzendanner
# University of Florida
# 
# version: 1.0: Initial release, July 8, 2015.
#
# Usage: python alignment_codon_parser.py -i alignment_file.phy
# Output datasets will be alignmnet_file.pos1.fna, alignmnet_file.pos2.fna, alignmnet_file.pos3.fna, and alignmnet_file.pos1and2.fna 

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input alignement file")
parser.add_argument("-f", help="Alignment format, fasta, phylip, nexus, etc. Use the file format string for BioPython. Default=phylip=relaxed.", default="phylip-relaxed")


args = parser.parse_args()

InFile = args.i
FileFormat= args.f

FileExtension=os.path.splitext(os.path.basename(InFile))[1]
FileBase=os.path.splitext(os.path.basename(InFile))[0]

for Record in SeqIO.parse(InFile, FileFormat):	
	SeqPos=[] 
	for x in range(3): #0,1,2
		SeqPos.append(Record.seq[int(x)::3]) #Get sequence for each position sequences
		#print ("Position %d is %s" %(x+1,SeqPos[x]))
			
		#Export the 3 single codon position datasets		
		OutFile=FileBase + ".pos" + str(x+1) + ".fna"
 		try:
 			OUT=open(OutFile, 'a')
 		except:
 			print("Can't open output file: %s" %(OutFile))
			
 		SeqIO.write(Record[int(x)::3], OUT, "fasta") #write this sequence
 		OUT.close #close the file.
		
	#Zip 1st and 2nd positions together
	FirstSecond=''.join([str(a)+b for a,b in zip(SeqPos[0],SeqPos[1])])
	#print ("Pos 1&2 is %s" %(FirstSecond))
	
	Record.seq=Seq(FirstSecond, IUPAC.ambiguous_dna)
	OutFile=FileBase + ".pos1and2.fna"
	try:
		OUT=open(OutFile, 'a')
	except:
		print("Can't open output file: %s" %(OutFile))
		
	SeqIO.write(Record, OUT, "fasta") #write this sequence
	OUT.close #close the file.
