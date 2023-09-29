#!/usr/bin/env python

from Bio import Entrez
import os
import argparse
import time

# =====================================================
#  Takes a list of species names and queries GenBank for 
#  that species. If any data are in GenBank, a file is 
#  written that has the GenBank IDs for that species.
#
#  Matt Gitzendanner
#  University of Florida
#  9/29/23 Version 1.1: Update prints for python3
#  3/07/16 Version 1.0
#
# =====================================================


#####################
# Options
#
# -i input file with list of species.
# -e email address used for Entrez
# -o Output folder
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with GenBank accession IDs")
parser.add_argument("-e", help="email address")
parser.add_argument("-o", help="output folder name")

args = parser.parse_args()

infile = args.i
Entrez.email = args.e #sets the email for Entrez.
OutputFolder= args.o


#Function to test for non-zero file size.
def is_non_zero_file(fpath):  
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False
    
    
    
try:
	IN=open(infile, 'r')
except IOError:
	print ("Can't open file", infile)

for Line in IN:
	Line=Line.strip('\n')
	Organism= Line + "[Orgn]"
	
	OutFile=os.path.join(OutputFolder, Line.replace(" ", "_"))
	#Check if we've already done this species--look for non-zero output file. 
	#This allows rerunning to catch failed runs without redoing all. 
	if is_non_zero_file(OutFile):
		pass
	else:
	
		for i in range(3, 0, -1):
			try:
				GBSeq = Entrez.esearch(db="nucleotide", term=Organism ) #Get the sequence
				
			except:
				if i == 1:
					raise
				print('Failed to connect. Retrying')
				time.sleep(5) #Wait 5 seconds and try again.
			else:
				break

		
		Record= Entrez.read(GBSeq)
	
		if int(Record["Count"]) > 0:
			print ("%s had %d records in GenBank" %(Line, int(Record["Count"])))
		
			try:
			
				OUT=open(OutFile, 'w')
			except:
				print ("Can't open file: %s" %(OutFile))
		
			for id in Record["IdList"]:
				OUT.write(id + "\n")
