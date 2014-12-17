#!/usr/bin/env python
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# =====================================================
#  Takes input of list of NCBI record IDs to fetch. 
#  Downloads these records and parses the CDSs out, putting
#  them into fasta files named by the annotations.
#
#  Matt Gitzendanner
#  University of Florida
#  12/17/14
# =====================================================


#####################
# Options
#
# -i input file with list of NCBI record IDs.
# -e email address used for Entrez
#
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with GenBank accession IDs")
parser.add_argument("-e", help="email address")

args = parser.parse_args()

infile = args.i
Entrez.email = args.e

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile


for Line in IN:
	Line.strip('\n')
	print "Getting %s" %(Line)
	GBSeq = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=Line) #Get the sequence
	
	for Sequence in SeqIO.parse(GBSeq, "gb"):			#Parse though each
		# Print some summary info about the sequence.
		print Sequence.id, Sequence.description[:50] + "..."
		print "Sequence length %i," % len(Sequence),
		print "%i features," % len(Sequence.features),
		print "from: %s" % Sequence.annotations["source"]
		
		print "Parsing CDSs...\n\n"
		
		# Look at the CDSs and extract them
		for Feature in Sequence.features:
			if Feature.type == 'CDS' :
				Gene= Feature.qualifiers["gene"]
				GeneFile= str(Gene[0]) + '.fna' #Name files by CDS name
				
				GeneFileAA= str(Gene[0]) + '.faa' #Name files by CDS name--Also get the translation
				
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print "Can't open file to append", GeneFile
					
				try:
					OUTAA=open(GeneFileAA, 'a')
				except IOError:
					print "Can't open file to append", GeneFileAA
					
					
				SeqNuc=Feature.extract(Sequence)		#Get the nucleotide sequence for the CDS
				
				SeqAA=SeqRecord(Seq(Feature.qualifiers['translation'][0], IUPAC.protein), id=Sequence.id, description=Sequence.description)
					
				#For some reason many are returning Unknown IDs and Descriptions. Fix these.
				if SeqNuc.id == '<unknown id>':
					SeqNuc.id= Sequence.id
					SeqNuc.description = Sequence.description
										
				SeqIO.write(SeqNuc, OUT, "fasta")	# Write that to the file in fasta format.
				SeqIO.write(SeqAA, OUTAA, "fasta")
