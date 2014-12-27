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
# -p parse orfs? Default= 0
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with GenBank accession IDs")
parser.add_argument("-e", help="email address")
parser.add_argument("-p", help="Parse ORFS? Default=0 (no)", default=0)

args = parser.parse_args()

infile = args.i
Entrez.email = args.e
ORFS= args.p

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
				try:  # Most CDS have a gene annotation.
					Gene= Feature.qualifiers["gene"]
					Gene= str(Gene[0]) # Gene is actually a list, so convert to string
					Gene = Gene[0].lower() + Gene[1:] #Convert 1st letter to lower case.
													# Had a lot of problems with some having atpA and others AtpA
													# This fixes that.
				except: # But if not, we need to handle them.
					try: # Many that are not annotated as gene are annotated with note
						Gene= Feature.qualifiers["note"]
						if type(Gene) is list:	#Some of these are lists, some are strings...
							Gene=str(Gene[0])
						Gene=Gene.replace("label: ","") #Clean up some where label: is part of the name.
						Gene = Gene[0].lower() + Gene[1:] #Convert 1st letter to lower case.
						
					except: #If we still can't get it, put in unknown.
						Gene= "unknown"
						print "Error parsing gene from CDS: ", Feature
				
				Gene=Gene.replace(" ", "_") #Clean up the name, replacing any spaces with underscores.
				
				if Gene[:3] != "orf" or ORFS==1: #if the gene name starts with orf, skip it unless user has opted to parse them.
					GeneFile= Gene + '.fna' #Name files by CDS name
				
					GeneFileAA= Gene + '.faa' #Name files by CDS name--Also get the translation
				
					try:
						OUT=open(GeneFile, 'a')
					except IOError:
						print "Can't open file to append", GeneFile
					
					try:
						OUTAA=open(GeneFileAA, 'a')
					except IOError:
						print "Can't open file to append", GeneFileAA
					
					try:	
						SeqNuc=Feature.extract(Sequence)		#Get the nucleotide sequence for the CDS
					except:
						print "Can't get sequence for ", Gene	#Handle problems.
					
					SeqAA=SeqRecord(Seq(Feature.qualifiers['translation'][0], IUPAC.protein), id=Sequence.id, description=Sequence.description)
					
					#For some reason many are returning Unknown IDs and Descriptions. Fix these.
					if SeqNuc.id == '<unknown id>':
						SeqNuc.id= Sequence.id
						SeqNuc.description = Sequence.description
										
					SeqIO.write(SeqNuc, OUT, "fasta")	# Write that to the file in fasta format.
					SeqIO.write(SeqAA, OUTAA, "fasta")
				else:
					print "Skipping ", Gene #Let the user know an orf was skipped.
