#!/usr/bin/env python
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
import time

# =====================================================
#  Takes input of list of NCBI record IDs to fetch. 
#  Downloads these records and sorts into files based on 
#  description
#
#  Matt Gitzendanner
#  University of Florida
#  09/29/23
#
# Version 1.2: Update for python 3
# Version 1.1: Added option to sort (try to sort by gene) or not sort (just download everything)
# Version 1.0: Initial version based on GenBank_parser.py
# =====================================================


#####################
# Options
#
# -i Input file with list of NCBI record IDs.
# -e Email address used for Entrez
# -o Output directory for files.
# -s Try to sort by gene?
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with GenBank accession IDs")
parser.add_argument("-e", help="email address")
parser.add_argument("-o", help="directory for output files, needs to exists")
parser.add_argument("-s", help="Try to sort sequences by gene? (default, 1), or just download everything (0)", default=1)
args = parser.parse_args()

infile = args.i
Entrez.email = args.e
OutDir= args.o
SortByGene= args.s

try:
	IN=open(infile, 'r')
except IOError:
	print ("Can't open file", infile)

			 
for Line in IN:
	GInumber=Line.strip('\n')
	print ("Getting %s" %(GInumber))
	
	
	for i in range(3, 0, -1):
		try:
			GBSeq = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=GInumber) #Get the sequence
		
		except:
			if i == 1:
				raise
			print('Failed to connect. Retrying')
			time.sleep(5) #Wait 5 seconds and try again.
		else:
			break

	
	
	
	#The default is to attempt to sort the sequences by gene, using the series of pattern matches against 
	#  the sequence description below. If matched, the sequences are added to a fasta file with that 
	#  gene's name. If not matched, the sequence is added to a fasta file call All_others.fna. 
	
	# 9/29/23: Wow this is ugly! Sorry for my old self...
	
	if SortByGene == 1:
		for Sequence in SeqIO.parse(GBSeq, "gb"):			#Parse though each
			# Print some summary info about the sequence.
			print (Sequence.id, Sequence.description[:50] + "...")
			print ("Sequence length %i," % len(Sequence))
			print ("%i features," % len(Sequence.features))
			print ("from: %s" % Sequence.annotations["source"])
	
			if "trnL-trnF" in Sequence.description:
				GeneFile=os.path.join(OutDir, "trnL-trnF_intergenic_spacer.fna")
		
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
			elif "internal transcribed spacer" in Sequence.description or "ITS" in Sequence.description:
				GeneFile=os.path.join(OutDir,"ITS.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "external transcribed spacer" in Sequence.description or "ETS" in Sequence.description:
				GeneFile=os.path.join(OutDir,"ETS.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "18S" in Sequence.description:
				GeneFile=os.path.join(OutDir,"18S.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
			elif "26S" in Sequence.description:
				GeneFile=os.path.join(OutDir,"26S.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "atpA" in Sequence.description or "atpa" in Sequence.description:
				GeneFile=os.path.join(OutDir,"atpA.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "atpB" in Sequence.description or "atpb" in Sequence.description:
				GeneFile=os.path.join(OutDir,"atpB.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
			elif "atpF" in Sequence.description or "atpf" in Sequence.description:
				GeneFile=os.path.join(OutDir,"atpF.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
			elif "rbcL" in Sequence.description or "rbcl" in Sequence.description or "ribulose-1,5-biophosphate" in Sequence.description:
				GeneFile=os.path.join(OutDir,"rbcL.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
	
			elif "trnL" in Sequence.description:		#Just trnL, since we already took care of trnL-trnF
				GeneFile=os.path.join(OutDir,"trnL.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "matK" in Sequence.description or "maturase K" in Sequence.description:
				GeneFile=os.path.join(OutDir,"matK.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
			
			elif "ndhF" in Sequence.description or "ndhf" in Sequence.description:
				GeneFile=os.path.join(OutDir,"ndhF.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")

			elif "ndhA" in Sequence.description or "ndha" in Sequence.description:
				GeneFile=os.path.join(OutDir,"ndhA.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
					
			elif "rpl16" in Sequence.description or "rpL16" in Sequence.description:
				GeneFile=os.path.join(OutDir,"rpl16.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
			elif "rps16" in Sequence.description or "rpS16" in Sequence.description:
				GeneFile=os.path.join(OutDir,"rps16.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")

			elif "ycf1" in Sequence.description:
				GeneFile=os.path.join(OutDir,"ycf1.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")

			elif "ycf2" in Sequence.description:
				GeneFile=os.path.join(OutDir,"ycf2.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
		
			elif "psbA-trnH" in Sequence.description:
				GeneFile=os.path.join(OutDir,"psbA-trnH_spacer.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
			elif "petB-petD" in Sequence.description:
				GeneFile=os.path.join(OutDir,"petB-petD_spacer.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "trnC-pet1N" in Sequence.description:
				GeneFile=os.path.join(OutDir,"trnC-pet1N_spacer.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
	
			elif "trnS-trnG" in Sequence.description:
				GeneFile=os.path.join(OutDir,"trnS-trnG_spacer.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
	
			elif "trnY-trnT" in Sequence.description:
				GeneFile=os.path.join(OutDir,"trnY-trnT_spacer.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "LEAFY" in Sequence.description:
				GeneFile=os.path.join(OutDir,"LEAFY.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			elif "NIA" in Sequence.description:
				GeneFile=os.path.join(OutDir,"NIA.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")	
		
			elif "microsatellite" in Sequence.description:
				GeneFile=os.path.join(OutDir,"microsatellite.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
			else:
				GeneFile=os.path.join(OutDir,"All_others.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append", GeneFile)
		
				SeqIO.write(Sequence, OUT, "fasta")
		
	#SortByGene not == 1
	#If we are not trying to sort by gene using the pattern matches above, download each sequence
	#  in GenBank format to a folder named by the species, naming each file by the source.  
	else: 
		
		TaxonOutDir=os.path.join(OutDir, infile)  #Make the directory for the sequences, name by Taxon.
		try:
			if not os.path.exists(TaxonOutDir):
				os.makedirs(TaxonOutDir)
			print ("Downloading sequences to %s" %(TaxonOutDir))		
			
		except :
			print ("Can't make directory at: %s" %(TaxonOutDir))
		
		
		for Sequence in SeqIO.parse(GBSeq, "gb"):
			
			OutName= Sequence.annotations["organism"].replace(" ", "_") + "_" + Sequence.id
			OutFile=os.path.join(TaxonOutDir, OutName)
			try:
				OUT=open(OutFile, 'w')
			except IOError:
				print ("Can't open file to write: %s" %(OutFile))
	
			SeqIO.write(Sequence, OUT, "genbank")
