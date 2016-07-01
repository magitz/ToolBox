#!/usr/bin/env python
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
import sys
import time
import collections
import re

# =====================================================
#  Takes input of list of directories containing GenBank
#  records for the taxon named in the directory.
#  Summarizes how many of the sequences are plastid in
#  origin (source contains "plastid").
#  Also puts genes into files if it can, using the parsing
#  rules developed for the GenBank_getByID_sort.py script.
#
#  Matt Gitzendanner
#  University of Florida
#  03/21/16
#
# =====================================================


#####################
# Options
#
# -i Input directory
# -o Output directory for files
# -s Summary file name
# -m Min length for a sequence (default, 0=no cutoff)
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input directory with sequences")
parser.add_argument("-c", help="Config file with gene/region definitions")
parser.add_argument("-o", help="directory for output files, needs to exists")
parser.add_argument("-s", help="Name of summary file with counts for plastid and each gene")
parser.add_argument("-m", help="Minimum length for a seqeunce to be counted, default=0, no cutoff", default=0)
parser.add_argument("-l", help="Log file name")
parser.add_argument("-v", help="Verbosity, 0=very little, mostly just errors; 1=default; 2=lots", default=1)

args = parser.parse_args()

InDir = args.i
ConfigFile= args.c
OutDir= args.o
SummaryFile= args.s
MinLength= int(args.m)
LogFile= args.l
Verbose=int(args.v)


def GeneOccurence(Gene, SearchList):
	if any(Term in Sequence.description for Term in SearchList):
		GeneFile=os.path.join(OutDir, Gene + ".fna")
		
		for i in range(3, 0, -1):			#We are planning to run this in parallel and may have file contention from multiple runs trying to append
											#Hopefully just waiting a second will take care of this...
			try:
				OUT=open(GeneFile, 'a')
			except:
				if i == 1:
					raise
				print('Failed to open for appending. Retrying')
				time.sleep(1) #Wait 5 seconds and try again.
			else:
				break

		SeqIO.write(Sequence, OUT, "fasta")
		
		try:			#if we've seen this gene in this taxon, add one
			TaxonGeneCount[Gene]+=1
		except:			#otherwise, set count to 1
			TaxonGeneCount[Gene]=1		
		
		return 1
	else:
		return 0
				


try:
	LOG=open(LogFile, 'a')
except IOError:
	print ("Can't open file to write: %s" %(LogFile))
	sys.exit("Can't open file to write: %s" %(LogFile))
try:
	CONFIG=open(ConfigFile, 'r')
except IOError:
	print ("Can't open file to read: %s" %(ConfigFile))
	sys.exit("Can't open file to read: %s" %(ConfigFile))

Taxon=os.path.basename(os.path.normpath(InDir)) #Get the taxon from the directory name, removing trailing slash if needed

TaxonGeneCount={} #Setup dictionary to hold the gene counts for this taxon.

GeneOrdDict= collections.OrderedDict()	# GeneOrdDitct is an ordered dictionary where the Keys are the gene/region names and the 
										#  values are a list of search strings to match a gene/region in the GenBank
										#  record description field. They are processed in order--thus the need to use and OrderedDict--
										#  Values are read from config file.
SourceOrdDict= collections.OrderedDict()	# SourceOrdDict is an ordered dictionary to tally from the source field.

Source=0
Description=0
for Line in CONFIG:			#Parse the config file.
	Line=Line.strip('\n')
	if len(Line.strip()) == 0:		#Skip blank lines
		continue
	elif Line[0] == "#":		#Skip lines starting with #
		continue
		
	elif "[Source Tally]" in Line:
		Source=1
	
	elif "[Description Tally and Sort]" in Line:
		Description=1
		
	elif Source == 1 and Description == 0:			#Parse the [Source Tally] portion of the config file
		LineBits=re.split(':', Line)
	
		if len(LineBits) == 2:
			SourceKey=LineBits[0]
			SourceValue=re.split(r"(?<!\\),",LineBits[1])	#Use negative lookbehind assertion to not split on , preceeded by \.
			SourceValue= [i.strip() for i in SourceValue]
			SourceValue= [i.replace("\\","") for i in SourceValue]  #The escaping \ ends up duplicated in the string, remove both.
			
			if Verbose > 0:
				LOG.write("Adding source: %s \t\t Search String: %s\n" %(SourceKey, SourceValue))
				
			SourceOrdDict[SourceKey]= SourceValue
		else:
			print ("Error: Molformed line in config file: %s" %(Line))
			sys.exit("Error: Molformed line in config file: %s" %(Line))
	
	elif Source == 1 and Description == 1:		#Parse the [Description Tally and Sort] portion of the config file
		LineBits=re.split(':', Line)
		
		if len(LineBits) == 2:
			GeneKey=LineBits[0]
			GeneValue=re.split(r"(?<!\\),",LineBits[1])	#Use negative lookbehind assertion to not split on , preceeded by \.
			GeneValue= [i.strip() for i in GeneValue]
			GeneValue= [i.replace("\\","") for i in GeneValue]  #The escaping \ ends up duplicated in the string, remove both.
			if Verbose > 0:
				LOG.write("Adding Gene/region: %s \t\t Search String: %s\n" %(GeneKey, GeneValue))
	
			GeneOrdDict[GeneKey]= GeneValue
		else:
			print ("Error: Molformed line in config file: %s" %(Line))
			sys.exit("Error: Molformed line in config file: %s" %(Line))
		
	else:
		print ("Error: Molformed line in config file: %s" %(Line))
		sys.exit("Error: Molformed line in config file: %s" %(Line))
	



for File in os.listdir(InDir):		#open each file in the input directory.
	try:
		IN=open(os.path.join(InDir,File), 'r')
	except IOError:
		print ("Can't open file: %s" %(File))

	for Sequence in SeqIO.parse(IN, "gb"):			#Parse though the sequence file
		if not len(Sequence) > MinLength:
			print ("Sequence is only %d, below minimum length threshold, skipping" %(len(Sequence)))
			
		else:
			# Print some summary info about the sequence.
			if Verbose > 1:
				LOG.write(Sequence.id + "\t" + Sequence.description[:50] + "...")
				LOG.write(" Sequence length %i," % len(Sequence))
				LOG.write(" %i features," % len(Sequence.features))
				LOG.write(" from: %s" % Sequence.annotations["source"] + "\n")
			
			for k,v in SourceOrdDict.items():
				if any(Term in Sequence.annotations["source"].lower() for Term in v):
					try:		#Keep track of plastid genes, if we've seen one before, add 1
						TaxonGeneCount[k]+=1
					except:
						TaxonGeneCount[k]=1
			
			FoundGene=0
			
				
			for k,v in GeneOrdDict.items():  #Iterate through the ordered dictionary and look for the search strings.
				FoundGene= GeneOccurence(k, v)
				if FoundGene > 0: break
			
			if FoundGene == 0:	
				#If FoundGene is still 0, we haven't found it above, so write to All_others
				GeneFile=os.path.join(OutDir,"All_others.fna")
				try:
					OUT=open(GeneFile, 'a')
				except IOError:
					print ("Can't open file to append: %s" %(GeneFile))
	
				SeqIO.write(Sequence, OUT, "fasta")

try:
	SUMMARY=open(SummaryFile, 'a')
except IOError:
	print ("Can't open file to write: %s" %(SummaryFile))
	sys.exit("Can't open file to write: %s" %(SummaryFile))




SUMMARY.write(Taxon + ",")

LOG.write("\n\n#########  Header Line for summary file ####################\n\n")
LOG.write("Taxon,")

#print SourceOrdDict
for k, v in SourceOrdDict.items():
	try:
		LOG.write(str(k) + ",")
		SUMMARY.write(str(TaxonGeneCount[k]) + ",")
		
	except KeyError:
		LOG.write("??,")
		SUMMARY.write("0,")


print GeneOrdDict	
for k,v in GeneOrdDict.items():
	try:
		SUMMARY.write(str(TaxonGeneCount[k]) + ",")
	except KeyError:
		SUMMARY.write("0,")
	
	LOG.write(str(k) + ",")

LOG.write("\n\n######### END Header Line for summary file ####################\n")
SUMMARY.write("\n")


	
