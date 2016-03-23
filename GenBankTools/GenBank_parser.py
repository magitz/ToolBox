#!/usr/bin/env python
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
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
#
# Version 1.2: 2/10/16
#   -Added option to include a certain number of base pairs upstream and downstream of the annotated gene.
#       Note that in this may include other annotated reagions.
# Version 1.1: 1/5/15
#	-Added handling of multiple sequences with the same annotation.
#		By default genes with the same name in a sample will have a _copy_# added to everything after the 1st copy.
#		Set to 1 to exclude subsequent copies from the results files. Only the first is reported.
#	-Added testing for different cases of gene names. Uses the first one
#		for all subsequent samples, so pick the first sample carefully.
#		e.g. if the first sample has ccmFC and others have ccmFc, they will
#		be added to the ccmFC file.
# =====================================================


#####################
# Options
#
# -i input file with list of NCBI record IDs.
# -e email address used for Entrez
# -n Name sample with description? Takes description field, removes spaces and adds to id field.
# -p parse orfs? Default= 0
# -d exclude duplicate copies? Default=0
# -f Include flanking region +/- # of bp. Default=0
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with GenBank accession IDs")
parser.add_argument("-e", help="email address")
parser.add_argument("-n", help="Name sample with description? Takes description field, removes spaces and adds to id field. Default= 1 (yes), use 0 for no.", default=1)
parser.add_argument("-p", help="Parse ORFS? Default=0 (no)", default=0)
parser.add_argument("-d", help="Exclude duplicate copies of genes? Default=0 (no)", default=0)
parser.add_argument("-f", help="Include flanking region +/- # of bp. Default=0", default=0)

args = parser.parse_args()

infile = args.i
Entrez.email = args.e
ORFS= args.p
AddDescription= int(args.n)
ExclDups= int(args.d)
Flanking= int(args.f)

# Function to get overlap in intervals. From: http://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
def getOverlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def addDescriptionToName(id, description):
	description_bits=description.split()[0:2] #get 1st 2 words from description (Geneus and species)
	id= "_".join(description_bits) + "_" + id

	return id

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

Gene_dict={} # Keep track of all the genes we've found. Keys are lowercase version, 
			 #  values are the case as found in the first sample where the gene is found.
			 
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
		
		Sample_gene_dict={} #Clear the dict for tracking genes found in each sample.
		# Look at the CDSs and extract them
		for Feature in Sequence.features:
			if Feature.type == 'CDS' :
				try:  # Most CDS have a gene annotation.
					Gene= Feature.qualifiers["gene"]
					Gene= str(Gene[0]) # Gene is actually a list, so convert to string
					
				except: # But if not, we need to handle them.
					try: # Many that are not annotated as gene are annotated with note
						Gene= Feature.qualifiers["note"]
						if type(Gene) is list:	#Some of these are lists, some are strings...
							Gene=str(Gene[0])
						Gene=Gene.replace("label: ","") #Clean up some where label: is part of the name.
						
					except: #If we still can't get it, put in unknown.
						Gene= "unknown"
						print "Error parsing gene from CDS: ", Feature
				
				Gene=Gene.replace(" ", "_") #Clean up the name, replacing any spaces with underscores.
				
				if Gene.lower() not in Gene_dict.keys():
					Gene_dict[Gene.lower()]=Gene		#if we haven't found this gene yet, add it ot the Gene_dict
														#with key as lowercase version and value as slightly cleaned up version we found.
				try:
					Sample_gene_dict[Gene.lower()] += 1 #Increment if we've already seen this gene in this sample.
				except KeyError:
					Sample_gene_dict[Gene.lower()] = 1 #If we haven't seen this gene in this sample yet, add it to the Sample_gene_dict
					
				if Gene[:3].lower() != "orf" or ORFS==1: #if the gene name starts with orf, skip it unless user has opted to parse them.
					GeneFile= Gene_dict[Gene.lower()] + '.fna' #Name files by CDS name using the case as found the first time the gene was found.
				
					GeneFileAA= Gene_dict[Gene.lower()] + '.faa' #Name files by CDS name--Also get the translation
				
					try:
						OUT=open(GeneFile, 'a')
					except IOError:
						print "Can't open file to append", GeneFile
					
					try:
						OUTAA=open(GeneFileAA, 'a')
					except IOError:
						print "Can't open file to append", GeneFileAA
					
					############
					# Manage getting flanking regions.
					############
					#print ("Feature location is: %s" %(Feature.location))
					
					Flanking_start=Feature.location.start.position - Flanking #Subtract from start
					if Flanking_start < 1:		#Handle going past the beginning of the sequence
						Flanking_start = 1
						
					Flanking_end=Feature.location.end.position + Flanking #Add to end
					if Flanking_end > len(Sequence):
						Flanking_end = len(Sequence)	#Handle going past the end of the Sequence
						
					Flanking_location= SeqFeature.FeatureLocation(Flanking_start, Flanking_end, Feature.location.strand)
					Feature.location = Flanking_location
					#print ("With flanking, feature location is: %s" %(Feature.location))
					
					try: del Add_to_name #Clear the Add_to_name variable
					except: pass
						
					# Check if the expanded Feature.location now includes another CDS
					for OtherGenes in Sequence.features:
						if OtherGenes.type == 'CDS' :
							try:  # Most CDS have a gene annotation.
								OtherGene= OtherGenes.qualifiers["gene"]
								OtherGene= str(OtherGene[0]) # Gene is actually a list, so convert to string
							except:
								pass

							if not Gene.lower() == OtherGene.lower():
													
								Overlap= getOverlap([Flanking_start, Flanking_end], [OtherGenes.location.start.position, OtherGenes.location.end.position])
							
								if Overlap > 0:
									#print ("Expanded flanking region for %s includes %s, adding to name of file" %(Gene.lower(), OtherGene))
									try: Add_to_name=Add_to_name + "_" + OtherGene.lower()
									except: Add_to_name= "_" + Gene + "_with_" + OtherGene.lower()
					
					#Report what genes were added to the resulting files by expanding the flanking region.
					try:
						print ("For gene %s, flanking region includes %s" %(Gene, Add_to_name))			
					except:
						pass
						
						
					try:	
						SeqNuc=Feature.extract(Sequence)		#Get the nucleotide sequence for the CDS	
					except:
						print "Can't get sequence for ", Gene	#Handle problems.
					
					SeqAA=SeqRecord(Seq(Feature.qualifiers['translation'][0], IUPAC.protein), id=Sequence.id, description=Sequence.description)
					
					#For some reason many are returning Unknown IDs and Descriptions. Fix these.
					if SeqNuc.id == '<unknown id>':
						SeqNuc.id= Sequence.id
						SeqNuc.description = Sequence.description
					
					if Gene != "unknown" and (Sample_gene_dict[Gene.lower()] > 1 and ExclDups == 0):	#If there's multiple copies, and we want them reported (skip unknowns since they aren't likely the same gene
						
						#Add _copy_# to the end of the name	
						SeqNuc.id = SeqNuc.id + "_copy_" + str(Sample_gene_dict[Gene.lower()])
						SeqAA.id = SeqAA.id	+ "_copy_" + str(Sample_gene_dict[Gene.lower()])
						
						SeqNuc.id= addDescriptionToName(SeqNuc.id, SeqNuc.description)
						SeqAA.id= addDescriptionToName(SeqAA.id, SeqAA.description)
						
						
						try:		#If we have added other regions vy including flanking region, add these to the name of the sequence so user knows.
									#Only do this for nuclear data as AA data is just the original gene.
							SeqNuc.id = SeqNuc.id + Add_to_name
						except: pass

						
						
						SeqIO.write(SeqNuc, OUT, "fasta")	# Write that to the file in fasta format.
						SeqIO.write(SeqAA, OUTAA, "fasta")
					
					elif Gene != "unknown" and (Sample_gene_dict[Gene.lower()] > 1 and ExclDups == 1):	#If there's multiple copies and we don't want them reported, let the user know.
						print "Copy %s of %s being excluded for %s" %(str(Sample_gene_dict[Gene.lower()]), Gene_dict[Gene.lower()], SeqNuc.id)
						
					else:
						
						SeqNuc.id= addDescriptionToName(SeqNuc.id, SeqNuc.description)
						SeqAA.id= addDescriptionToName(SeqAA.id, SeqAA.description)
						
						try:		#If we have added other regions vy including flanking region, add these to the name of the sequence so user knows.
									#Only do this for nuclear data as AA data is just the original gene.
							SeqNuc.id = SeqNuc.id + Add_to_name
						except: pass

						SeqIO.write(SeqNuc, OUT, "fasta")	# Write that to the file in fasta format.
						SeqIO.write(SeqAA, OUTAA, "fasta")
				else:
					print "Skipping ", Gene #Let the user know an orf was skipped.
