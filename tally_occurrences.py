#!/usr/bin/env python

import argparse
import re


# Takes a list of taxa, goes through a GBIF occurrence file and tallies the occurrences of each taxon. 


#Parse commandline options.

parser = argparse.ArgumentParser(description='Tally occurrences from GBIF.')
parser.add_argument("-i", help="Input GBIF occurrence.txt file")
parser.add_argument("-t", help="Taxon list file")
parser.add_argument("-o", help="Output file name")
args = parser.parse_args()


in_file= args.i
taxon_file= args.t
out_file= args.o

#Open some files for reading and output
try:
	IN=open(in_file, 'r')
except IOError:
	print "Can't open input occurrence file", in_file

try:
	TAXON=open(taxon_file, 'r')
except IOError:
	print "Can't open taxon file", name_file


try:
	OUT=open(out_file, 'w')
except IOError:
	print "Can't open output file", out_file

print "Infile is: %s" %(in_file)
print "Namefile is: %s" %(taxon_file)
print "Output will go to: %s" %(out_file)

	
name_hash={}
missing_hash={}

for Line in TAXON:		#Read the new names into a hash with key being the species name and value set to 0
	Line = Line.strip('\n')
	name_hash[Line] = 0
	
for Line in IN: #Read through the occurrence.txt file, get species (column 219)
	Line = Line.strip('\n')
	Line_bits=re.split('\t', Line)
	
	try:
		name_hash[Line_bits[219]]+=1
	except:
		try:
			missing_hash[Line_bits[219]]+=1
		except:
			missing_hash[Line_bits[219]]=1

for key in name_hash:
	OUT.write(key + '\t' + str(name_hash[key]) + '\n')
	
for key in missing_hash:
	print "Not in list %s \t %d" %(key, missing_hash[key])
	
	

