#!/usr/bin/env python
import os
import re
import string

import argparse

parser = argparse.ArgumentParser(description='Remove taxa from phylip file that have no real sequence data. i.e. The sequence is all -,N, or space characters.')
parser.add_argument('-i', action="store", dest="infile")
parser.add_argument('-o', action="store", dest="outfile")
args = parser.parse_args()

myfile=open(args.infile, 'r')
line=myfile.readline()

header=re.search('(\d+)\s*(\d+)', line)
tax_count=int(header.group(1))
chars=header.group(2)
keep_list=['']

for line in myfile:
	
	tax_line=re.search('(\S+)\s*(.*)', line)
	seq=tax_line.group(2)
	
	seq=seq.strip('\n')
	seq=re.sub('-','', seq)
 	seq=re.sub('N','', seq)
 	seq=re.sub('\s','', seq)
 	if (len(seq) > 0):
 		
 		keep_list.append(line)
 	else:
 		tax_count=tax_count-1
 		print 'No data for ', tax_line.group(1), ', removing from ', args.infile, '\n'
		
out=open(args.outfile, 'w')

myline=str(tax_count) + " " + chars + "\n"

out.write(myline)
for item in keep_list:
	out.write(item)

print args.infile, 'complete with ', tax_count, ' remaining taxa\n\n\n\n'

