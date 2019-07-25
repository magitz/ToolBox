#!/usr/bin/env python3

# Takes a fasta file with a set of names and another file with the same names
# plus additional information. This additional information is combined with 
# the name in the output fasta file.
#
# The additional information should be in the format:
# name <TAB> additional info
#
# Matt Gitzendanner
# University of Florida

import argparse
from Bio import SeqIO
import os


#Parse commandline options.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input fasta file")
parser.add_argument("-n", help="Input file with names to add to the current fasta name")
parser.add_argument("-o", help="Output file name")

args = parser.parse_args()


#Open the input fasta file
try:
    FASTA = open (args.i, 'r')
except OSError as exception:
    pass

#Open the input names file
try:
    NAMES = open (args.n, 'r')
except OSError as exception:
    pass

#Open the output fasta file
try:
    OUT = open (args.o, 'w')
except OSError as exception:
    pass

Name_dict={}
for Line in NAMES:
    Line.strip()
    Line_bits=Line.split('\t')
    Name_dict[Line_bits[0]]=Line_bits[1]

for record in SeqIO.parse(FASTA, "fasta") : #Read each record from the in_file
    thisID=record.id  #The name of the sequence.

    #try to find this name in the Name_dict
    try:
        newID = record.id + " " + Name_dict[thisID]
        record.id = newID
    except:
        print("No long description for {} ".format(thisID))

    SeqIO.write(record, OUT, "fasta")

