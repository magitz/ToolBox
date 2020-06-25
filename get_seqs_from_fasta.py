#!/usr/bin/env python

import argparse
from Bio import SeqIO
import sys

# Pull specified sequences out of a fasta file.
#
# Matt Gitzendanner
# University of Florida
#
# get_seqs_from_fasta.py -i file.fa -o file.phy -s headers_of_sequences_to_get
#

__version__= 1.0   #  June 25, 2020

def open_files(in_file, out):
    """Opens the input and output files and returns file handles"""
    try:
        in_fh=open(in_file, 'r')
    except IOError:
        print (f"Can't open file for reading: {in_file}")
        sys.exit()
    try:
        out_fh=open(out, 'a')
    except IOError:
        print (f"Can't open file for writing: {out}")
        sys.exit()
    
    return in_fh, out_fh

def pull_seq(seqs_in, in_format, seqs_out, out_format, seq_list, verbose):
    """Pulls specified sequences from a file and write to a new file"""
    for record in SeqIO.parse(seqs_in, in_format):
        if record.id in seq_list:
            if verbose:
                print(f"Writing sequence for {record.id}")

            SeqIO.write(record, seqs_out, out_format)

def main():
    parser = argparse.ArgumentParser(prog='get_seqs_from_fasta.py')
    parser.add_argument("-i", "--in_file", help="input file")
    parser.add_argument("-f", "--in_format", help="File format of input file.",
                        default='fasta')
    parser.add_argument("-o", "--out_file", help="output file")
    parser.add_argument("-x", "--out_format", help="File format of output file.",
                        default='fasta')
    parser.add_argument('-s', '--sequences', nargs='+', default=[],
                         help="List of sequence IDs to keep")
    parser.add_argument('-v', '--verbose',action='store_true', default=False)
    parser.add_argument('--version', action='version', version=f'%(prog)s  Version: {__version__}:')

    args = parser.parse_args()

    in_fh, out_fh = open_files(args.in_file, args.out_file)

    pull_seq(in_fh, args.in_format, out_fh, args.out_format, args.sequences, args.verbose)

if __name__ == "__main__":
    main()



		
