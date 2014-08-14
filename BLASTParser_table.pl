#!/usr/bin/perl

use strict;
use Getopt::Std;
use Data::Dumper;
#use Bio::SearchIO;
#use Bio::Tools::Run::Alignment::Clustalw;


#####################
# Options
#
# -i input file
# -t taxonomy file with species of each hit.
# -o otput file
#####################

my %options;
getopts('i:o:t:', \%options);

open INFILE, $options{i} or die "Problem opening input file";

open OUTFILE, ">$options{o}" or die "problem opening output file";
print OUTFILE "Query Name","\t","#of letters in query","\t","Target sequence","\t","Target gi", "\t","Target Species","\t","Hit length","\t","Fraction conserved","\t","Evalue","\n";

open TAXFILE, $options{t} or die "Problem opening taxonomy file";

my %taxonomy=();

while (<TAXFILE>)
{
	#line looks like this: scaffold-AEKF-2015794-Penium_margaritaceum	4097	root	Eukaryota	Viridiplantae	n	Streptophyta	n	n	n	asterids	n	Solanales	n	n	Solanaceae	Nicotianoideae	Nicotianeae	n	Nicotiana	n	Nicotiana tabacum	n	170340
	#									0								1		2			3		4				5	6				7	8	9	10			11	12			13	14	15			16				17			18	19			20	21					22	23
	# So we want the order, which should be element 13
	my @line_bits=split('\t', $_);
	if (($line_bits[13] ne "n") and ($line_bits[13] ne " "))
	{
		$taxonomy{$line_bits[0]}=$line_bits[13];
	}
	elsif (($line_bits[14] ne "n")  and ($line_bits[13] ne " "))#in case 13 is an n, move down to 14
	{
		$taxonomy{$line_bits[0]}=$line_bits[14];
	}
		elsif (($line_bits[15] ne "n")  and ($line_bits[13] ne " "))# or 15
	{
		$taxonomy{$line_bits[0]}=$line_bits[15];
	}
		elsif (($line_bits[16] ne "n")  and ($line_bits[13] ne " "))#or 16, which should be family
	{
		$taxonomy{$line_bits[0]}=$line_bits[16];
	}
	else
	{
		if (($line_bits[3] ne "n") and (defined($line_bits[3])))
		{
			$taxonomy{$line_bits[0]}=$line_bits[3]; # fall back to 3, which is kingdom
		}
		else {$taxonomy{$line_bits[0]}="Unknown";}
	}

}

while (<INFILE>)
{

	#line is: -outfmt '6 qseqid sseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore' 
	#						0		1	2	3		4		5		6		7		8	9		10	11		12
	#Also need to get taxonomy from all_taxonomy/CODE.taxonomy file with format
	
	my @line_bits=split('\t', $_);
	
	print OUTFILE $line_bits[0], "\t",
			$line_bits[4], "\t",
			$line_bits[1], "\t",
			$line_bits[2], "\t",
			$taxonomy{$line_bits[0]}, "\t",
			"\t", #no hsp_length--we could calculate from start-send, but we're not using anyway
			$line_bits[3], "\t",
			$line_bits[11], "\n";
}
