#!/usr/bin/perl

use strict;
use Getopt::Std;
use Data::Dumper;
use Bio::SearchIO;
#use Bio::Tools::Run::Alignment::Clustalw;


#####################
# Options
#
# -i input file
# -o otput file
#####################

my %options;
getopts('i:o:', \%options);

my $sequence_number = 1;

open INFILE, $options{i} or die "Problem opening input file";

open OUTFILE, ">$options{o}" or die "problem opening output file";
print OUTFILE "Query Name","\t","#of letters in query","\t","Target sequence","\t","Target gi", "\t","Target Species","\t","Hit length","\t","Fraction conserved","\t","Evalue","\n";

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $options{i});
while( my $result = $in->next_result) {
 my $hit=$result->next_hit;
 if (! defined ($hit)) { next;}
 my $hsp=$hit->next_hsp;
 
 my $hit_descr=$hit->description;
 my $hit_spp=$hit_descr;
 $hit_spp=~s/.*\[(.*)\].*/$1/;
 
 my $hit_gi=$hit->name;
 $hit_gi=~s/gi\|(\d*).*/$1/;
 
 print OUTFILE $result->query_name, "\t",
 		$result->query_length, "\t",
 		$hit_descr, "\t",
 		$hit_gi,"\t",
 		$hit_spp, "\t",
 		$hsp->hsp_length, "\t",
 		$hsp->frac_conserved, "\t",
 		$hsp->evalue,"\n";
 		
}
print $sequence_number;

