ToolBox
=======

Various scripts for bioinformatics/phyloinformatics/etc that I've put together/stollen/borrowed over time

I've provided attribution where I can remember I took code bits from somewhere. If there's something that's your's I'm happy to add attribution or take it down if I accidentally violated the license.

See the [GenBankTools](GenBankTools/) folder for info on those scripts.

Some description of the scripts:

* **alignment_codon_parses.py** -- Create datasets by codon and 1st and 2nd position from an input nucleotide alignment.

* **BLASTPaser.pl** -- is a simple parser for BLAST output. The script is based on the [BioPerl SearchIO tutorial](http://www.bioperl.org/wiki/HOWTO:SearchIO).

* **BLASTParser_table.pl** -- is the same as above, but for tabular BLAST output.

* **clean_blank_taxa.pl** -- a script to remove taxa from a phylip file that have no data. Often when a multigene dataset is split into the individual genes, there are taxa which don't have data for some genes. These cause problems in the analysis, so this removs them.

* **compute_consensus_from_mafft.pl** -- take aligned fasta file (like a MAFFT ouput) and creats a single consensus sequence, using ambiguity codes in regions of overlap. Most of the code was stolen from part of a similar script by [Joseph Hughes](https://github.com/josephhughes/Sequence-manipulation/blob/master/Consensus.pl).

* **fasta_to_phy.pl** -- Convert a fasta file (aligned hopefully!) to (relaxed) phylip format.
 
* **MFAtoPHY.pl** -- a handy script from [Yu-Wei Wu](http://yuweibioinfo.blogspot.com/2009/01/fasta-to-phylip-converter.html) to convert a fasta file to a phylip file.

* **nex_to_fasta.py** -- A simple converter from Nexus to Fasta format using BioPython.

* **nex_to_phy.py** -- A simple converter from Nexus to relaxed-Phylip format using BioPython.

* **tally_occurrences.py** -- Tallies the occurrences of a list of species in an occurrence.txt file downloaded from GBIF.

* **seperate_genes_to_files.py** -- Takes an input fasta file with multiple genes from a species and splits that into gene files. Fasta records are assumed to be named by the gene name. Input file is assumed to be the species name. Appends record to gene file using the species name.


