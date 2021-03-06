# ToolBox

=======

Various scripts for bioinformatics/phyloinformatics/etc that I've put together/stollen/borrowed over time.

I've provided attribution where I can remember I took code bits from somewhere. If there's something that's yours I'm happy to add attribution or take it down if I accidentally violated the license.

The [GenBankTools](GenBankTools/) folder has scripts for interacting with GenBank records.

The [converters](converters/) folder has tools to convert files from one format to another. Most of these are sequence and alignment file converters.

The [BLASTTools](BLASTTools/) folder has tool for parsing or analyzing BLAST output.

The [SLURM_tools](SLURM_tools) folder has some scripts for managing and generating jobs under SLURM.

## Some description of the scripts in this directory

* **alignment_codon_parses.py** -- Create datasets by codon and 1st and 2nd position from an input nucleotide alignment.

* **clean_blank_taxa.pl** -- a script to remove taxa from a phylip file that have no data. Often when a multigene dataset is split into the individual genes, there are taxa which don't have data for some genes. These cause problems in the analysis, so this removs them.

* **compute_consensus_from_mafft.pl** -- take aligned fasta file (like a MAFFT ouput) and creats a single consensus sequence, using ambiguity codes in regions of overlap. Most of the code was stolen from part of a similar script by [Joseph Hughes](https://github.com/josephhughes/Sequence-manipulation/blob/master/Consensus.pl).

* **get_seqs_from_fasta.py** -- Get a set of sequences from an input fasta file and write to a new fasta file (can actually use different formats too, doesn't have to be fasta formatted).

* **remove_gaps.py** -- removes all gaps from sequences in an input fasta file.

* **remove_taxa.py** -- Removes taxa from a fasta file. Multiple taxa cna be designated with a comma separated list (no spaces).

* **separate_genes_to_files.py** -- Takes an input fasta file with multiple genes from a species and splits that into gene files. Fasta records are assumed to be named by the gene name. Input file is assumed to be the species name. Appends record to gene file using the species name.

* **tally_occurrences.py** -- Tallies the occurrences of a list of species in an occurrence.txt file downloaded from GBIF.

