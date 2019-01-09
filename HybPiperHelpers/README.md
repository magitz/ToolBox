# HybPiperHelpers

We are making more and more use of Matt Johnson et al.'s [HybPiper](https://github.com/mossmatters/HybPiper) pipeline.

This set of scripts is meant to help with running that pipeline.

* **prepTargetFile.py** Prepare fasta files of individual genes to be used in a HybPiper target file
   * Shortens name to first two underscore seperated items in record id (hopefully, genus_species)
   * Adds a gene number as passed in on command line: >genus_speces-G###
   * Removes reference sequence as indicated on command line
   * Resulting files can be concatinated into a single file for the target file.