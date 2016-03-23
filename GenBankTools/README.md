ToolBox: GenBankTools
=======

Various scripts mostly scripts for querying GenBank or analyzing GenBank files to extract data.

* **Check_for_species_in_GenBank.py** -- Takes a list of species names and queries GenBank for that species. If any data are in GenBank, a file is written that has the GenBank IDs for that species. 

* **GenBank_getByID_sort.py** -- Takes input of list of NCBI record IDs to fetch. Downloads these records and sorts into files based on description.

* **GenBank_parser.py** --Takes input of list of NCBI record IDs to fetch. Downloads these records and parses the CDSs out, putting them into fasta files named by the annotations.
 
* **GenBank_summarize.py** -- Takes input of list of directories containing GenBank records for the taxon named in the directory. Summarizes how many of the sequences are plastid in origin (source contains "plastid"). Also puts genes into files if it can, using the parsing rules developed for the GenBank_getByID_sort.py script. See also the SummarizeExampleDictionary.txt as an example of the input config file.