### Information about source files downloaded on 06 June 2021
File | Source 
:----- | :----: 
assembly_summary_refseq.txt | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
atcc2nctc.csv | adapted from https://www.phe-culturecollections.org.uk/products/bacteria/atcc-equivalents.aspx
### Scripts
#### gather_atcc_nctc_assemblies.py
This script processes the Bacterial RefSeq assembly summary file and keeps all records that contained the keyword atcc or nctc. 
#### gather_metadata.py 
This script processes the output from the above script, downloads each assembly's assembly_summary.txt report, and pulls out the metadata from that report. 
#### parse_metadata.py 
This script parses the pickled metadata object created in the above script and creates a tab-delimited file with one line for each assembly containing all available metadata.
