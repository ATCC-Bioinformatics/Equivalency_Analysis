#Collection of scripts involved in the assessment of small variants (SNPs/Indels) between refseq references of ATCC strains relative to the ATCC genome portal strain fastq files.

run_referenceAssembly_compile_data.py will run ReferenceAssembly.py in order to generate small variants between ACC read data and related RefSeq reference data and provide annotations such as synonymous/nonsynonymous mutation.

#File paths:#

Input file:

atcc_product_in_RefSeq.csv contains columns Item Name - the ATCC Product Catalog Number, GCF - the RefSeq accession number, Illumina Filename - the basename for both R1/R2 fastq files, and ONT Filename - the OxfordNanopore filename.

Raw reads: Available from Azure Blob

#Tools:# - Python 3.6.8 - Qualimap v.2.2.1 - BWA 0.7.17 - samtools 1.9 - GATK 4.0.8.1 - Requires Java - VEP 95.1 - tabix 1.10.2 - bcftools 1.9
