#Collection of scripts involved in the assessment of small variants (SNPs/Indels) between refseq references of ATCC strains relative to the ATCC genome portal strain.

run_referenceAssembly_compile_data.py will run ReferenceAssembly.py in order to generate small variants between ACC read data and related RefSeq reference data and provide annotations such as synonymous/nonsynonymous mutation.

#File paths:#

Input files:

atcc_product_in_RefSeq.csv 

june2021_gcf_list.csv 

Raw reads: Available form Azure Blob

#Tools:# - Python 3.6.8 - Qualimap v.2.2.1 - BWA 0.7.17 - samtools 1.9 - GATK 4.0.8.1 - Requires Java - VEP 95.1 - tabix 1.10.2 - bcftools 1.9
