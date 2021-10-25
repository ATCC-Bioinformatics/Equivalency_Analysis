## Assessment of small variants (SNPs/Indels) between RefSeq sequences labeled as ATCC strains (reference) relative to the ATCC genome portal strain fastq files (query). 

Scripts | Description | Software Dependencies
--------|-------------|----------------------
run_referenceAssembly_compile_data.py | will run ReferenceAssembly.py in order to generate small variants between ATCC read data and related RefSeq reference data and provide annotations such as synonymous/nonsynonymous mutation | `python version 3.6.8`, `Qualimap v.2.2.1`, `BWA version 0.7.17`, `samtools version 1.9`, `GATK version 4.0.8.1`, `Java`, `VEP version 95.1`, `tabix version 1.10.2`, `bcftools version 1.9`

Additional Files | Notes
-----------------|------
| ReferenceAssembly.py, atcc_product_in_RefSeq.csv | atcc_product_in_RefSeq.csv contains columns: 1. Item Name: the ATCC Product Catalog Number, 2. GCF: the RefSeq accession number, 3. Illumina Filename: the basename for both R1/R2 fastq files, 4. and ONT Filename: the OxfordNanopore filename.
