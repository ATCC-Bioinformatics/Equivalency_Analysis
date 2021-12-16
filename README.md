# Equivalency Analysis
This document provides additional documentation and supplemental data supporting our comparative analysis of bacterial reference genomes in NCBI's RefSeq database, as detailed in our recent pre-print:

> David A Yarmosh, Juan G Lopera, Nikhita P Puthuveetil, Patrick Ford Combs, Amy L Reese, Corina Tabron, Amanda E Pierola, James Duncan, Samuel R Greenfield, Robert Marlow, Stephen King, Marco A Riojas, John Bagnoli, Briana Benton, Jonathan L Jacobs ***"Comparative Analysis and Data Provenance for 1,113 Bacterial Genome Assemblies"*** bioRxiv 2021.12.14.472616; doi: https://doi.org/10.1101/2021.12.14.472616 (currently in peer review)

Additional details on our methods can be found in our [ASM Resource Announcements paper](https://journals.asm.org/doi/10.1128/MRA.00818-21), and in the [associated Git repo](https://github.com/ATCC-Bioinformatics/AGP-Resource-Announcement).

## Data Usage Agreement
All whole genome assemblies from the ATCC Genome Portal are avialable for Research use Only. Researchers who are interested in downloading individual genomes must
1. create an account on ATCC's main website (https://www.atcc.org)
2. agree to ATCC's Data Use Agreement when downloading data from the AGP. A copy of the Data Use Agreement can be found here: https://www.atcc.org/policies/product-use-policies/data-use-agreement.

## Raw FASTA/Q Data
Raw data for all ATCC Genome Portal assemblies is subject to the same Data Use Agreement described above. Bulk download of Illumina and Nanopore data is available from our MS Azure Blog storage space. The base URL is https://atccbioinformatics.blob.core.windows.net/publicdata, to which you would add a specific (case-sensitive) suffix for fastqs for both the ONT and Illumina from each assembly or the assembly suffix for a compressed fasta. The specific suffix list for all assemblies is maintained in a separate GitHub repo ([AGP-Raw-Data](https://github.com/ATCC-Bioinformatics/AGP-Raw-Data)) to be shared upon request.


## Code

*Please see individual Readme files within each folder for more details.*

Folder | Description | Scripts Within Folder
-------|-------------|----------------------
Assembly_statistics | contains scripts related to the creation of supplementary tables S3 and S4 where tableS3 represents assembly statistic to assembly statistic comparisons (i.e., Scaffold Count, N50, GC content) and table S4 represents structural variations between RefSeq assemblies and ATCC Standard Reference Genomes | SVs_and_ANI.sh, atcc_in_refseq_downloader_and_GC.sh, compare.all.levels.py, quality_puller.py, run.list
Database_survey | contains scripts related to the creation of table S2 (and all other tables). These download and parse metadata from RefSeq, ENSEMBL, and JGI. Only data derived from RefSeq was used in the publication | gather_atcc_nctc_assemblies.py, gather_metadata.py, parse_metadata.py
RefSeq_to_ATCC_small_variation | contains scripts related to the creation of table S5 that will compare ATCC Standard Reference Genome fastq files to their corresponding RefSeq assemblies and then parse the output | run_referenceAssembly_compile_data.py, ReferenceAssembly.py
