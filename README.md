## Data Usage Agreement
All whole genome assemblies from the ATCC Genome Portal are avialable for Research use Only. Researchers who are interested in downloading individual genomes must 1/ create an account on ATCC's main website (https://www.atcc.org), and 2/ agree to ATCC's Data Use Agreement when downloading data from the AGP. A copy of the Data Use Agreement can be found here: https://www.atcc.org/policies/product-use-policies/data-use-agreement.

## Raw FASTQ Data
Raw data for all ATCC Genome Portal assemblies is subject to the same Data Use Agreement described above. Bulk download of Illumina and Nanopore data is available from our MS Azure Blog storage space. The base URL is https://atccbioinformatics.blob.core.windows.net/publicdata, to which you would add a specific suffix for a compressed tarball for both the ONT and Illumina from each assembly. The specific suffix list for all assemblies is maintained in a separate GitHub repo ([AGP-Raw-Data](https://github.com/ATCC-Bioinformatics/AGP-Raw-Data)) to be shared upon request.

This repository contains the various scripts and input files used in ATCC's Equivalency Analysis (name to be updated with final publicaiton name).

There are three folders containing scripts related to the compilation and generation of various data.

1. Assembly_statistics contains scripts related to the creation of supplementary tables S3 and S4 where tableS3 represents assembly statistic to assembly statistic comparisons (i.e., Scaffold Count, N50, GC content) and table S4 represents structural variations between RefSeq assemblies and ATCC Standard Reference Genomes.

2. Database_survey contains scripts related to the creation of table S2 (and all other tables). These download and parse metadata from RefSeq, ENSEMBL, and JGI. Only data derived from RefSeq was used in the publication.

3. RefSeq_to_ATCC_small_variation contains scripts related to the creation of table S5 that will compare ATCC Standard Reference Genome fastq files to their corresponding RefSeq assemblies and then parse the output.

Please see individual Readme files within each folder for more details.

DOI:
