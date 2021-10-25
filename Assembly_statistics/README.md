# Assembly statistic compilation scripts as of July 22nd 2021

`GNU parallel 20161222` is frequently used throughout these scripts. Information related to GNU parallel can be found [here](https://www.gnu.org/software/parallel/).

Scripts | Description | Software Dependencies | Additional Files | Notes
--------|-------------|-----------------------|------------------|------
SVs_and_ANI.sh | NA | `dnadiff version 1.3`, `fastANI version 1.33`, `seqtk version 1.3-r106` | run.list, (assembly summary refseq) **NEEDS TO BE ADDED** | pathways to files are hardcoded in run.list. See below.*
atcc_in_refseq_downloader_and_GC.sh | pulls data directly from RefSeq  | `BBMap version 38.18` | (assembly summary refseq) **NEEDS TO BE ADDED** | 
compare.all.levels.py | compares assembly statistics (Length, N50, GC Content, Contig Count) between ATCC assemblies and RefSeq assemblies | `python >=3.5` and `pandas` | tableS4.Refseq.assembly.summary.tsv and tableS4.atcc.assembly.summary.tsv | tableS4\* were made from internal product tracking or information from Refseq (with GC content added from bbmap)
quality_puller.py | downloads several assembly statistics for RefSeq assemblies of ATCC strains | `python >=3.5` with packages `Bio, xmltodict, pandas, argparse` | tableS4.Refseq.assembly.summary.tsv | *expect error message for header, will work for other rows*

### run.list Pathway Information
- ATCC assembly pathway: $pwd/fastas/atcc_fastas/$ref.fasta where $ref is the catalog number associated with the atcc assembly
- RefSeq assembly pathway: $pwd/fastas/refseq_fastas/$ext\*.fna where $ext is the GCF identified for the RefSeq assembly
