# Assembly statistic compilation scripts as of July 22nd 2021

`GNU parallel 20161222` is frequently used throughout these scripts. Information related to GNU parallel can be found [here](https://www.gnu.org/software/parallel/).

Scripts | Description | Software Dependencies | Additional Files
--------|-------------|-----------------------|-----------------
atcc_in_refseq_downloader_and_GC.sh | NA? | `BBMap version 38.18` | NA
quality_puller.py | downloads several assembly statistics for RefSeq assemblies of ATCC strains | `python >=3.5` with packages `Bio, xmltodict, pandas, argparse` | tableS4.Refseq.assembly.summary.tsv
SVs_and_ANI.sh | NA | `dnadiff version 1.3`, `fastANI version 1.33`, `seqtk version 1.3-r106` | atcc_in_refseq.tsv
compare.all.levels.py | compares assembly statistics (Length, N50, GC Content, Contig Count) between ATCC assemblies and RefSeq assemblies | `python >=3.5` and `pandas` | NA


All atcc assemblies in $pwd/fastas/atcc_fastas/$ref.fasta where $ref is the catalog number associated with the atcc assembly

All RefSeq assemblies in $pwd/fastas/refseq_fastas/$ext*.fna where $ext is the GCF identified for the RefSeq assembly

It runs on run.list, provided in this repo, which links ATCC assembly fastas (identified by a hash string) to the appropriate Refseq assemblies.


tableS4.atcc.assembly.summary.tsv and tableS4.Refseq.assembly.summary.tsv are available in this repo and were made from internal product tracking or information from Refseq (with GC content added from bbmap)
