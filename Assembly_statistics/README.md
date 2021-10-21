#Assembly statistic compilation scripts as of July 22nd 2021

GNU parallel is often employed - we used GNU parallel 20161222

#atcc_in_refseq_downloader_and_GC.sh requires:#
bbmap, specifically stats.sh - we used the conda installation of BBMap version 38.18

#quality_puller.py downloads several assembly statistics for RefSeq assemblies of ATCC strains. It requires:
ATCC-filtered assembly_summary_refseq.tsv
python >=3.5 with packages Bio, xmltodict, pandas, argparse

#SVs_and_ANI.sh requires:#

dnadiff installed - we used version 1.3

fastANI installed - we used version 1.33

seqtk installed - we used version 1.3-r106

atcc_in_refseq.tsv as input

All atcc assemblies in $pwd/fastas/atcc_fastas/$ref.fasta where $ref is the catalog number associated with the atcc assembly

All RefSeq assemblies in $pwd/fastas/refseq_fastas/$ext*.fna where $ext is the GCF identified for the RefSeq assembly

It runs on run.list, provided in this repo, which links ATCC assembly fastas (identified by a hash string) to the appropriate Refseq assemblies.

#compare.all.levels.py will compare assembly statistics (Length, N50, GC Content, Contig Count) between ATCC assemblies and RefSeq assemblies.#

It requres:

python >=3.5 and pandas

tableS4.atcc.assembly.summary.tsv and tableS4.Refseq.assembly.summary.tsv are available in this repo and were made from internal product tracking or information from Refseq (with GC content added from bbmap)
