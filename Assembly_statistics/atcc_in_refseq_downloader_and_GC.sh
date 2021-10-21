#Given assembly_summary_refseq.txt that has been filtered for ATCC genomes:
jobs=32 #define how many threads to use in download
cut -f 20 ATCC_summary_refseq.tsv > atcc.summary.tmp #get ftp address field
rev atcc.summary.tmp | cut -d\/ -f1 | rev > atcc.names.tmp #get filename prefix from address
while read line; do echo "_genomic.fna.gz"; done < atcc.names.tmp > atcc.genomic.tmp #set filename suffix
paste -d'\0' atcc.names.tmp atcc.genomic.tmp > atcc.addr.tmp #combine filename prefix and suffix with no delimiter
paste -d\/ atcc.summary.tmp atcc.addr.tmp | sed 's|ftp://|https://|g' > atcc.dl.txt #combine ftp address with filename by '/' delimiter to create address. Replace ftp with https
rm atcc*tmp #get rid of temporary files

cat atcc.dl.txt | parallel --jobs $jobs "wget {}" #download all files in parallel according to how many threads chosen at the top

mkdir -p gc_content
#use stats.sh from bbmap to generate gc content on many files, parse results into single file
for i in *fna; do stats.sh $i > gc_content/$(basename $i .fna); done #fast enough that parallelization adds little
#get prefix based on refseq name -- GCF_000006805.1_ASM680v1_genomic.fna, pull assembly gc content from bbmap output file. Print all to single file
for i in gc_content/*; do name=$(basename $i | cut -d\_ -f1,2); gc=$(head -n2 $i | tail -n1 | cut -f8); echo -e "$name\t$gc"; done > gc.content.bbmap.refseq