analyze(){
    ref=$(echo $1 | cut -d\  -f1)
	echo -e "$ref \n\n"
    ext=$(echo $1 | cut -d\  -f2)

    seqtk seq -L 1000 fastas/atcc_fastas/$ref.fasta > fastas/atcc_greater_than_1kb/$ref.fasta
    seqtk seq -L 1000 fastas/refseq_fastas/$ext*.fna > fastas/refseq_greater_than_1kb/$ext.fna

    fastANI -r ./fastas/greater_than_1kb/$ref -q ./fastas/greater_than_1kb/$ext -o ./ANIresults/$ext.ani.txt
    dnadiff -p ./mummer_results/$ext ./fastas/atcc_greater_than_1kb/$ref.fasta ./fastas/refseq_greater_than_1kb/$ext.fna
}
export -f analyze

mkdir -p $pwd/ANIresults
mkdir -p $pwd/mummer_results/
mkdir -p $pwd/fastas/atcc_greater_than_1kb/
mkdir -p $pwd/fastas/refseq_greater_than_1kb/

cat run.list | parallel --jobs 30 "analyze {}"

echo -e "GenbankID\tRefBreakpoints\tQueryBreakpoints\tRefRelocations\tQueryRelocations\tRefTranslocations\tQueryTranslocations\tRefInversions\tQueryInversions\tReferenceBases\tQueryBases\tReferenceAlignedBases\tQueryAlignedBases\tRefTotalSeqs\tQueryTotalSeqs\tRefAlignedSeqs\tQueryAlignedSeqs\tRefUnalignedSeqs\tQueryUnalignedSeqs"  > dnadiff_greater_than_1kb_summary.tsv
for i in $(find mummer_results -name *report); do
    name=$(basename $i .report)
    rbases=$(grep TotalBases $i | sed 's/ \+/\t/g' | cut -f2)
    qbases=$(grep TotalBases $i | sed 's/ \+/\t/g' | cut -f3)
    rbp=$(grep Breakpoint $i | sed 's/ \+/\t/g'|cut -f2)
    rreloc=$(grep Relocation $i | sed 's/ \+/\t/g'|cut -f2)
    rtrans=$(grep Translocation $i | sed 's/ \+/\t/g'|cut -f2)
    rinv=$(grep Inversion $i | sed 's/ \+/\t/g'|cut -f2)
    qbp=$(grep Breakpoint $i | sed 's/ \+/\t/g'|cut -f3)
    qreloc=$(grep Relocation $i | sed 's/ \+/\t/g'|cut -f3)
    qtrans=$(grep Translocation $i | sed 's/ \+/\t/g'|cut -f3)
    qinv=$(grep Inversion $i | sed 's/ \+/\t/g'|cut -f3)
    rseq=$(grep TotalSeq $i | sed 's/ \+/\t/g' | cut -f2)
    qseq=$(grep TotalSeq $i | sed 's/ \+/\t/g' | cut -f3)
    raln=$(grep AlignedSeqs $i | sed 's/ \+/\t/g' | cut -f2| cut -d\( -f2 | cut -d\% -f1)
    runaln=$(grep UnalignedSeqs $i | sed 's/ \+/\t/g' | cut -f2| cut -d\( -f2 | cut -d\% -f1)
    ralnb=$(grep AlignedBases $i | sed 's/ \+/\t/g' | cut -f2| cut -d\( -f2 | cut -d\% -f1)
    runalnb=$(grep UnalignedBases $i | sed 's/ \+/\t/g' | cut -f2| cut -d\( -f2 | cut -d\% -f1)
    qaln=$(grep AlignedSeqs $i | sed 's/ \+/\t/g' | cut -f3| cut -d\( -f2 | cut -d\% -f1)
    qunaln=$(grep UnalignedSeqs $i | sed 's/ \+/\t/g' | cut -f3| cut -d\( -f2 | cut -d\% -f1)
    qalnb=$(grep AlignedBases $i | sed 's/ \+/\t/g' | cut -f3| cut -d\( -f2 | cut -d\% -f1)
    qunalnb=$(grep UnalignedBases $i | sed 's/ \+/\t/g' | cut -f3| cut -d\( -f2 | cut -d\% -f1)
    echo -e "$name\t$rbp\t$qbp\t$rreloc\t$qreloc\t$rtrans\t$qtrans\t$rinv\t$qinv\t$rbases\t$qbases\t$ralnb\t$qalnb\t$rseq\t$qseq\t$raln\t$qaln\t$runaln\t$qunaln"
done >> dnadiff_greater_than_1kb_summary.tsv
