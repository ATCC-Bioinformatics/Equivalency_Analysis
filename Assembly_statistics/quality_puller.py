from Bio import Entrez
import xmltodict
import pandas as pd
import argparse as ap
parser = ap.ArgumentParser(description="Pull assorted assembly statistics for a list of assembly accessions of the form: GCA/GCF")
parser.add_argument('-l',required=True,help="File containing list of genbank or refseq assembly accessions -- GCF/GCA")
parser.add_argument('-p',required=True,help='Output Prefix')
parser.add_argument('-e',required=True,help='Email for Entrez to track')
parser.add_argument('-a',required=False,help='NCBI API Key - speeds up downloads by removing N queries per second limit')
myargs=parser.parse_args()

api_key = myargs.a
acc_list = list(pd.read_csv(myargs.l,header=Noneskiprows=1)[0])
prefix = myargs.p

Entrez.email=myargs.e


#initialize lists pre-dataframe
taxid = []
name = []
assembly_level = []
cn50 = []
sn50 = []
coverage = []
sub_date = []
up_date = []
source = []
contig_count = []
scaf_count = []
chrom_count = []
tot_len = []

for acc in acc_list:
    acc = acc.upper()
    if api_key is not None:
        id = Entrez.read(Entrez.esearch(db='assembly',term=acc,api_key=api_key))['IdList'][0] #should only have one id per assembly accession
        summary = xmltodict.parse(Entrez.esummary(db='assembly',id=id,api_key=api_key).read())['eSummaryResult']['DocumentSummarySet']['DocumentSummary']
    else:
        id = Entrez.read(Entrez.esearch(db='assembly',term=acc))['IdList'][0] #should only have one id per assembly accession
        summary = xmltodict.parse(Entrez.esummary(db='assembly',id=id).read())['eSummaryResult']['DocumentSummarySet']['DocumentSummary']
    print(summary['Taxid']) #progress bar
    taxid.append(summary['Taxid'])
    name.append(summary['Organism'])
    assembly_level.append(summary['AssemblyStatus'])
    cn50.append(summary['ContigN50'])
    sn50.append(summary['ScaffoldN50'])
    coverage.append(summary['Coverage'])
    sub_date.append(summary['SubmissionDate'])
    up_date.append(summary['LastUpdateDate'])
    source.append(summary['Biosource'])

    #get contig counts and scaffold counts ## <Stat category="contig_count" sequence_tag="all">24<
    contstart = summary['Meta'].find('contig_count') #find start position of contig count field
    count_end = summary['Meta'][contstart:10000].find('<')  #find end of contig count #number provided is just intended to be larger than it ever has to be. find returns position of first instance
    count_start = summary['Meta'][contstart:contstart+count_end].find('>') #find beginning of contig count
    contig_count.append(int(summary['Meta'][contstart+count_start+1:contstart+count_end])) #add 1 for indexing reasons

    scafstart = summary['Meta'].find('scaffold_count')
    scaf_end = summary['Meta'][scafstart:10000].find('<')
    scaf_count_start = summary['Meta'][scafstart:scafstart+scaf_end].find('>')
    scaf_count.append(int(summary['Meta'][scafstart+scaf_count_start+1:scafstart+scaf_end]))

    chromstart = summary['Meta'].find('chromosome_count')
    chrom_end = summary['Meta'][chromstart:10000].find('<')
    chrom_count_start = summary['Meta'][chromstart:chromstart+chrom_end].find('>')
    chrom_count.append(int(summary['Meta'][chromstart+chrom_count_start+1:chromstart+chrom_end]))

    lenstart = summary['Meta'].find('total_length')
    len_end = summary['Meta'][lenstart:10000].find('<')
    len_start = summary['Meta'][lenstart:lenstart+len_end].find('>')
    tot_len.append(int(summary['Meta'][lenstart+len_start+1:lenstart+len_end]))

df = pd.DataFrame()
df['Assembly Accession'] = acc_list #could initialize the dataframe this way from the beginning and loop over the column instead. Probably insignificant difference
df['TaxID'] = taxid
df['Organism'] = name
df['Assembly Level'] = assembly_level
df['Chromsome Count'] = chrom_count
df['Total Length'] = tot_len
df['Contig Count']  = contig_count
df['Contig N50'] = cn50
df['Scaffold Count'] = scaf_count
df['Scaffold N50'] = sn50
df['Read Coverage'] = coverage
df['Submission Date'] = sub_date
df['Latest Update'] = up_date
df['Bio Source'] = source

df.to_csv('{0}.assembly.summary.tsv'.format(prefix),sep='\t',index=False)

