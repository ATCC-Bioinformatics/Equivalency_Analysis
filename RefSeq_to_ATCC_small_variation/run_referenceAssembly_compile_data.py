from Bio import Entrez
import pandas as pd
import ssl
import os
import urllib
from urllib.error import HTTPError, URLError
from pandas.errors import EmptyDataError
import subprocess
import time
import glob
import io
import pysam
import xmltodict

def download_refseq(gcf_id):
    # Downloads RefSeq Assemblies based on a given GCF ID
    ids = Entrez.read(Entrez.esearch(db="assembly", term=gcf_id, retmax='200'))['IdList']
    print(f"found {len(ids)} ids")

    links = []
    for id in ids:
        #get summary
        summary = Entrez.read(Entrez.esummary(db="assembly", id=id, report="full"), validate=False)

        #get ftp link
        refseq_url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']

        if refseq_url == '':
            continue
    
        label = os.path.basename(refseq_url)
        link = os.path.join(refseq_url,label+'_genomic.fna.gz')

        file_path = f"public_assemblies/{label}.fna.gz"

        if os.path.exists(file_path):
            print("Reference already downloaded")
        else:
            print("Downloading...")
            while True:
                try:
                    print("Attempting to download Refseq fna...")
                    urllib.request.urlretrieve(link, file_path)
                    print(f'{label} downloaded sucessfully')
                    print('Unzipping file...')
                    unzip = subprocess.run(['gunzip', file_path, "-f"])
                    break
                except URLError as error:
                    print(f'Timeout error. Attempting to download {label} again')
                    time.sleep(10)
                except HTTPError as error:
                    print(f'Timeout error. Attempting to download {label} again')
                    time.sleep(10)

def download_gbk(gcf_id):
    # Downloads GenBank GBKs based on RefSeq GCF ID
    ids = Entrez.read(Entrez.esearch(db='assembly',term=gcf_id))['IdList']
    for id in ids:
        summary =  xmltodict.parse(Entrez.esummary(db='assembly',id=id,report='full').read())['eSummaryResult']['DocumentSummarySet']['DocumentSummary']
        ftp_main = summary['FtpPath_RefSeq']
        address = '{0}/{1}_genomic.gff.gz'.format(ftp_main.replace('ftp://', 'https://'), ftp_main.split('/')[-1])
        output_path = "gbk"

        filename = f"{ftp_main.split('/')[-1]}_genomic.gff"
        if os.path.exists(os.path.join(output_path, filename)):
            print("Gbk already downloaded")
        else:
            print(f"Downloading {filename}")
            os.system('wget {0} -P {1}'.format(address, output_path))
            os.system(f'gunzip gbk/{filename}.gz')
        
        return filename


list_161 = pd.read_csv(r"atcc_product_in_RefSeq.csv")
gcfs = list_161[["Illumina Filename", "GCF"]]

results = pd.DataFrame(columns = ["ATCC Catalog Number", "RefSeq Reference", "Average Coverage", "Variant Coverage", "# of SNPs", "# of Indels",
                                    "synonymous_variant", "missense_variant", "frameshift_variant", "frameshift_variant,start_lost", "frameshift_variant,stop_lost",
                                    "intergenic_variant", "start_lost", "stop_lost", "stop_gained", "stop_retained_variant", "-", "incomplete_terminal_codon_variant,coding_sequence_variant",
                                    "upstream_gene_variant", "downstream_gene_variant", "non_coding_transcript_exon_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "num_duplicates"])

# Run through ReferenceAssembly.py
assembly_script_path = "ReferenceAssembly.py"

all_variant_annotations = []

for i in range(len(gcfs)):
    acc_name = gcfs["GCF"].iloc[i].upper()
    download_refseq(acc_name)
    gbk_name = download_gbk(acc_name)

    ngs = gcfs["Illumina Filename"].iloc[i]
    file_name = ngs.split(".")[0] #Retrieve ngs_id without the file extension

    run_id = file_name[0:11]
    item_num = file_name.split("_")[2]
    print(item_num)

    r1_dir = f"{file_name}_R1.*"
    r2_dir = f"{file_name}_R2.*"
    

    r1_path = glob.glob(r1_dir)[0]
    r2_path = glob.glob(r2_dir)[0]    
    print(acc_name)
    ref = glob.glob(f"public_assemblies/{acc_name}*")[0]

    ###### Run referenceassembly.py
    os.system("mkdir -p results")
    output_path = f"results/{item_num}_{acc_name}"

    print("Running ReferenceAssembly.py....")
    run = os.system(f"python3 {assembly_script_path} -pe-1 {r1_path} -pe-2 {r2_path} -ref {ref} -out {output_path} -prefix {item_num} -gff gbk/{gbk_name} -threads 20")

    ## Run qualimap on samples
    os.system(f"qualimap bamqc -bam {output_path}/tmp/{item_num}.final.bam --outdir {output_path}/reads-to-reference_bamqc -outformat PDF")
    

    ###### Compile data
    bamqc_report_path = f"{output_path}/reads-to-reference_bamqc/genome_results.txt"
    qc_report_path = f"{output_path}/QCreport.txt"
    variant_ann_path = f"{output_path}/variant-ann.csv"
    
    cov = 0
    ## Read qualimap rSesults
    with open(bamqc_report_path, "r") as qc_reader:
        for lines in qc_reader: 
            if "mean coverageData" in lines:
                cov = lines.split("=")[1].strip()

    ## Parse variant annotations
    variant_ann = pd.read_csv(variant_ann_path) # Read in annotation csv generated by ReferenceAssembly.py
    variant_ann = variant_ann[variant_ann["FILTER"] == "PASS"] # Filter annotations that have passed
    num_duplicates = variant_ann.duplicated(subset="POS").sum()

    consequence = variant_ann["Consequence"].value_counts().to_dict()
    
    variant_types = ["synonymous_variant", "missense_variant", "frameshift_variant", "frameshift_variant,start_lost", "frameshift_variant,stop_lost",
                    "intergenic_variant", "start_lost", "stop_lost", "stop_gained", "stop_retained_variant", "-", "incomplete_terminal_codon_variant,coding_sequence_variant",
                    "upstream_gene_variant", "downstream_gene_variant", "non_coding_transcript_exon_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "num_duplicates"]

    # If annotations for the above variant types don't exist, set the count of that annotation type to 0
    ann_results = {}
    for t in variant_types:
        if t not in consequence.keys():
            ann_results[t] = 0
        else:
            ann_results[t] = consequence[t]
            
    try:        
    # Read bamqc qcreport.txt
        with open(qc_report_path, "r") as reader:
            file_list = reader.readlines()

            snp_line = file_list[0].strip()
            snp_count = snp_line.split(":")[1]

            indel_line = file_list[1].strip()
            indel_count = indel_line.split(":")[1]

            v_cov_line = file_list[2].strip()
            v_cov = v_cov_line.split(":")[1]

    except FileNotFoundError as error:            
        v_cov = "N/A"
        snp_count = "N/A"
        indel_count = "N/A"
  
    results.loc[len(results)] = [item_num, acc_name, cov, v_cov, snp_count, indel_count, 
                                ann_results["synonymous_variant"], ann_results["missense_variant"], ann_results["frameshift_variant"], ann_results["frameshift_variant,start_lost"], ann_results["frameshift_variant,stop_lost"], 
                                ann_results["intergenic_variant"], ann_results["start_lost"], ann_results["stop_lost"], ann_results["stop_gained"], ann_results["stop_retained_variant"], ann_results["-"], 
                                ann_results["incomplete_terminal_codon_variant,coding_sequence_variant"], ann_results["upstream_gene_variant"], ann_results["downstream_gene_variant"], ann_results["non_coding_transcript_exon_variant"], 
                                ann_results["3_prime_UTR_variant"], ann_results["5_prime_UTR_variant"], num_duplicates]


results.to_csv("list_of_161_ReferenceAssembly_results_with_gbk.csv", index=False)

