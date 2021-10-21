#! /usr/bin/python

# The prokaryotes table contains replicon information not found in the assembly
# summary tables
prokaryotes_accession_replicon = {}
prokaryotes_accession_level = {}
prokaryotes_accession_line = {}
with open('prokaryotes.csv','r') as f:
    c = 0
    for line in f:
        if c == 0:
            c+=1
            pass
        else:
            line = line.lower()
            line_arr = line.split(',')
            accession = line_arr[5].replace("\"","")
            # line_arr 9 contains replicon information (plasmid)
            prokaryotes_accession_replicon[accession] = line_arr[9]
            prokaryotes_accession_level[accession] = line_arr[6].replace("\"","")
            prokaryotes_accession_line[accession] = line

# genbank data collected for use with Ensembl and JGI tables
# refseq is collected separately below
genbank_accession_level = {}
genbank_accession_line = {}
with open('assembly_summary_genbank.txt','r',encoding='utf-8') as f:
    c = 0
    for line in f:
        if c < 2:
            c+=1
            pass
        else:
            line = line.lower()
            line_arr = line.split('\t')
            gca = line_arr[0]
            gcf = line_arr[17]
            level = line_arr[11]
            # Store the genbank data for use with other databases
            genbank_accession_level[gca] = line_arr[11]
            genbank_accession_line[gca] = line

#############################      NCBI Data       #############################

def search_ncbi_refseq():
    # Store data to search for ATCC keyword later
    refseq_accession_level = {}
    refseq_accession_line = {}
    genbank2refseq = {}
    # Counts for whole dataset
    total = 0
    complete_or_chromosome = 0
    contigs_or_scaffolds = 0
    complete_or_chromosome_and_plasmids = 0
    # Counts for atcc strains in dataset
    atcc_total = 0
    atcc_complete_or_chromosome = 0
    atcc_contigs_or_scaffolds = 0
    atcc_complete_or_chromosome_and_plasmids = 0

    with open('assembly_summary_refseq.txt','r',encoding='utf-8') as f:
        c = 0
        for line in f:
            if c < 2:
                c+=1
                pass
            else:
                line = line.lower()
                line_arr = line.split('\t')
                gcf = line_arr[0]
                gca = line_arr[17]
                level = line_arr[11]

                # Store the refseq data for use with other databases
                refseq_accession_level[gcf] = line_arr[11]
                refseq_accession_line[gcf] = line
                genbank2refseq[gca] = gcf
                prokaryotes_line = ''
                if gca in prokaryotes_accession_line.keys():
                    prokaryotes_line = prokaryotes_accession_line[gca]
                replicon_info = ''
                if gca in prokaryotes_accession_replicon.keys():
                    replicon_info = prokaryotes_accession_replicon[gca]

                # If the accession number is one included in the Prokaryotic RefSeq
                if level == "complete genome":
                    if "plasmid" in replicon_info:
                        complete_or_chromosome_and_plasmids+=1
                        if "atcc" in line or 'atcc' in prokaryotes_line:
                            # Increment assembly-type count
                            atcc_complete_or_chromosome_and_plasmids+=1
                            atcc_total+=1
                    else:
                        complete_or_chromosome+=1
                        if "atcc" in line or 'atcc' in prokaryotes_line:
                            # Increment assembly-type count
                            atcc_complete_or_chromosome+=1
                            atcc_total+=1
                    total+=1
                elif level == "chromosome":
                    complete_or_chromosome+=1
                    if "atcc" in line or 'atcc' in prokaryotes_line:
                        # Increment assembly-type count
                        atcc_complete_or_chromosome+=1
                        atcc_total+=1
                    total+=1
                elif level == "contig" or level == "scaffold":
                    contigs_or_scaffolds+=1
                    if "atcc" in line or 'atcc' in prokaryotes_line:
                        # Increment assembly-type count
                        atcc_contigs_or_scaffolds+=1
                        atcc_total+=1
                    total+=1

    ## Data for table
    print("RefSeq Data:")
    print("Total:",total,"\nContig or Scaffod:",contigs_or_scaffolds,"-",round(contigs_or_scaffolds/total*100,2),"%",\
          "\nComplete or Chromosome:",complete_or_chromosome,"-",round(complete_or_chromosome/total*100,2),"%",\
         "\nComplete with Chromosomes and Plasmid(s):",complete_or_chromosome_and_plasmids,"-",\
          round(complete_or_chromosome_and_plasmids/total*100,2),"%")
    print("ATCC Total:",atcc_total,round(atcc_total/total*100,2),"%\nContig or Scaffod:",atcc_contigs_or_scaffolds,"-",round(atcc_contigs_or_scaffolds/atcc_total*100,2),"%",\
          "\nComplete or Chromosome:",atcc_complete_or_chromosome,"-",round(atcc_complete_or_chromosome/atcc_total*100,2),"%",\
         "\nComplete with Chromosomes and Plasmid(s):",atcc_complete_or_chromosome_and_plasmids,"-",\
          round(atcc_complete_or_chromosome_and_plasmids/atcc_total*100,2),"%")

    return refseq_accession_level, refseq_accession_line, genbank2refseq

#############################     Ensembl Data     #############################
def search_ensembl():
    ensembl_accessions = []
    with open('species_EnsemblBacteria.txt','r') as f:
        c = 0
        for line in f:
            if c < 2:
                c+=1
                pass
            else:
                for line in f:
                    line = line.lower()
                    line = line.split('\t')
                    ensembl_accessions.append(line[5])
    # print('ensembl accession list length:',len(ensembl_accessions))

    # Counts for whole dataset
    total = 0
    complete_or_chromosome = 0
    contigs_or_scaffolds = 0
    complete_or_chromosome_and_plasmids = 0
    # Counts for atcc strains in dataset
    atcc_total = 0
    atcc_complete_or_chromosome = 0
    atcc_contigs_or_scaffolds = 0
    atcc_complete_or_chromosome_and_plasmids = 0

    # Use for examining if any accessions are missing
    # prokaryotes_missing = 0
    # genbank_missing = 0
    # both_missing = 0
    for accession in ensembl_accessions:
        prokaryotes_line = ''
        genbank_line = ''
        if accession in prokaryotes_accession_line.keys():
            prokaryotes_line = prokaryotes_accession_line[accession]
        if accession in genbank_accession_line.keys():
            genbank_line = genbank_accession_line[accession]
        # if replicon information is available (this is true for all but ~800)
        if accession in prokaryotes_accession_replicon.keys():
            total+=1
            level = prokaryotes_accession_level[accession]
            if level == "complete":
                if "plasmid" in prokaryotes_accession_replicon[accession]:
                    if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                        atcc_complete_or_chromosome_and_plasmids+=1
                        atcc_total+=1
                    complete_or_chromosome_and_plasmids+=1
                else:
                    if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                        atcc_complete_or_chromosome+=1
                        atcc_total+=1
                    complete_or_chromosome+=1
            elif level == "chromosome":
                if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                    atcc_complete_or_chromosome+=1
                    atcc_total+=1
                complete_or_chromosome+=1
            elif level == "contig" or level == "scaffold":
                if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                    atcc_contigs_or_scaffolds+=1
                    atcc_total+=1
                contigs_or_scaffolds+=1
        elif accession in genbank_accession_level.keys():
            total+=1
            level = genbank_accession_level[accession]
            if level == "complete genome":
                if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                    atcc_complete_or_chromosome+=1
                    atcc_total+=1
                atcc_complete_or_chromosome+=1
            elif level == "chromosome":
                if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                    atcc_complete_or_chromosome+=1
                    atcc_total+=1
                complete_or_chromosome+=1
            elif level == "contig" or level == "scaffold":
                if "atcc" in prokaryotes_line or 'atcc' in genbank_line:
                    atcc_contigs_or_scaffolds+=1
                    atcc_total+=1
                contigs_or_scaffolds+=1

        ## This section keeps track of accessions that are not found in either the prokaryotes
        ## CSV or the assembly_summary_genbank.txt file. An NCBI search of some of these missing
        ## accessions suggests that they have all been "updated"

    #     if accession not in genbank_accession_level.keys() and accession not in prokaryotes_accession_replicon.keys():
    #         # print(accession)
    #         both_missing+=1
    #     if accession not in genbank_accession_level.keys():
    #         genbank_missing+=1
    #     if accession not in prokaryotes_accession_replicon.keys():
    #         prokaryotes_missing+=1
    # print("Missing from prokaryotes table:",prokaryotes_missing,\
    #      "\nMissing from genbank summary:",genbank_missing,\
    #      "\nMissing from both:",both_missing)


    print('\nEnsembl data:')
    print("\nTotal:",total,"\nContig or Scaffod:",contigs_or_scaffolds,"-",round(contigs_or_scaffolds/total*100,2),"%",\
          "\nComplete of Chromosome:",complete_or_chromosome,"-",round(complete_or_chromosome/total*100,2),"%",\
         "\nComplete with Chromosomes and Plasmid(s):",complete_or_chromosome_and_plasmids,"-",\
          round(complete_or_chromosome_and_plasmids/total*100,2),"%")
    print("ATCC Total:",atcc_total,round(atcc_total/total*100,2),"%\nContig or Scaffod:",atcc_contigs_or_scaffolds,"-",round(atcc_contigs_or_scaffolds/atcc_total*100,2),"%",\
          "\nComplete of Chromosome:",atcc_complete_or_chromosome,"-",round(atcc_complete_or_chromosome/atcc_total*100,2),"%",\
         "\nComplete with Chromosomes and Plasmid(s):",atcc_complete_or_chromosome_and_plasmids,"-",\
          round(atcc_complete_or_chromosome_and_plasmids/atcc_total*100,2),"%")

#############################       JGI Data       #############################
## Notes on process:
# > 0 chromo = chromosome level completion
# if > 0 chromos and > 0 plasmids then = chromos/plasmid level completion

def search_jgi():
    # Counts for whole dataset
    total = 0
    complete_or_chromosome = 0
    contigs_or_scaffolds = 0
    complete_or_chromosome_and_plasmids = 0
    # Counts for atcc strains in dataset
    atcc_total = 0
    atcc_complete_or_chromosome = 0
    atcc_contigs_or_scaffolds = 0
    atcc_complete_or_chromosome_and_plasmids = 0
    with open('jgi_download_16jun2021_simplified.txt','r',encoding='latin-1') as f:
        c=0
        for line in f:
            if c == 0:
                c+=1
            else:
                line = line.lower()
                line_arr = line.replace("\n","").split("\t")
                # Count chromosomes and plasmid
                if line_arr[4] == '':
                    chrom_count = 0
                else:
                    chrom_count = int(line_arr[4])
                if line_arr[5] == '':
                    plas_count = 0
                else:
                    plas_count = int(line_arr[5])
                # incremenet groups
                if chrom_count > 0 and plas_count > 0:
                    complete_or_chromosome_and_plasmids+=1
                    if "atcc" in line:
                        atcc_complete_or_chromosome_and_plasmids+=1
                        atcc_total +=1
                elif chrom_count > 0:
                    complete_or_chromosome +=1
                    if "atcc" in line:
                        atcc_complete_or_chromosome+=1
                        atcc_total+=1
                else:
                    contigs_or_scaffolds+=1
                    if "atcc" in line:
                        atcc_contigs_or_scaffolds+=1
                        atcc_total+=1
                total+=1
    print("\nJGI data")
    print("Total:",total,"\nContig or Scaffod:",contigs_or_scaffolds,"-",round(contigs_or_scaffolds/total*100,2),"%",\
          "\nComplete of Chromosome:",complete_or_chromosome,"-",round(complete_or_chromosome/total*100,2),"%",\
         "\nComplete with Chromosomes and Plasmid(s):",complete_or_chromosome_and_plasmids,"-",\
          round(complete_or_chromosome_and_plasmids/total*100,2),"%")
    print("ATCC Total:",atcc_total,round(atcc_total/total*100,2),"%\nContig or Scaffod:",atcc_contigs_or_scaffolds,"-",round(atcc_contigs_or_scaffolds/atcc_total*100,2),"%",\
          "\nComplete of Chromosome:",atcc_complete_or_chromosome,"-",round(atcc_complete_or_chromosome/atcc_total*100,2),"%",\
         "\nComplete with Chromosomes and Plasmid(s):",atcc_complete_or_chromosome_and_plasmids,"-",\
          round(atcc_complete_or_chromosome_and_plasmids/atcc_total*100,2),"%")


refseq_accession_level, refseq_accession_line, genbank2refseq = search_ncbi_refseq()
search_ensembl()
search_jgi()
