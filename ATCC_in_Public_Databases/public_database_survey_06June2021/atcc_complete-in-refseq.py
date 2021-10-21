#! /usr/bin/python
import re
# The prokaryotes table contains replicon information not found in the assembly
# summary tables
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
            prokaryotes_accession_line[accession] = line

def search_ncbi_refseq():
    # Store data to search for ATCC keyword later
    genbank2refseq = {}

    atcc_lines = []
    atcc_lines_complete = []

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

                prokaryotes_line = ''
                if gca in prokaryotes_accession_line.keys():
                    prokaryotes_line = prokaryotes_accession_line[gca]

                # If the accession number is one included in the Prokaryotic RefSeq
                if level == "complete genome":
                    if "atcc" in line or 'atcc' in prokaryotes_line:
                        # Append to atcc lines and atcc complete
                        atcc_lines_complete.append([line,prokaryotes_line])
                        atcc_lines.append([line,prokaryotes_line])
                elif level == "chromosome":
                    if "atcc" in line or 'atcc' in prokaryotes_line:
                        # Append to atcc lines and atcc complete
                        atcc_lines_complete.append([line,prokaryotes_line])
                        atcc_lines.append([line,prokaryotes_line])
                elif level == "contig" or level == "scaffold":
                    if "atcc" in line or 'atcc' in prokaryotes_line:
                        # Append to atcc lines and atcc complete
                        atcc_lines.append([line,prokaryotes_line])


    return atcc_lines_complete, atcc_lines, genbank2refseq
atcc_lines_complete, atcc_lines, genbank2refseq = search_ncbi_refseq()



## This chunk of code is very ad hoc and customized to catch all of the atcc ids
atcc_ids = []
atcc_ids_lines = {}
total = 0
found = 0
for both_lines in atcc_lines_complete:
# for line in atcc_lines: # atcc_lines is the set of all atcc products in RefSeq
    line = both_lines[0].split('\t')
    atcc = [e for e in line if "atcc" in e]
    if len(atcc) > 1:
        atcc = [e for e in atcc if "strain" in e]
        try:
            atcc = atcc[0]
            atcc = re.findall(r'atcc\s*\d+|atcc\:*\s*\-*\w*\-*\d+|atcc\(b\) \d+',atcc)
            atcc_id = atcc[0]
            atcc_id = atcc_id.replace('atcc','')
    #         print(atcc_id)
            atcc_ids.append(atcc_id.replace(' ',''))
            found +=1
        except:
            atcc = [e for e in line if "atcc" in e][0]
            atcc_id = atcc_id.replace('atcc','')
#             print(atcc_id)
            atcc_ids.append(atcc_id.replace(' ',''))
            found += 1

    else:
        try:
            atcc = atcc[0]
            atcc = re.findall(r'atcc\s*\d+|atcc\:*\s*\-*\w*\-*\d+|atcc\(b\) \d+',atcc)
            atcc_id = atcc[0]
            atcc_id = atcc_id.replace('atcc','')
            atcc_ids.append(atcc_id.replace(' ',''))
#             print(atcc_id)
            found += 1
        except:
            if re.search('strain=phila',''.join(line)) or re.search('philadelphia_1_atcc',''.join(line)):
                atcc_id = '33152'
#                 print(atcc_id)
                atcc_ids.append(atcc_id.replace(' ',''))
                found += 1
            else:
                prok_line = both_lines[1].split(',')
                atcc = [e for e in prok_line if "atcc" in e]
                atcc = atcc[0]
                atcc = re.findall(r'atcc\s*\d+|atcc\:*\s*\-*\w*\-*\d+|atcc\(b\) \d+',atcc)
                atcc_id = atcc[0]
                atcc_id = atcc_id.replace('atcc','')
                atcc_id = atcc_id.replace(' ','')
                atcc_ids.append(atcc_id.replace(' ',''))
                found += 1

    if atcc_id.replace(' ','') in atcc_ids_lines.keys():
        atcc_ids_lines[atcc_id.replace(' ','')].append(line)
    else:
        atcc_ids_lines[atcc_id.replace(' ','')] = [line]
#     if "13880" in atcc_id:
#         print(line)
    total+=1
# print(found,'of',total)


## Check if product is published by examining genome status page
onecodex_published_list = []
for line in open("genome_status.csv","r"):
    line_arr = line.replace('\n','').split(',')
    if re.search('Unicycler Hybrid',line_arr[10]) and line_arr[8] == 'published':
        assemblyID = line_arr[5]
        onecodex_published_list.append(assemblyID)
# print(len(onecodex_published_list))
onecodex2atccid = {}
atccid2onecodex = {}
atccid2illumina_nanoporeids = {}
for line in open("sbc_dashboard.csv","r"):
    line_arr = line.replace('\n','').replace('\"','').split(',')
    onecodex2atccid[line_arr[0]] = line_arr[1]
    atccid2onecodex[line_arr[1]] = line_arr[0]
    atccid2illumina_nanoporeids[line_arr[1]] = [line_arr[2],line_arr[3]]
published_list = []  # contains all the IDs of the currently published ATCC prodcuts
for assemblyid in onecodex_published_list:
    if assemblyid in onecodex2atccid.keys():
        atcc_id = onecodex2atccid[assemblyid]
        published_list.append(atcc_id)
    else:
        pass
# print(len(published_list))


# ATCC product published on portal with assembly found in RefSeq
c=0
atcc_refseq_published = {}
for atcc_id in atcc_ids:
    if atcc_id.upper() in published_list:
        atcc_refseq_published[atcc_id.upper()] = list(atcc_ids_lines[atcc_id])
        c+=1
print(c)
print(c,' GCF assemblies linked to ',len(atcc_refseq_published), ' atcc published products')
atcc_refseqpublished_gcfs = {}
for key in atcc_refseq_published.keys():
    atcc_refseqpublished_gcfs[key] = []
    for line in atcc_refseq_published[key]:
        line_arr = line[0].split(",")
        gcf = [e for e in line_arr if "gcf" in e][0]
        gcf = gcf[gcf.index("gcf_"):]
        gcf = gcf[0:gcf.index(".")+2]
        atcc_refseqpublished_gcfs[key].append(gcf)
c=0
k=0
for key in atcc_refseqpublished_gcfs.keys():
    k+=1
    for gcf in atcc_refseqpublished_gcfs[key]:
        c+=1
print(c,k, 'should be same as above')

with open('atcc_complete-in-refseq.txt','w') as f:
    for key in atcc_refseqpublished_gcfs.keys():
        ilm = atccid2illumina_nanoporeids[key][0]
        ont = atccid2illumina_nanoporeids[key][1]
        for gcf in atcc_refseqpublished_gcfs[key]:
            f.write('\t'.join([key,ilm,ont,gcf]))
            f.write('\n')
