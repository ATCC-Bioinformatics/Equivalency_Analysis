#! /usr/bin/python
import re
import string
# NCTC 2 ATCC
nctc2atcc={}
for line in open('atcc2nctc.csv','r'):
    line = line.replace('\n','').split(',')
    nctc = line[0].replace('NCTC ','')
    atcc = line[-1].replace('ATCC ','')
    nctc2atcc[nctc]= atcc

##### find nctc id #######
def find_nctc(line):
    line = line.replace('\n','').lower().replace('atcc(b)','atcc ')
    if "atcc sd5219" in line:
        return ''
    else:
        id = re.findall(r'atcc\s*\d+|atcc\s*[=:_-]\s*\d+|atcc\s*[=:_-]*\s*baa\s*[=:_-]*\s*\d+|atcc\s*[=:_-]*\s*pta\s*[=:_-]*\s*\d+|atcc\s*[=:_-]*\s*vr\s*[=:_-]*\s*\d+',line)
        if len(id) > 0:
            return id[0]
        else:
            id = re.findall(r'nctc\s*\d+|nctc\s*[=:_-]\s*\d+',line.lower())
            if len(id) > 0:
                return id[0]
            else:
                if re.search('strain=phila',''.join(line)) or re.search('philadelphia_1_atcc',''.join(line)):
                    return '33152'
                else:
                    print('Not found')
                    print(id)
                    print(line)
                    return 'error'

#############################      NCBI Data       #############################
def search_ncbi_refseq():
    found_count = 0
    out_lines = []
    with open('assembly_summary_refseq.txt','r',encoding='utf-8') as f:
        c = 0
        for line in f:
            if line[0] == '#':
                pass
            else:
                line_original = line.split('\t')
                line = line.lower()
                line_arr = line.split('\t')
                gcf = line_arr[0]
                gca = line_arr[17]
                level = line_arr[11]
                url =line_original[-4]
                # If the accession number is one included in the Prokaryotic RefSeq
                product_id = ''
                # nctc_id = find_nctc([line,prokaryotes_line])
                if level == "complete genome":
                    if "atcc" in line or "nctc" in line:
                        # Increment assembly-type count
                        product_id = find_nctc(line)
                        found_count+=1
                elif level == "chromosome":
                    if "atcc" in line or "nctc" in line:
                        # Increment assembly-type count
                        product_id = find_nctc(line)
                        found_count+=1
                elif level == "contig" or level == "scaffold":
                    if "atcc" in line or "nctc" in line:
                        # Increment assembly-type count
                        product_id = find_nctc(line)
                        found_count+=1
                if product_id == 'error':
                    pass
                elif product_id != '' :#and 'vr' not in nctc_id:
                    out_lines.append([product_id,gcf,level,url])
    return out_lines,found_count

out_lines,found_count = search_ncbi_refseq()

# with open('atcc_nctc_in_refseq_id-gcf-level-url.txt','w') as f:
#     for line in out_lines:
#         f.write('\t'.join(line))
#         f.write('\n')
