#! /usr/bin/python
import re
import string
# function to extract atcc or nctc id from refseq record
# regular expressions and ad hoc code were used to capture all the IDs 
def find_id(line):
    line = line.replace('\n','').lower().replace('atcc(b)','atcc ')
    if "atcc sd5219" in line:
        return 'atcc sd5219'
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
                # if neither atcc or nctc id could be extracted at this point, the following block
                # prints out information that could be used to update this function in order
                # to catch all examples
                else:
                    print('Not found')
                    print(id)
                    print(line)
                    return 'error'

# Go through NCBI RefSeq Bacteria assembly summary file and find all assemblies of ATCC products
def search_ncbi_refseq():
    out_lines = []
    # iterate through each line
    with open('assembly_summary_refseq.txt','r',encoding='utf-8') as f:
        c = 0
        # skip header lines
        for line in f:
            if line[0] == '#':
                pass
            else:
                # convert to lowercase and split on tab
                line = line.lower()
                line_arr = line.split('\t')
                # pull out GCF, GCA, assembly level, and url
                gcf = line_arr[0]
                gca = line_arr[17]
                level = line_arr[11]
                url =line_original[-4]
                # initialize product id (will become atcc or nctc id)
                product_id = ''
                if "atcc" in line or "nctc" in line:
                    product_id = find_id(line)
                # if a product id is found, append the data to the output lines
                if product_id == 'error':
                    pass
                elif product_id != '' :
                    out_lines.append([product_id,gcf,level,url])
    return out_lines
out_lines = search_ncbi_refseq()

# print the resulting lines to file

# with open('atcc_nctc_in_refseq_id-gcf-level-url.txt','w') as f:
#     for line in out_lines:
#         f.write('\t'.join(line))
#         f.write('\n')
