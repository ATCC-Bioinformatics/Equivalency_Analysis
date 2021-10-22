import pickle

# parse atcc2nctc.csv and create conversion dictionary
nctc2atcc={}
for line in open('atcc2nctc.csv','r'):
    line = line.replace('\n','').split(',')
    nctc = line[0].replace('NCTC ','')
    atcc = line[-1].replace('ATCC ','')
    nctc2atcc[nctc]= atcc

# unpickle data from gather_metadata.py
with open('assembly_metadata.pkl','rb') as f:
    data = pickle.load(f)

# collect complete set of metadata fields from metadata dictionary
keys_set = set()
for k in data.keys():
    for sub_key in data[k].keys():
        keys_set.add(sub_key)
# open file for writing out metadata
with open('atcc_nctc_gcf-level-metadata.txt','w') as f:
    # header contains each metadata field
    f.write('\t'.join(list(keys_set)))
    f.write('\n')
    # add each metadata field value for each assembly to the line variable
    # if assembly doesn't have that field, then add "n/a"
    for assembly in data.keys():
        # line will hold the tab delimited metadata fields
        line = ''
        # keep is for each assembly, unless it has an nctc product id with no atcc id equivalent
        keep = True
        # go through each metadata field in the entire set of available metadata fields collected 
        # in lines 16 through 20
        for key in keys_set:
            # if the field is Product ID, convert nctc id to atcc id
            if key == 'Product ID':
                if 'atcc' in data[assembly][key]:
                    line += data[assembly][key] + '\t'
                elif 'nctc' in data[assembly][key]:
                    prev_id = data[assembly][key]
                    nctc = data[assembly][key].replace('nctc','').replace(' ','')
                    if nctc in nctc2atcc.keys():
                        atcc_id = nctc2atcc[nctc]
                        line += atcc_id +' (' + prev_id + ')' + '\t'
                    # keep is only false if there is no nctc to atcc conversion (see comment on line 30)
                    else:
                        keep = False
                else:
                    print('problem')
            # if the assembly_report contained the metadata field, then add the value to the line
            elif key in data[assembly].keys():
                line += data[assembly][key] + '\t'
            # if the assembly_report did not contain the metadata field, then add "n/a" to the line
            else:
                line += "n/a\t"
        # print line
        if keep == True:
            f.write(line)
            f.write('\n')
