import pickle

# NCTC 2 ATCC
nctc2atcc={}
for line in open('atcc2nctc.csv','r'):
    line = line.replace('\n','').split(',')
    nctc = line[0].replace('NCTC ','')
    atcc = line[-1].replace('ATCC ','')
    nctc2atcc[nctc]= atcc
# print(nctc2atcc)

with open('assembly_metadata.pkl','rb') as f:
    data = pickle.load(f)

keys_set = set()
for k in data.keys():
    for sub_key in data[k].keys():
        keys_set.add(sub_key)
# print(keys_set)

with open('atcc_nctc_gcf-level-metadata.txt','w') as f:
    f.write('\t'.join(list(keys_set)))
    f.write('\n')
    for assembly in data.keys():
        line = ''
        keep = True
        for key in keys_set:
            if key == 'Product ID':
                if 'atcc' in data[assembly][key]:
                    line += data[assembly][key] + '\t'
                elif 'nctc' in data[assembly][key]:
                    prev_id = data[assembly][key]
                    nctc = data[assembly][key].replace('nctc','').replace(' ','')
                    if nctc in nctc2atcc.keys():
                        atcc_id = nctc2atcc[nctc]
                        line += atcc_id +' (' + prev_id + ')' + '\t'
                    else:
                        keep = False
                else:
                    print('problem')
            elif key in data[assembly].keys():
                line += data[assembly][key] + '\t'
            else:
                line += 'n\\a\t'
        if keep == True:
            f.write(line)
            f.write('\n')
