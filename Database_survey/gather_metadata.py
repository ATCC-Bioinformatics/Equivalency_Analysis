import requests
import pickle
# go through the file created from gather_atcc_nctc_assemblies.py
assemblies = []
for line in open('atcc_nctc_in_refseq_id-gcf-level-url.txt','r'):
        assemblies.append(line)
# create a dictionary to store all of the data
c = 1
metadata = {}
for line in assemblies:
    print(c,'of ',len(assemblies))
    line = line.replace('\n','').split('\t')
    # pull out url and generate path to the assembly_report.txt
    url = line[-1].replace('ftp://','https://')
    name = url.split('/')[-1]
    new_url = url + '/' + name + '_assembly_report.txt'
    # pull out basic metadata field for each assembly and store in dictionary
    data = {}
    data['Product ID'] = line[0]
    data['GCF'] = line[1]
    data['Assembly Level'] = line[2]
    # download the assembly_report
    r = requests.get(new_url)
    # go through assembly report, pull out each metadata field, and store in dictionary
    for line in r.text.split('\n'):
        # metadata fields start with a # and contain a : delimiter
        if line[:2] == "# " and ':' in line:
            line = line[2:].replace('\r','').split(':')
            data[line[0]] = line[1].strip()
    # store assembly metadata dictionary in main dictionary
    metadata[data['GCF']] = data
    c+=1

# pickle the metadata dictionary to be parsed in parse_metadata.py
with open('assembly_metadata.pkl', 'wb') as file:
    pickle.dump(metadata, file)
