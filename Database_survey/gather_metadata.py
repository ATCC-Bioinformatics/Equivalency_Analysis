import requests
import pickle
assemblies = []
for line in open('atcc_nctc_in_refseq_id-gcf-level-url.txt','r'):
        assemblies.append(line)
c = 1
metadata = {}
for line in assemblies:
    print(c,'of ',len(assemblies))
    line = line.replace('\n','').split('\t')
    url = line[-1].replace('ftp://','https://')
    name = url.split('/')[-1]
    new_url = url + '/' + name + '_assembly_report.txt'
    data = {}
    data['Product ID'] = line[0]
    data['GCF'] = line[1]
    data['Assembly Level'] = line[2]
    r = requests.get(new_url)
    for line in r.text.split('\n'):
        if line[:2] == "# " and ':' in line:
            line = line[2:].replace('\r','').split(':')
            data[line[0]] = line[1].strip()
    metadata[data['GCF']] = data
    c+=1


with open('assembly_metadata.pkl', 'wb') as file:
    pickle.dump(metadata, file)
