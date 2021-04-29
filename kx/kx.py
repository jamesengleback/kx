import re
import os
import requests
import json
import pandas as pd
from tqdm import tqdm
from kx.ec_numbers import ec_numbers



# kx(1.1.1.1)
# ec -> sequence ids -> done/data
#   |               \-> to do
#   \-> compound ids -> smiles
#
# ec regex
# output dir
# concurrency
# restart / wait / daemon

def kx(ec, savedir=None):
    assert ec in ec_numbers # todo - grep

    # file setup
    if savedir is None:
        savedir = os.path.expanduser('~/.kx/')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    out_csv = os.path.join(savedir, f'{ec}-sequences.csv')
    ecjson = os.path.join(savedir,f'{ec}.json')

    # get those gene numbers & chems
    if not os.path.exists(ecjson):
        kegg_ex = get_ec(ec)
        substrates = re.findall('CPD:(C\d+)',''.join(kegg_ex['SUBSTRATE'])) # i know :/
        products = re.findall('CPD:(C\d+)',''.join(kegg_ex['PRODUCT']))
        genes = {re.findall('([A-Z]+):',i)[0]:re.findall('(\S+)',i.split(':')[1]) for i in kegg_ex['GENES']}
        d = {'ec':ec,
             'substrates':substrates,
             'products':products,
             'genes':genes}
        with open(ecjson,'w') as f:
            json.dump(d,f)
    else:
        with open(ecjson, 'r') as f:
            d = json.load(f)

    ## get sequences
    # track done
    # to csv
    pd.DataFrame([],columns = ['species','gene','sequence', 'structure']).to_csv(out_csv, index=False)
    todo = [[i,j] for i in d['genes'].keys() for j in d['genes'][i]]

    def done():
        df = pd.read_csv(out_csv, columns=['species','gene','sequence', 'structure'])
        return [[i,j] for i, j  in zip(df['species'], df['gene'])]

    # split into new fn
    def extract(species,gene):
        gene_data = get_gene(species, gene)
        if isinstance(gene_data, dict):
            aaseq = re.findall('[A-Z]+', ''.join(gene_data['AASEQ']))[0]
            if 'STRUCTURE' in gene_data.keys():
                struc = ','.join(re.findall('PDB:(.+)', ''.join(gene_data['STRUCTURE']))[0].split())
            else:
                struc = None
            data = pd.DataFrame([[species, gene, aaseq, struc]], 
                columns=['species','gene','sequence', 'structure'])
            return data

    for i in tqdm(todo):
        data = extract(*i)
        if data is not None:
            data.to_csv(out_csv, mode = 'a', header=None, index=False)
            todo.remove(i)


def get_gene(species, geneid):
    gene = re.sub("\([^)]*\)",'',geneid) # parenthes = bad
    r = requests.get(f'http://rest.kegg.jp/get/{species.lower()}:{gene}')
    if r.status_code == 200:
        data = r.text
        feilds = re.findall('\n([A-Z]+) ',data)
        idx = lambda feild : re.search(feild, data).span()
        chunks = {i:data[max(idx(i)):min(idx(j))].replace('  ','').splitlines() for i, j in zip(feilds, feilds[1:])}
        return chunks
    else:
        return r.status_code

def get_ec(ec_num):
    # just gets KEGG data as dict
    r = requests.get(f'http://rest.kegg.jp/get/ec:{ec_num}')
    if r.status_code == 200:
        data = r.text#.splitlines()
        feilds = re.findall('\n([A-Z]+) ',data)
        idx = lambda feild : re.search(feild, data).span()
        chunks = {i:data[max(idx(i)):min(idx(j))].replace('  ','').splitlines() for i, j in zip(feilds, feilds[1:])}
        return chunks
    else:
        return r.status_code
