#!/usr/bin/env python
import sys
import os
import re
import json
from pprint import pprint
import requests
from tqdm import tqdm
import pymongo as mongo 

def mongodb():
    MONGODB_ADDR="mongodb://localhost:27017/"
    DBNAME = 'ec_nums'
    client = mongo.MongoClient(MONGODB_ADDR)
    #assert DBNAME in client.list_database_names()
    db = client.get_database(DBNAME)
    return db


def get_kegg_cpd(cpd):
    url = f"http://rest.kegg.jp/get/cpd:{cpd}"
    r = requests.get(url)
    if r.status_code == 200:
        return r.text

def parse_kegg(entry):
    fields = re.findall('\n([A-Z]+)\s', entry)
    fmt = lambda s : s.replace('\n','').split()
    return {i:fmt(entry[re.search(i, entry).end():re.search(j, entry).start()]) \
            for i, j in zip(fields, fields[1:] + ['///'])}

def get_pubchem(pubchem_id):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/record/json'
    r = requests.get(url)
    if r.status_code == 200:
        return json.loads(r.text)['PC_Compounds']

def main():
    db = mongodb()
    collection = db['ec']
    
    ecs = collection.find({})
    cpds = []
    for doc in tqdm(ecs):
        substrates_kegg = doc['substrate_ids'] if 'substrate_ids' in doc else []
        products_kegg   = doc['product_ids']  if 'product_ids' in doc else []
        cpds += substrates_kegg + products_kegg
    cpds = list(set(cpds))

    collection = db['smiles']

    for cpd in tqdm(cpds):
        if len(list(collection.find({'cpd':cpd}).limit(1))) == 0:
            data = get_kegg_cpd(cpd)
            if data is not None:
                doc = parse_kegg(data)
                if 'DBLINKS' in doc:
                    dblinks_ = doc['DBLINKS']
                    dblinks = {i.replace(':',''):j for i,j in zip(dblinks_[::2], dblinks_[1::2])}
                    if 'PubChem' in dblinks:
                        pubchem_id = dblinks['PubChem']
                        pubchem_data = get_pubchem(pubchem_id) # list of dicts
                        if pubchem_data is not None:
                            for i in pubchem_data:
                                if i is not None:
                                    props = i['props']
                                    for j in props:
                                        if j['urn']['label'] == 'SMILES':
                                            smiles = j['value']['sval']
                                            doc = {'cpd':cpd, 'smiles':smiles}
                                            db['smiles'].insert_one(doc)


if __name__ == '__main__':
    main()
