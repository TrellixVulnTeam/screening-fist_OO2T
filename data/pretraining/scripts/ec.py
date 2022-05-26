#!/usr/bin/env python
import sys
import os
import re
import requests
from tqdm import tqdm
import pymongo as mongo 

def mongodb():
    MONGODB_ADDR="mongodb://localhost:27017/"
    DBNAME = 'ec_nums'
    client = mongo.MongoClient(MONGODB_ADDR)
    db = client.get_database(DBNAME)
    return db


def get_kegg_ec(ec):
    url = f"http://rest.kegg.jp/get/ec:{ec}"
    r = requests.get(url)
    if r.status_code == 200:
        return r.text

def parse_kegg(entry):
    fields = re.findall('\n([A-Z]+)\s', entry)
    fmt = lambda s : s.replace('\n','').split()
    return {i:fmt(entry[re.search(i, entry).end():re.search(j, entry).start()]) \
            for i, j in zip(fields, fields[1:] + ['///'])}

def main_stdin(stdin):
    db = mongodb()
    collection = 'ec'
    ecs = [i for i in stdin.read().split('\n') if i != '']
    for ec in tqdm(ecs):
        data = get_kegg_ec(ec)
        if data is not None:
            doc = parse_kegg(data)
            doc['ec'] = ec
            if len(list(db[collection].find({'ENTRY':{'$regex':ec}}).limit(1))) == 0:
                db[collection].insert_one(doc)

def main(args):
    db = mongodb()
    collection = 'ec'
    for arg in args:
        with open(arg) as f:
            ecs = [i for i in f.read().split('\n') if i != '']
        for ec in tqdm(ecs):
            data = get_kegg_ec(ec)
            if data is not None:
                doc = parse_kegg(data)
                doc['ec'] = ec
                if 'SUBSTRATE' in doc:
                    substrates_kegg = re.findall('CPD:(C[0-9]+)', ''.join(doc['SUBSTRATE']))
                    doc['substrate_ids'] = substrates_kegg
                if 'PRODUCT' in doc:
                    products_kegg = re.findall('CPD:(C[0-9]+)', ''.join(doc['PRODUCT']))
                    doc['product_ids'] = products_kegg
                if len(list(db[collection].find({'EC':ec}).limit(1))) == 0:
                    db[collection].insert_one(doc)

if __name__ == '__main__':
    if not sys.stdin.isatty():
        main_stdin(sys.stdin)
    else:
        main(sys.argv[1:])
